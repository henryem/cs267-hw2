#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/for_each.h>
#include <thrust/scatter.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include "common.h"

#define NUM_THREADS 256

struct GridMetadata {
  int num_particles;
  // Approximately, but not exactly, equal to square_size*side_count.
  double side_size;
  double square_size;
  int side_count;
  int count;

  GridMetadata(int num_particles_v, double side_size_v, double square_size_v, int side_count_v):
    num_particles(num_particles_v),
    side_size(side_size_v),
    square_size(square_size_v),
    side_count(side_count_v),
    count(side_count_v*side_count_v) { }

  __host__ __device__ int particle_to_flat_square_idx(const particle_t& p) const {
    int x_idx = (int) (p.x / square_size);
    int y_idx = (int) (p.y / square_size);
    return square_to_flat_square_idx(x_idx, y_idx);
  }

  __host__ __device__ int square_to_flat_square_idx(int square_x, int square_y) const {
    const int idx = square_x + side_count*square_y;
    return idx;
  }
};

// Functor to map a particle to its grid location.
struct GridSquareCmIndex : public thrust::unary_function<particle_t&, int> {
  const GridMetadata g;

  GridSquareCmIndex(const GridMetadata g_v) : g(g_v) { }

  __host__ __device__ int operator()(particle_t& p) const {
    return g.particle_to_flat_square_idx(p);
  }
};

// Functor to move a particle.
struct MoveParticle : public thrust::unary_function<particle_t&, void> {
  const GridMetadata g;

  MoveParticle(const GridMetadata g_v) : g(g_v) { }

  __host__ __device__ void operator()(particle_t& p) {
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while(p.x < 0 || p.x > g.side_size) {
      p.x  = p.x < 0 ? -(p.x) : 2*g.side_size-p.x;
      p.vx = -(p.vx);
    }
    while(p.y < 0 || p.y > g.side_size) {
      p.y  = p.y < 0 ? -(p.y) : 2*g.side_size-p.y;
      p.vy = -(p.vy);
    }
  }
};

struct IgnoreZeroPredicate : public thrust::unary_function<int, bool> {
  __host__ __device__ bool operator()(int i) {
    return i != 0;
  }
};

__device__ void apply_force_gpu(particle_t& particle, particle_t& neighbor) {
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__device__ int num_particles_in_square(int square_idx, const int* grid_offsets, const GridMetadata grid) {
  int square_start = grid_offsets[square_idx];
  // For all but the last grid square, the number of particles in the square
  // equals the next square's offset minus its own offset.  The last square
  // just has all the remaining particles.
  int square_end = (square_idx < (grid.count - 1)) ? grid_offsets[square_idx+1] : grid.num_particles;
  return square_end - square_start;
}

__global__ void compute_forces_gpu (particle_t* particles, const int* grid_offsets, const GridMetadata grid) {
  // Get grid square ID.  Each call to this function computes forces
  // for all of the particles in one grid square.
  const int square_idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (square_idx >= grid.count) return;

  const int square_x = square_idx % grid.side_count;
  const int square_y = square_idx / grid.side_count;
  const int first_particle_idx = grid_offsets[square_idx];
  const int num_ps = num_particles_in_square(square_idx, grid_offsets, grid);
  for (int particle_idx = first_particle_idx; particle_idx < first_particle_idx + num_ps; particle_idx++) {
    particles[particle_idx].ax = particles[particle_idx].ay = 0;
    // Iterate only over neighboring grid squares.  Note that we could reorder
    // this to access each neighbor only once, which would probably be better.
    for (int y_offset = -1; y_offset <= 1; y_offset++) {
      const int neighbor_y = square_y + y_offset;
      if (neighbor_y < 0 || neighbor_y >= grid.side_count) {
        continue;
      }
      for (int x_offset = -1; x_offset <= 1; x_offset++) {
        const int neighbor_x = square_x + x_offset;
        if (neighbor_x < 0 || neighbor_x >= grid.side_count) {
          continue;
        }
        const int neighbor_square_idx = grid.square_to_flat_square_idx(neighbor_x, neighbor_y);
        // Now we iterate over all the particles in the neighbor and apply
        // forces to our particle.
        const int first_neighbor_particle_idx = grid_offsets[neighbor_square_idx];
        const int num_neighbor_ps = num_particles_in_square(neighbor_square_idx, grid_offsets, grid);
        for (int neighbor_particle_idx = first_neighbor_particle_idx;
            neighbor_particle_idx < first_neighbor_particle_idx + num_neighbor_ps;
            neighbor_particle_idx++) {
          apply_force_gpu(particles[particle_idx], particles[neighbor_particle_idx]);
        }
      }
    }
  }
}

/*
 * @param particle_square_idx_storage is a preallocated device array of size
 *   at least grid.num_particles.  Its value can be arbitrary.  It will be
 *   clobbered.
 * @param grid_offsets will be populated with offsets corresponding to grid
 *  squares.  [grid_offsets[i], grid_offsets[i+1]) (note the inclusiveness!)
 *  is the set of indices in particles for grid square i, after this function
 *  returns.
 * @param grid_idx_storage is a preallocated device array of size at least
 *   grid.count.  Its value can be arbitrary.  It will be clobbered.
 * @param grid_count_storage is a preallocated device array of size at least
 *   grid.count.  Its value can be arbitrary.  It will be clobbered.
 */
void sort_to_bins(
    thrust::device_vector<particle_t>& particles,
    thrust::device_vector<int>& particle_square_idx_storage,
    thrust::device_vector<int>& grid_offsets,
    thrust::device_vector<int>& grid_idx_storage,
    thrust::device_vector<int>& grid_count_storage,
    const GridMetadata& grid) {
  thrust::fill(particle_square_idx_storage.begin(), particle_square_idx_storage.end(), 0);
  thrust::transform(
      particles.begin(),
      particles.end(),
      particle_square_idx_storage.begin(),
      GridSquareCmIndex(grid));

  // Sort the particles by column-major order in the grid.  Thrust offers no
  // sort function that leaves the keys in place, so we must actually allocate
  // memory for the square indices.  Another option is to use a comparison
  // functor; that would avoid memory allocation but end up invoking a less
  // efficient sorting algorithm that works on user-defined functors (rather
  // than the fast int-sorter).  It could be worth trying that, too.
  //
  // I think we could do much better than sort_by_key here.  This makes the
  // algorithm O(n log n) when we can do O(n).  Here is a sketch of a better
  // algorithm that seems to equal the asymptotic efficiency of a graph-based
  // or message-passing method, while having very low constant factors when
  // implemented on the GPU:
  //  * Associate with each particle its grid square coordinates (x,y).
  //  * Use a LSH algorithm to hash the grid square coordinates.  LSH means
  //    that the absolute difference between the hashes of two points is (with
  //    high probability) close to the L2 distance of the two points.
  //  * Sort the particles according to the hash.  The sorting algorithm we use
  //    is important:
  //    - Notice that, since particles may move
  //      only a bounded L2 distance in an iteration, there is a bound on the
  //      distance any particle needs to travel in this sorting step.  That is,
  //      the particles are already "almost sorted."
  //    - Therefore we want a sorting algorithm that can be implemented quickly
  //      on the GPU and that takes only O(k N) time, where k is the maximum
  //      (or average) distance moved by any element.
  //    - There are fast algorithms for nearly-sorted lists, but I haven't
  //      found on that is parallelizable yet.  That seems very possible,
  //      though.  So I'm leaving this for future work, due to time
  //      constraints.

  // Unrelated note: Thrust's sort_by_key seems to have a bug (or exercises a
  // memory leak in my code?) that kicks in for particle counts between 18,000
  // and 26,000 and causes this line to crash with a bad_alloc error.  For
  // larger or smaller particle counts, everything seems fine, which is
  // terrifying.
  thrust::sort_by_key(
      particle_square_idx_storage.begin(),
      particle_square_idx_storage.end(),
      particles.begin());

  // Compute the starting offset of each grid square (that is, the index of the
  // first particle in particles contained in each square).
  // I cannot find a totally natural way to do this using Thrust, so we have to
  // jump through some hoops.  We reduce_by_key to compute the number of
  // particles in each nonempty grid square and the corresponding square
  // indices.  Then we scatter into a vector that is initially filled with 0s,
  // mapping each count according to its index.  Then we scan with plus, so
  // that the value at the ith location is the number of particles in grid
  // squares preceding the ith grid square.
  thrust::fill(grid_idx_storage.begin(), grid_count_storage.end(), 0);
  thrust::fill(grid_count_storage.begin(), grid_count_storage.end(), 0);
  thrust::reduce_by_key(
      particle_square_idx_storage.begin(),
      particle_square_idx_storage.end(),
      thrust::make_constant_iterator(1),
      grid_idx_storage.begin(),
      grid_count_storage.begin());
  // Here we set grid_offsets[grid_idx_storage[i]] = grid_count_storage[i]
  // for all i.  Since some grid squares might have been empty, there might
  // be several trailing elements of grid_idx_storage that are zero.  So we
  // use scatter_if() to ignore those.
  thrust::fill(grid_offsets.begin(), grid_offsets.end(), 0);
  thrust::scatter_if(
      grid_count_storage.begin(),
      grid_count_storage.end(),
      grid_idx_storage.begin(),
      grid_count_storage.begin(),
      grid_offsets.begin(),
      IgnoreZeroPredicate());
  thrust::exclusive_scan(
      grid_offsets.begin(),
      grid_offsets.end(),
      grid_offsets.begin());
}

void simulate_forces(thrust::device_vector<particle_t>& particles, const thrust::device_vector<int>& grid_offsets, const GridMetadata& grid) {
  // The communication pattern is not simple, so we have to resort to writing
  // device code ourselves here.
  particle_t* d_particles = thrust::raw_pointer_cast(particles.data());
  const int* d_grid_offsets = thrust::raw_pointer_cast(grid_offsets.data());
  int num_blocks = div_round_up(grid.count, NUM_THREADS);
  compute_forces_gpu<<<num_blocks, NUM_THREADS>>>(d_particles, d_grid_offsets, grid);
}

void simulate_movement(thrust::device_vector<particle_t>& particles, const GridMetadata& grid) {
  thrust::for_each(particles.begin(), particles.end(), MoveParticle(grid));
}



int main( int argc, char **argv )
{    
  // This takes a few seconds to initialize the runtime
  cudaThreadSynchronize();

  if (find_option( argc, argv, "-h" ) >= 0) {
    printf( "Options:\n" );
    printf( "-h to see this help\n" );
    printf( "-n <int> to set the number of particles\n" );
    printf( "-o <filename> to specify the output file name\n" );
    printf( "-no turns off all correctness checks and particle output\n");
    return 0;
  }

  const int n = read_int( argc, argv, "-n", 1000 );
  const bool fast = (find_option( argc, argv, "-no" ) != -1);
  const char *savename = read_string( argc, argv, "-o", NULL );
  const char *sumname = read_string( argc, argv, "-s", NULL );

  FILE *fsave = ((!fast) && savename) ? fopen( savename, "w" ) : NULL;
  const double size = set_size( n );

  // Particles are stored in a flattened array of squares.  Each square
  // is large enough that particles can only move across 1 square per
  // simulated time step, but small enough that the expected number of
  // particles (and hopefully the maximum number) per square is a small
  // constant.
  // Following NVIDIA's example particle simulator, we store the particles
  // in a flattened array and sort them when we need to change the grid
  // structure.  That is, we use one vector
  // of size n to store the particles in sorted order (column major by grid
  // square, with arbitrary order within each square), and a second vector
  // to store the starting offsets for the
  // particles contained in each square.  Then rebuilding the grid involves
  // an in-place sort and recomputing the index.
  const double square_size = cutoff + 0.0001;
  const int side_count = div_round_up_f(size, square_size);
  const GridMetadata grid(n, size, square_size, side_count);

  double init_on_host_time = read_timer();
  //TODO: This part of the initialization is serial.  It is embarrassingly
  // parallel and could easily be done on the GPU.
  particle_t* particles = init_particles(n);
  thrust::host_vector<particle_t> ps(particles, particles+n);
  init_on_host_time = read_timer() - init_on_host_time;

  double init_on_device_time = read_timer();
  // Copy the particles to the GPU.
  thrust::device_vector<particle_t> d_ps = ps;
// Allocate the structure that maps grid locations to offsets in
  // d_ps.  Like d_ps, this will be populated inside the simulation loop,
  // and for now is uninitialized.
  thrust::device_vector<int> d_grid_offsets(grid.count);
  // Allocate scratch space that the algorithm will need.
  thrust::device_vector<int> d_particle_square_idx_storage(grid.num_particles);
  thrust::device_vector<int> d_grid_idx_storage(grid.count);
  thrust::device_vector<int> d_grid_count_storage(grid.count);
  init_on_device_time = read_timer() - init_on_device_time;

  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer();

  for (int step = 0; step < NSTEPS; step++) {
    // First, we must build the grid.
    sort_to_bins(d_ps, d_particle_square_idx_storage, d_grid_offsets, d_grid_idx_storage, d_grid_count_storage, grid);

    // Now we can simulate forces and movement.
    simulate_forces(d_ps, d_grid_offsets, grid);
    simulate_movement(d_ps, grid);

    if( fsave && (step%SAVEFREQ) == 0 ) {
      // Copy the particles back to the CPU.
      ps = d_ps;
      save( fsave, n, ps.data());
    }
  }
  cudaThreadSynchronize();
  simulation_time = read_timer( ) - simulation_time;

  printf( "CPU-GPU copy time = %g seconds\n", init_on_device_time);
  printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );

  free(particles);
  if( fsave )
    fclose( fsave );

  return 0;
}
