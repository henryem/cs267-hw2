#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "common.h"

#define NUM_THREADS 256

extern double size;
//
//  benchmarking program
//

struct GridMetadata {
  //FIXME
  double grid_square_size;
  int num_grid_squares_per_side;
  int num_grid_squares;

  GridMetadata() { } //FIXME

  GridMetadata(double grid_square_size_v, int num_grid_squares_per_side_v):
    grid_square_size(grid_square_size_v),
    num_grid_squares_per_side(num_grid_squares_per_side_v),
    num_grid_squares(num_grid_squares_per_side_v*num_grid_squares_per_side_v) { }
};

inline int flat_square_idx(const particle_t& p, const GridMetadata& g) {
  return 0; //FIXME
}

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

__global__ void compute_forces_gpu(particle_t* particles, int n) {
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particles[tid].ax = particles[tid].ay = 0;
  for(int j = 0 ; j < n ; j++) {
    apply_force_gpu(particles[tid], particles[j]);
  }
}

__global__ void move_gpu(particle_t * particles, int n, double size) {

  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t * p = &particles[tid];
  //
  //  slightly simplified Velocity Verlet integration
  //  conserves energy better than explicit Euler method
  //
  p->vx += p->ax * dt;
  p->vy += p->ay * dt;
  p->x  += p->vx * dt;
  p->y  += p->vy * dt;

  //
  //  bounce from walls
  //
  while(p->x < 0 || p->x > size) {
    p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
    p->vx = -(p->vx);
  }
  while(p->y < 0 || p->y > size) {
    p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
    p->vy = -(p->vy);
  }

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
        return 0;
    }
    
    const int n = read_int( argc, argv, "-n", 1000 );
    const bool fast = (find_option( argc, argv, "-no" ) != -1);
    const char *savename = read_string( argc, argv, "-o", NULL );
    const char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    const double size = set_size( n );

    // Particles are stored in an array of blocks of squares.  Each square
    // is large enough that particles can only move across 1 square per
    // simulated time step, but small enough that the expected number of
    // particles (and hopefully the maximum number) per square is a small
    // constant.
    // The block structure and squares are constant, but the square locations
    // of particles changes across iterations and is not fixed.  We use an
    // expandable Thrust vector for convenience.
    const double grid_square_size = cutoff + 0.0001;
    const int num_grid_squares_per_side = div_round_up_f(size, grid_square_size);
    const int num_grid_squares = num_grid_squares_per_side * num_grid_squares_per_side;
    //FIXME
    const GridMetadata grid(grid_square_size, num_grid_squares_per_side);
    thrust::host_vector<thrust::host_vector<particle_t> > ps(num_grid_squares);

    particle_t* particles = init_particles(n);
    // Generate the grid initially, serially, on the host.  This is slow, but
    // it cannot really be done efficiently on the GPU.  Note that the
    // particles are _copied_ into the grid so that there will be good
    // cache locality in the GPU code; the grid is the authoritative store
    // of the particles.
    for (unsigned int i = 0; i < n; i++) {
      particle_t& p = particles[i];
      ps[flat_square_idx(p, grid)].push_back(p);
    }

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
//    cudaMemcpy(d_particles, (*particles).data(), (*particles).size() * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ ) {
      // First, we must build the grid.



      //
      //  compute forces
      //

      int blks = div_round_up(n, NUM_THREADS);
      compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n);

      //
      //  move particles
      //
      move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);

      //
      //  save if necessary
      //
      if( fsave && (step%SAVEFREQ) == 0 ) {
        // Copy the particles back to the CPU
        cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
        save( fsave, n, particles);
      }
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    cudaFree(d_particles);
    free(particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
