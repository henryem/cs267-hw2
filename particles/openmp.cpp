#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include "SimpleIterator.h"
#include "Grid.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int numthreads;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        printf( "-p <int> to set the (maximum) number of threads used\n");
        return 0;
    }

    const int n = read_int( argc, argv, "-n", 1000 );
    const bool fast = (find_option( argc, argv, "-no" ) != -1);
    const char *savename = read_string( argc, argv, "-o", NULL );
    const char *sumname = read_string( argc, argv, "-s", NULL );
    const int num_threads_override = read_int( argc, argv, "-p", 0);

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      
    const double size = set_size( n );
    // We need to set the size of a grid square so that the average number of
    // particles per grid square is constant.  The simulation already ensures
    // that the average number of particles in an arbitrary region is constant
    // and proportional to the area.  So this is just a constant.
    const double grid_square_size = sqrt(0.0005) + 0.000001;
    const int num_grid_squares_per_side = size / grid_square_size;
    printf("Using %d grid squares of side-length %f for %d particles.\n", num_grid_squares_per_side*num_grid_squares_per_side, grid_square_size, n);
    std::unique_ptr<std::vector<particle_t> > particles = init_particles(n);

    if (num_threads_override > 0) {
      omp_set_dynamic(0);
      omp_set_num_threads(num_threads_override);
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    int max_num_threads = omp_get_max_threads();

    // User-defined reductions aren't available in the version of OMP we're
    // using.  Instead, we accumulate per-thread stats in this global array
    // and reduce manually when we're done.
    Stats per_thread_stats[max_num_threads];

    // Shared across threads.
#ifdef PARALLEL_GRID
    std::unique_ptr<Grid> g0(new Grid(size, num_grid_squares_per_side));
    std::unique_ptr<Grid> g1(new Grid(size, num_grid_squares_per_side));
#else
    std::unique_ptr<Grid> g(new Grid(size, num_grid_squares_per_side));
#endif

    #pragma omp parallel
    {
    numthreads = omp_get_num_threads();
    for (int step = 0; step < 1000; step++) {
      //TODO: Does this need to be declared private?
      int thread_idx;

      //FIXME: If this is the first step, initialize the grid without
      // respecting cache locality here.  Since we cannot use the existing
      // grid, we have to just divide the particles arbitrarily.  This
      // means that the subsequent code for simulating forces and movement
      // will have almost no cache locality on the first iteration: Each thread
      // has picked up an arbitrary subset of the particles to insert into the
      // grid, and then the threads are responsible for simulating a different,
      // mostly-disjoint subset of the particles.  On subsequent iterations,
      // only the particles that have moved will cause cache misses, so we
      // should have much better locality.  If we want to really optimize,
      // it may be worth rethinking how we store particles and communicate among
      // threads.  But at that point we might as well write distributed-memory
      // code.

      // Here we are building the grid that maps locations to sets of
      // particles.  This step does O(n) work, so it is a bottleneck if done
      // serially.  For performance comparisons, we have two versions of the
      // grid-formation code.  The second simply forms the grid serially, in a
      // single arbitrary thread.  The first is parallel and attempts
      // some cache locality.  Each thread is responsible for re-inserting
      // the grid elements that previously lay in its subgrid.
      //NOTE: We could instead re-insert each particle right after moving it.
      // This would be faster, but it would require us to think about
      // simultaneous parallel delete and add, while the current scheme needs
      // only support parallel add.  (Deleting the entire grid at once is an
      // O(1) operation, so we can do it in one thread with a barrier.)
#ifdef PARALLEL_GRID
      //FIXME: Explain.  We need the old Grid around to find the particles in
      // it.
      #pragma omp single
      if (step % 2 == 0) {
        g1.reset(new Grid(size, num_grid_squares_per_side));
      } else {
        g0.reset(new Grid(size, num_grid_squares_per_side));
      }

      thread_idx = omp_get_thread_num();
      //FIXME: Breaks abstraction.  The Grid should divide itself up and give
      // us an iterator over particle_t& for each thread.
      std::unique_ptr<SimpleIterator<int> > grid_squares_for_thread = get_squares_for_thread(thread_idx);

      // Now insert each guy into the new grid
#else
      #pragma omp single
      g.reset(new Grid(size, num_grid_squares_per_side, *particles));
#endif

      //TODO: Could improve data locality by blocking according to the block
      // structure of the grid.  That would require keeping track, dynamically,
      // of the locations of each particle.  It would be interesting to test
      // whether manually allocating sub-blocks (as in the distributed memory
      // code) to threads improves things further.
      #pragma omp for
      for (int i = 0; i < n; i++) {
        thread_idx = omp_get_thread_num();
        particle_t& p = (*particles)[i];
        p.ax = p.ay = 0;
        std::unique_ptr<SimpleIterator<particle_t&> > neighbors = (*g).neighbor_iterator(p);
        while (neighbors->hasNext()) {
          particle_t& neighbor = neighbors->next();
          apply_force(p, neighbor, per_thread_stats[thread_idx]);
        }
      }

      // There is an implicit barrier here, which is important for correctness.
      // (Technically, some asynchrony could be allowed: A thread's sub-block
      // can be moved once it receives force messages from its neighboring
      // sub-blocks.)

      //
      //  move particles
      //
      #pragma omp for
      for (int i = 0; i < n; i++) {
        move((*particles)[i]);
      }

      if (!fast) {
        //
        //  save if necessary
        //
        #pragma omp master
        if( fsave && (step%SAVEFREQ) == 0 ) {
          save( fsave, n, (*particles).data() );
        }
      }
    }
    }
    simulation_time = read_timer( ) - simulation_time;

    // Could do a tree reduce here, but it seems unnecessary.
    Stats overall_stats;
    for (int thread_idx = 0; thread_idx < max_num_threads; thread_idx++) {
      overall_stats.aggregate_left(per_thread_stats[thread_idx]);
    }

    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if (!fast) {
      //
      //  -the minimum distance absmin between 2 particles during the run of the simulation
      //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
      //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
      //
      //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
      //
      printf( ", absmin = %lf, absavg = %lf", overall_stats.min, overall_stats.avg);
      if (overall_stats.min < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
      if (overall_stats.avg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    if( fsave )
        fclose( fsave );
    
    return 0;
}
