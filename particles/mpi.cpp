#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include "SimpleIterator.h"
#include "Grid.h"
#include "Stats.h"

void add_stats(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
  Stats* in = (Stats *) invec;
  Stats* out = (Stats *) inoutvec;
  for (int i = 0; i < *len; i++) {
    out[i].aggregate_left(in[i]);
  }
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    const int n = read_int( argc, argv, "-n", 1000 );
    const bool fast = (find_option( argc, argv, "-no" ) != -1);
    const char *savename = read_string( argc, argv, "-o", NULL );
    const char *sumname = read_string( argc, argv, "-s", NULL );


    const double size = set_size( n );
    // We need to set the size of a grid square so that the average number of
    // particles per grid square is constant.  The simulation already ensures
    // that the average number of particles in an arbitrary region is constant
    // and proportional to the area.  So this is just a constant.
    const double grid_square_size = sqrt(0.0005) + 0.000001;
    const int num_grid_squares_per_side = size / grid_square_size;

    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


    particle_t *particles = (particle_t*) malloc(n*sizeof(particle_t));
    
    // Create an MPI type for particle_t.
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    // Create an MPI type and aggregation function for Stats.
    MPI_Datatype STATS;
    MPI_Datatype STATS_TYPES[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    int STATS_BLOCKLENS[3] = {1, 1, 1};
    MPI_Aint STATS_DISPLACEMENTS[3];
    STATS_DISPLACEMENTS[0] = offsetof(Stats, min);
    STATS_DISPLACEMENTS[1] = offsetof(Stats, avg);
    STATS_DISPLACEMENTS[2] = offsetof(Stats, n);
    MPI_Type_create_struct(3, STATS_BLOCKLENS, STATS_DISPLACEMENTS, STATS_TYPES, &STATS);
    MPI_Type_commit(&STATS);

    MPI_Op ADD_STATS;
    MPI_Op_create(add_stats, 1, &ADD_STATS);

    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 ) {
        particles = init_particles(n).release()->data();
    }
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    
    Stats local_stats;
    Stats global_stats;

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        // 
        //  collect all global data locally (not good idea to do)
        //
        MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if(!fast)
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
        for( int i = 0; i < nlocal; i++ )
        {
            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force(local[i], particles[j], local_stats);
        }

        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );
    }

    // Aggregate stats accumulated on local processes.
    MPI_Reduce(&local_stats, &global_stats, 1, STATS, ADD_STATS, 0, MPI_COMM_WORLD);

    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

      if(!fast)
      {
        //
        //  -the minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", global_stats.min, global_stats.avg);
        if (global_stats.min < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (global_stats.avg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      printf("\n");     
        
      //  
      // Printing summary data
      //  
      if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
