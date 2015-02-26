#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "Grid.h"
#include "SimpleIterator.h"

void simulate_step(int num_particles, particle_t* particles, double size, int num_grid_squares_per_side, Stats &s, FILE *fsave, bool fast, bool save_output) {
  Stats step_stats;
  Grid g(size, num_grid_squares_per_side, num_particles, particles);
  //
  //  compute forces
  //
  for( int i = 0; i < num_particles; i++ ) {
    particle_t& p = particles[i];
    p.ax = p.ay = 0;
    std::unique_ptr<SimpleIterator<particle_t&> > neighbors = g.neighbor_iterator(p);
    while (neighbors->hasNext()) {
      particle_t& neighbor = neighbors->next();
      apply_force(p, neighbor, step_stats);
    }
  }

  //
  //  move particles
  //
  for( int i = 0; i < num_particles; i++ ) {
    move( particles[i] );
  }

  if (!fast) {
    s.aggregate_left(step_stats);
    if (save_output) {
      save(fsave, num_particles, particles);
    }
  }
}

//
//  benchmarking program
//
int main( int argc, char **argv ) {
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
  
  int n = read_int( argc, argv, "-n", 1000 );
  bool fast = (find_option( argc, argv, "-no" ) == -1);

  char *savename = read_string( argc, argv, "-o", NULL );
  char *sumname = read_string( argc, argv, "-s", NULL );
  
  FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
  FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

  Stats overall_stats;
  particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
  const double size = set_size( n );
  // We need to set the size of a grid square so that the average number of
  // particles per grid square is constant.  The simulation already ensures
  // that the average number of particles in an arbitrary region is constant
  // and proportional to the area.  So this is just a constant.
  const double grid_square_size = sqrt(0.0005) + 0.000001;
  const int num_grid_squares_per_side = size / grid_square_size;
  init_particles( n, particles );
  
  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer( );

  for( int step = 0; step < NSTEPS; step++ ) {
    simulate_step(n, particles, size, num_grid_squares_per_side, overall_stats, fsave, fast, fsave && (step%SAVEFREQ) == 0);
  }
  simulation_time = read_timer( ) - simulation_time;
  
  printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
      fprintf(fsum,"%d %g\n",n,simulation_time);

  //
  // Clearing space
  //
  if( fsum ) {
    fclose(fsum);
  }
  free(particles);
  if (fsave) {
    fclose(fsave);
  }
  
  return 0;
}
