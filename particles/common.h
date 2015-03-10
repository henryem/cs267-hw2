#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#ifndef __CUDACC__ // CUDA doesn't support things in <vector> or <memory>.
  #include <memory>
  #include <vector>
#endif
#include <stdio.h>
#include "Stats.h"

#ifndef __CUDACC__ // Already defined in CUDA code.
  inline int min( int a, int b ) { return a < b ? a : b; }
  inline int max( int a, int b ) { return a > b ? a : b; }
  inline double min( double a, double b ) { return a < b ? a : b; }
  inline double max( double a, double b ) { return a > b ? a : b; }
#endif
inline int div_round_up(int dividend, int divisor) { return (dividend - 1 + divisor) / divisor; }
inline int div_round_up_f(double dividend, double divisor) { return (int) ((dividend - 1.0 + divisor) / divisor); }

//
//  tuned constants
//
const double density = 0.0005;
const double mass = 0.01;
const double cutoff = 0.01;
const double min_r = (cutoff/100);
const double dt = 0.0005;

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
struct particle_t {
public:
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;

  particle_t():
      x(0.0), y(0.0), vx(0.0), vy(0.0), ax(0.0), ay(0.0) { }
  particle_t(double x_v, double y_v, double vx_v, double vy_v, double ax_v, double ay_v):
      x(x_v), y(y_v), vx(vx_v), vy(vy_v), ax(ax_v), ay(ay_v) { }
  particle_t(double x_v, double y_v):
      x(x_v), y(y_v), vx(0.0), vy(0.0), ax(0.0), ay(0.0) { }
};

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
double set_size(int n);
#ifndef __CUDACC__ // CUDA doesn't support things in <vector>.
  std::unique_ptr<std::vector<particle_t> > init_particles(int n);
#else
  particle_t* init_particles(int n);
#endif
void apply_force(particle_t &particle, particle_t &neighbor, Stats &stats);
void move( particle_t &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//

void dump_particle(particle_t *p, int count, int rank);
int find_boundry_proc(particle_t p, int reduced_n_proc, int * partition_grids, double grid_size, int * find_boundry_proc);

int find_proc_no(particle_t p, int reduced_n_proc, int * partition_grids, double grid_size);
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
