#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <memory>
#include <vector>
#include <stdio.h>
#include "Stats.h"

inline int min( int a, int b ) { return a < b ? a : b; }
inline double min( double a, double b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }
inline double max( double a, double b ) { return a > b ? a : b; }
inline int div_round_up(int dividend, int divisor) { return (dividend - 1 + divisor) / divisor; }

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
std::unique_ptr<std::vector<particle_t>> init_particles(int n);
void apply_force( particle_t &particle, particle_t &neighbor, Stats &stats);
void move( particle_t &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
