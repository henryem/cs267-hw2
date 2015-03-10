#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include "Stats.h"
#include <limits.h>

double size;

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
double set_size( int n )
{
    size = sqrt( density * n );
    return size;
}

//
//  Initialize the particle positions and velocities
//
#ifndef __CUDACC__
  std::unique_ptr<std::vector<particle_t>> init_particles(int n) {
    std::vector<particle_t>* ps = new std::vector<particle_t>(n);

    srand48( time( NULL ) );

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
      shuffle[i] = i;

    for( int i = 0; i < n; i++ )
    {
      particle_t& p = (*ps)[i];
      //
      //  make sure particles are not spatially sorted
      //
      int j = lrand48()%(n-i);
      int k = shuffle[j];
      shuffle[j] = shuffle[n-i-1];

      //
      //  distribute particles evenly to ensure proper spacing
      //
      p.x = size*(1.+(k%sx))/(1+sx);
      p.y = size*(1.+(k/sx))/(1+sy);

      //
      //  assign random velocities within a bound
      //
      p.vx = drand48()*2-1;
      p.vy = drand48()*2-1;
    }
    free( shuffle );

    return std::unique_ptr<std::vector<particle_t> >(ps);
  }
#else
  particle_t* init_particles(int n) {
    particle_t* p = (particle_t*) malloc(n*sizeof(particle_t));

    srand48( time( NULL ) );

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
      shuffle[i] = i;

    for( int i = 0; i < n; i++ )
    {
      //
      //  make sure particles are not spatially sorted
      //
      int j = lrand48()%(n-i);
      int k = shuffle[j];
      shuffle[j] = shuffle[n-i-1];

      //
      //  distribute particles evenly to ensure proper spacing
      //
      p[i].x = size*(1.+(k%sx))/(1+sx);
      p[i].y = size*(1.+(k/sx))/(1+sy);

      //
      //  assign random velocities within a bound
      //
      p[i].vx = drand48()*2-1;
      p[i].vy = drand48()*2-1;
    }
    free( shuffle );
    return p;
  }
#endif

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor, Stats &s)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
	  if (r2 != 0) {
      s.add_left(sqrt(r2) / cutoff);
    }
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
	
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
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
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}


int find_boundry_proc(particle_t p, int reduced_n_proc, int * partition_grids, double grid_size, int * boundry_proc)
{
    int grid_x = p.x / grid_size;
    int grid_y = p.y / grid_size;
    //printf("the grid_x is %d, the grid_y is %d\n",grid_x, grid_y);  


    int count = 0;
    int row = -1;
    int column = -1;
    int proc_per_side = sqrt(reduced_n_proc);

    // Find neiborhood proc
    for (int i = 1; partition_grids[i]< INT_MAX; i++){
        if (grid_x == partition_grids[i]-1){
            row = grid_y/partition_grids[1];
            boundry_proc[count] = i + row *proc_per_side; 
            count++;
        }
        else if (grid_x == partition_grids[i]){
            row = grid_y/partition_grids[1];
            boundry_proc[count] = (i-1) + row *proc_per_side; 
            count++;
        }
        if (grid_y == partition_grids[i]-1){
            column = grid_x/partition_grids[1];
            boundry_proc[count] = i *proc_per_side+ column ; 
            count++;
        }
        else if (grid_y == partition_grids[i]){
            column = grid_x/partition_grids[1];
            boundry_proc[count] = (i-1) *proc_per_side + column; 
            count++;
        }
    }
    // Find Digno proc only when there are already two proc
    if (count == 2){
        if ( (grid_x % partition_grids[1] == partition_grids[1] -1) && (grid_y % partition_grids[1] == partition_grids[1] -1 ))
        {
            boundry_proc[count] = (grid_x/partition_grids[1]+1) + (grid_y/partition_grids[1]+1)*proc_per_side;
            count++;
        }
        else if ( (grid_x % partition_grids[1] == 0) && (grid_y % partition_grids[1] == partition_grids[1] -1))
        {
            boundry_proc[count] = (grid_x/partition_grids[1]-1) + (grid_y/partition_grids[1]+1)*proc_per_side;
            count++;
        }
        else if ( (grid_x % partition_grids[1] == partition_grids[1] -1) && (grid_y % partition_grids[1] == 0))
        {
            boundry_proc[count] = (grid_x/partition_grids[1]+1) + (grid_y/partition_grids[1]-1)*proc_per_side;
            count++;
        }
        else if ( (grid_x % partition_grids[1] == 0) && (grid_y % partition_grids[1] == 0 ))
        {
            boundry_proc[count] = (grid_x/partition_grids[1]-1) + (grid_y/partition_grids[1]-1)*proc_per_side;
            count++;
        }
    }
    /*
    for (int i =0; i<count; i++){
        printf("The partition is at boundry of processor: %d\n", boundry_proc[i]);
    }
    */
    return count;

}

// Using for debugging
void dump_particle(particle_t *p, int count, int rank)
{
    for (int i = 0; i < count; i++){
        printf("p%d: (%f,%f)\n",rank, p[i].x, p[i].y);
    }
}

//
// given a particle, return the processor number it has been assigned.
//
int find_proc_no(particle_t p, int reduced_n_proc, int * partition_grids, double grid_size){
    
    //printf("the grid_size is %g\n", grid_size);
    // printf("the x is %g, the y is %g\n",p.x, p.y);    
    int grid_x = p.x / grid_size;
    int grid_y = p.y / grid_size;
    //printf("the grid_x is %d, the grid_y is %d\n",grid_x, grid_y);  


    int i,j;
    for (i = 0; i< reduced_n_proc; i++){
        if (grid_x < partition_grids[i+1])
            break;
    }
    for (j = 0; j< reduced_n_proc; j++){
        if (grid_y < partition_grids[j+1])
            break;
    }

    //printf("i=%d, j=%d\n", i,j);

    return i+j*sqrt(reduced_n_proc);

}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}


int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
