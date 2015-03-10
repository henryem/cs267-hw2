#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
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
        printf( "-p <int> to set the (maximum) number of threads used\n");
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
    printf("Using %d grid squares of side-length %f for %d particles.\n", num_grid_squares_per_side*num_grid_squares_per_side, grid_square_size, n);
    //
    //  set up MPI
    //
    int n_proc, 
        rank, 
        reduced_n_proc, 
        n_proc_per_side, 
        num_grid_squares_per_side_per_proc;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    

    // Reduced processors:
    // In order to equally divide area into processors, number of processors
    // must be a square number. Here we recalculate the number of processors.
    int i;
    //printf("Using %d of processors\n", n_proc);
    for (i=0; ;i++){
        if (i*i > n_proc)
            break;
    }
    n_proc_per_side = i-1;
    reduced_n_proc = (i-1)*(i-1);
    //printf("Num of processor per side is %d\n", n_proc_per_side);
    //printf("Reduce processors to %d\n", reduced_n_proc);


    // Decide num of grids per side per processor
    if (num_grid_squares_per_side % n_proc_per_side ==0)
        num_grid_squares_per_side_per_proc = num_grid_squares_per_side / n_proc_per_side;
    else
        num_grid_squares_per_side_per_proc = num_grid_squares_per_side / n_proc_per_side + 1;
//    printf("Num of grids per side per processor is %d\n", num_grid_squares_per_side_per_proc);

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




    // Calculate the parition_grids, a map to allocate pariticle
    int partition_grids[n_proc_per_side+1];
    partition_grids[n_proc_per_side] = INT_MAX;

    for (int i = 0; i < n_proc_per_side; i++){
        partition_grids[i] = i * num_grid_squares_per_side_per_proc;
        //printf("partition_grids %d is %d\n", i, partition_grids[i] );
    }
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //

    if( rank == 0 ) {
        
        particles = init_particles(n).release()->data();
        
        // p_buf: particle buffer that store particles for all processor
        // p_counter: indicate how many particles are currently in processor i
        particle_t *p_buf = (particle_t*) malloc (reduced_n_proc* n * sizeof(particle_t));
        int p_counter[reduced_n_proc];
        // Initiation
        for (int i = 0; i< reduced_n_proc; i++){
            p_counter[i] = 0;
        }

        for (int i = 0; i < n; i++){
            int proc_no = find_proc_no(particles[i], reduced_n_proc, partition_grids, grid_square_size);
            p_buf[proc_no*n+p_counter[proc_no]] = particles[i];
            p_counter[proc_no] +=1;
        }


        // b_buf: particle buffer that store particles in boundry region of each processor.
        // b_counter: indicate how many particles are currently in processor boundry region.
        particle_t *b_buf = (particle_t*) malloc (reduced_n_proc * n * sizeof(particle_t));
        int b_counter[reduced_n_proc];
        int max_case = 3;
        int *boundry_proc = (int *) malloc (max_case* sizeof(int));

        // Initiation
        for (int i = 0; i< reduced_n_proc; i++){
            b_counter[i] = 0;
        }

        for (int i = 0; i < n; i++){
            int count;
            count = find_boundry_proc(particles[i], reduced_n_proc, partition_grids, grid_square_size, boundry_proc);
            for (int j = 0; j<count; j++){
                int proc_no = boundry_proc[j];
                //printf("add partical into proc_no: %d\n", proc_no);
                b_buf[proc_no*n+b_counter[proc_no]] =particles[i];
                b_counter[proc_no] +=1;
            }
        } 

        MPI_Request request, request2;
        for (int i = 0; i< reduced_n_proc; i++){
            MPI_Isend(&p_buf[i*n],p_counter[i],PARTICLE,i,0,MPI_COMM_WORLD, &request);   // Using tag value 0 to indicate it's particle inside the processor
            MPI_Isend(&b_buf[i*n],b_counter[i],PARTICLE,i,1,MPI_COMM_WORLD, &request2);  // Using tag value 1 to indicate it's particle at the boundry of processor
        }

        free(p_buf);
        free(b_buf);
    }


    // Need to remove later
//    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    
    Stats local_stats;
    Stats global_stats;

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );


    // Each processor receive particles from processor zero.
    int count_inside, count_boundry;


    Grid * p_grid = new Grid(size, num_grid_squares_per_side);

    particle_t *p_inside = (particle_t *) malloc(n*sizeof(particle_t));
    particle_t *p_boundry = (particle_t *) malloc(n*sizeof(particle_t));

    if (rank < reduced_n_proc){
        MPI_Status status, status2;

        MPI_Recv(p_inside,n,PARTICLE,0,0,MPI_COMM_WORLD,&status);
        MPI_Recv(p_boundry,n,PARTICLE,0,1,MPI_COMM_WORLD,&status2);
        
        MPI_Get_count(&status, PARTICLE,&count_inside);
        MPI_Get_count(&status2, PARTICLE,&count_boundry);

        //Debug
        //printf("I'm processor %d, I received %d inside particles \n",rank,count_inside);
        //printf("I'm processor %d, I received %d boundry particles \n",rank,count_boundry);
        //printf("I'm processor %d, the first particle I received is (%f,%f)\n",rank,p_boundry[0].x,p_boundry[0].y);

    }

/*
    for (int i = 0; i < count_inside; i++)
        p_grid->add(p_inside[i]);
*/

    for( int step = 0; step < 1000/*NSTEPS*/; step++ )
    {
        if (rank < reduced_n_proc){

            //Calculate force
            for (int i = 0; i< count_inside; i++){
                //printf("p%d:(%f,%f)\n",rank,p_inside[i].x, p_inside[i].y);
                particle_t& p = p_inside[i];
                p.ax = p.ay = 0;

/*
                std::unique_ptr<SimpleIterator<particle_t&> > neighbors = p_grid->neighbor_iterator(p);

                while (neighbors->hasNext()) {
                    particle_t& neighbor = neighbors->next();
                    apply_force(p, neighbor, local_stats);
                }
*/

                for (int j = 0; j < count_inside; j++){
                    if (p_inside[j].x == p.x && p_inside[j].y == p.y)
                        continue;
                    else{
                        apply_force(p,p_inside[j], local_stats);
                    }
                }


                int boundry_never_use[3];
                if (find_boundry_proc(p, reduced_n_proc, partition_grids, grid_square_size, boundry_never_use) >0){
                    for (int j = 0; j < count_boundry; j++){
                        apply_force(p, p_boundry[j], local_stats);
                    }
                }

            }
            // Move the particles
            for( int i = 0; i < count_inside; i++ ){
                move( p_inside[i] );
                //printf("p%d: after moving, dump.\n", rank);
                //dump_particle(p_inside,count_inside, rank);
            }

            // Scan if particle outside boundry

            particle_t *p_buf = (particle_t*) malloc (reduced_n_proc* n * sizeof(particle_t));
            int p_counter[reduced_n_proc];
            for (int i = 0; i< reduced_n_proc; i++){
                p_counter[i] = 0;
            }

            particle_t *b_buf = (particle_t*) malloc (reduced_n_proc * n * sizeof(particle_t));
            int b_counter[reduced_n_proc];
            int max_case = 3;
            int *boundry_proc = (int *) malloc (max_case* sizeof(int));
            for (int i = 0; i< reduced_n_proc; i++){
                b_counter[i] = 0;
            }

            //Calculate which particles are going to send out to other processor
            for (int i = 0; i < count_inside; i++){
                // Check if the particle is still inside the processor
                int current_proc_no = find_proc_no(p_inside[i], reduced_n_proc, partition_grids, grid_square_size); 
                // The particle is outside the processor
                //printf("Current_proc_no: %d, rank: %d\n", current_proc_no, rank);
                if (current_proc_no != rank){
                    p_buf[current_proc_no*n+p_counter[current_proc_no]] = p_inside[i];
                    p_counter[current_proc_no] +=1;

                    //printf("p%d:move particle (%f,%f) into proc_no: %d\n", rank,p_inside[i].x, p_inside[i].y, current_proc_no);
 
                    // Mark p_inside[i] as removed (-1,-1), remove in the future
                    p_inside[i].x = -1;
                    p_inside[i].y = -1;
                }
                // The particle is inside the processor, check if boundry particle
                else{
                    // Notify other processor it's boundry particle now
                    int count;
                    count = find_boundry_proc(p_inside[i], reduced_n_proc, partition_grids, grid_square_size, boundry_proc);
                    for (int j = 0; j<count; j++){
                        int proc_no = boundry_proc[j];
                        
//                        printf("p%d: throw boundry (%f,%f) to %d\n", rank, p_inside[i].x, p_inside[i].y, proc_no);
                        
                        b_buf[proc_no*n+b_counter[proc_no]] = p_inside[i];
                        b_counter[proc_no] +=1;
                    }
                }
            }

            // Removing particles
            for (int i = 0; i< count_inside; i++){
                if (p_inside[i].x == -1 && p_inside[i].y == -1){
                    p_inside[i] = p_inside[count_inside-1];
                    count_inside -=1;
                    i -=1;

                }
            }

            // Sending out the particles
            MPI_Request request, request2;
            for (int i = 0; i< reduced_n_proc; i++){
//                printf("p%d: send %d inside particles to p%d\n",rank, p_counter[i], i);
//                printf("p%d: send %d boundry particles to p%d\n",rank, b_counter[i], i);
                MPI_Isend(&p_buf[i*n],p_counter[i],PARTICLE,i,0,MPI_COMM_WORLD, &request);   // Using tag value 0 to indicate it's particle inside the processor
                MPI_Isend(&b_buf[i*n],b_counter[i],PARTICLE,i,1,MPI_COMM_WORLD, &request2);  // Using tag value 1 to indicate it's particle at the boundry of processor
            }            

            // Receiving particles
            MPI_Status status, status2;
            // Each processor receive particles from processor zero.
            int total_count_inside = 0;
            int total_count_boundry = 0;

            particle_t *receive_inside = (particle_t *) malloc(n*sizeof(particle_t));
            particle_t *receive_boundry = (particle_t *) malloc(n*sizeof(particle_t));

            for (int i = 0; i < reduced_n_proc; i++){
                int receive_count_inside = 0; 
                int receive_count_boundry = 0;
                MPI_Recv(&receive_inside[total_count_inside],n,PARTICLE,i,0,MPI_COMM_WORLD,&status);
                MPI_Recv(&receive_boundry[total_count_boundry],n,PARTICLE,i,1,MPI_COMM_WORLD,&status2);
            
                MPI_Get_count(&status, PARTICLE,&receive_count_inside);
                MPI_Get_count(&status2, PARTICLE,&receive_count_boundry);
             
                total_count_inside += receive_count_inside;
                total_count_boundry += receive_count_boundry;

            }

            //Update p_inside
            for (int i = 0; i < total_count_inside; i++){
                p_inside[count_inside] = receive_inside[i];
                count_inside +=1;
            }

            count_boundry = 0;
            for (int i = 0; i < total_count_boundry; i++){
                p_boundry[count_boundry] = receive_boundry[i];
                count_boundry+=1;
            }

            free(p_buf);
            free(b_buf);
            free(receive_inside);
            free(receive_boundry);
            free(boundry_proc);

        }

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if(!fast)
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //

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


    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
