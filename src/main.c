/*********************
* Date of creation: 02/11/17
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors:
**************************************
* Change History
**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <mpi.h>

#include "ising_functions.h"


#define MASTER  0

const double gBoltzmann = 1.38064852E-23;
const double gPi = 3.14159;

int main(int argc,char **argv)
{
    int numOfNodes=1;
    int debug = 0;  // Debug
    int x = 200; // x - dimension lattice sites
    int y = 200; // y - dimension lattice sites
    int z = 1; // z - dimension latties sites
    double temperature = 0.0001; // Kelvin
    double beta = 1/temperature;
    int J = 1; // anti/ferromagnetic
    int MPI_error = 0;
    //
    // Handle input arguments
    //
    int numOfThreads;
    for(int i = 0; i < argc; i++)
    {
        switch (i)
        {
        /*    case 0: numOfThreads = omp_get_max_threads(); break;
            case 1: numOfThreads =  atoi(argv[i]); break;
            case 2: debug = atoi(argv[i]); break;*/
            default:  numOfThreads = omp_get_max_threads(); break;
        }
    }

    if(debug == 1)
    {
        printf("Warning! Debug mode entered. Press any key to continue...\n");
        getchar();
    }

    MPI_Status status;
    /*****************INITIALISE****************************/

    //
    // Initialise MPI and get number of nodes.
    // Determine if even computation is possible
    // If possible, get rank of process and calculate chunksize
    // Set number of threads each node will use
    //
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfNodes);
    if( (x * y * z) % (numOfThreads * numOfNodes) != 0 )
    {
        printf("Exit: Grid size not evenly divisable by number of threads times nodes\n");
        MPI_Abort(MPI_COMM_WORLD, MPI_error);
        return 0;
    }

    int taskid;

    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    int chunksize = x * y * z / ( numOfNodes * numOfThreads );

    omp_set_num_threads(numOfThreads);
    /*****************INITIALISE END****************************/


    printf("Ising model begin rank %d.\n", taskid);
    fflush(stdout);



    /*****************INITIALISE FILE****************************/
    FILE *output = NULL ;
    if(taskid == MASTER)
    {
        output=fopen("output.txt","w");

        fprintf(output,"%d\n",x);
        fprintf(output,"%d\n",y);

        if(output == NULL)
        {
            /*int fail = 1;
            for(int task = 1; task<numOfNodes; task++)
            {
                MPI_Send(&fail,1,MPI_INT,task,messageTag[0],MPI_COMM_WORLD);
            }*/
            MPI_Abort(MPI_COMM_WORLD,MPI_error);
            exit(EXIT_FAILURE);
        }
    /*    else
        {
            int fail = 0;
            for(int task = 1; task<numOfNodes; task++)
            {
                MPI_Send(&fail,1,MPI_INT,task,messageTag[0],MPI_COMM_WORLD);
            }
        }*/
    }
/*    else
    {
        int fail = 0;
        MPI_Recv(&fail,1,MPI_INT,messageSource,messageTag[0],MPI_COMM_WORLD,&status);

        if(fail != 0)
        {
            MPI_Abort();
            exit(EXIT_FAILURE);
        }
    }*/

    /*****************INITIALISE FILE END****************************/




    /*****************MEMORY MANAGEMENT****************************/

    //
    // Allocate lattice memory
    //
    int *ising_lattice = NULL;

    if(( ising_lattice = calloc(x * y * z, sizeof(*ising_lattice) )) == NULL)
    {
        printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD,MPI_error);
        exit(EXIT_FAILURE);
    }

    //
    //  Allocate random number generators
    //

    gsl_rng *rndarray[numOfThreads];

    #pragma omp parallel for
    for(int i = 0; i<numOfThreads; i++)
    {
        rndarray[i] = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rndarray[i], i * taskid);
    }


    int messageTag[2]=
    {
        0,  // File status message
        1   // Lattice transfer message
    };
    int messageSource = 0;


    if(taskid == MASTER)
    {
        //
        // Assign random spins to the lattice, ie 1 or -1
        //
        randLattice(ising_lattice,x,y,z);
        for(int task = 1; task<numOfNodes; task++)
        {
            MPI_Send(&ising_lattice[0],x*y*z,MPI_INT,task,messageTag[1],MPI_COMM_WORLD);
            printf("Sent %d elements to rank %d\n",x*y*z,task);
            fflush(stdout);

        }
    }
    else
    {
        MPI_Recv(&ising_lattice[0],x*y*z,MPI_INT,messageSource,messageTag[1],MPI_COMM_WORLD,&status);
        printf("Rank %d recieved %d elements from task %d\n",taskid,x*y*z,messageSource);

    }
    /*****************MEMORY MANAGEMENT END****************************/

    long int loop = 0;


    //
    //
    //

    double time = omp_get_wtime();

	//While application is running
    int loopmax=1000/(numOfNodes*numOfThreads);
	while((loop/10000 <= loopmax) )
	{
        if(loop == 0 && taskid == MASTER)
        {
            printf("Approximate Kbytes : %lf  Mbytes : %lf  Gbytes : %lf\n",(double)x*y*loopmax*8/(double)1024, (double)x*y*loopmax*8/(double)(1024*1024), (double)x*y*loopmax*8/(double)(1024*1024*1024));
            fflush(stdout);
        }
        if(loop%10000 == 0)
        {
            printf("rank %d loop %ld\n",taskid,loop/10000);
            fflush(stdout);
        }
        #pragma omp parallel
        {
            gsl_rng *r = rndarray[omp_get_thread_num()];
            #pragma omp single
            {
    		//Handle events on queue

                if(loop%10000 == 0 && taskid == MASTER)
                {
                    int sum = 0;
                    for(int i = 0; i < x*y; i++)
                    {
                        fprintf(output,"%d ",ising_lattice[i]);
                        sum += ising_lattice[i];
                    }
                    fprintf(output,"\n");

                    fprintf(output,"%d\n",sum);
                    fflush(output);
                }
            }
            #pragma omp barrier

            int location = (double)gsl_rng_uniform(r)* x * y * (omp_get_thread_num()+1) / (double)numOfThreads;

            int E1=0, E2=0;
            //
            // Calculate spin energy
            //
            E1 = energy_calculation(ising_lattice, location, x,  y, J);
            //
            // Flip spin
            //
            ising_lattice[location] *=-1 ;
            //
            // Calcualte new spin-flipped energy
            //
            E2= energy_calculation(ising_lattice, location, x,  y, J);
            //
            // Keep spin flip with probability determined in spin_flip_check()
            //
            ising_lattice[location]*=spin_flip_check(E1, E2, beta, r);

        }

        loop ++;
    }
    time = omp_get_wtime() - time;
    printf("Time = %lf\n",time);

    MPI_Finalize();
    return 0;
}
