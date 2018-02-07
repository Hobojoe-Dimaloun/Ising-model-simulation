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
int gNumOfNodes ;
int gNumOfthreads ;

int main(int argc,char **argv)
{
    int numOfSpinFlips = 1E5;
    int debug = 0;  // Debug
    int x = 32; // x - dimension lattice sites
    int y = 32; // y - dimension lattice sites
    int z = 1; // z - dimension latties sites
    double temperature =0.01; // Kelvin
    double tempmax =6;
    double tempdelta =0.01;
    int J = 1; // anti/ferromagnetic
    int MPI_error = 0;
    //
    // Handle input arguments
    //
    for(int i = 0; i < argc; i++)
    {
        switch (i)
        {
        /*    case 0: numOfThreads = omp_get_max_threads(); break;
            case 1: numOfThreads =  atoi(argv[i]); break;
            case 2: debug = atoi(argv[i]); break;*/
            default:  gNumOfthreads = 2;//omp_get_max_threads();
            break;
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
    MPI_Comm_size(MPI_COMM_WORLD, &gNumOfNodes);

    //
    // Find out if the grid is 2^n
    //
    if((x*y*z)%gNumOfNodes != 0 || (((x*y*z)%gNumOfNodes))%gNumOfthreads!=0)
    {
        printf("Exit: Grid size not evenly divisable by number of threads times nodes\n");
        MPI_Abort(MPI_COMM_WORLD, MPI_error);
        return 0;
    }

    int taskid;

    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    //
    // Determine number of cells per core
    //

    int chunksize = x*y*z  / gNumOfNodes ;

    omp_set_num_threads(gNumOfthreads);
    /*****************INITIALISE END****************************/


    //printf("Ising model begin rank %d.\n", taskid);
//    fflush(stdout);

    int messageTag[4]=
    {
        0,  // File status message
        1,  // Full lattice transfer message
        2,   // Partial lattice transfer message
        3   // Boundary lattice transfer message
    };



        /*************************MASTER TASK*********************************/
        if(taskid == MASTER)
        {
            printf("T = %lf\n",temperature );
            /*****************INITIALISE FILE****************************/
            FILE *output = NULL ;

            output=fopen("output.txt","a");

            if(output == NULL)
            {
                MPI_Abort(MPI_COMM_WORLD,MPI_error);
                exit(EXIT_FAILURE);
            }

            fprintf(output,"%d\n",x);
            fprintf(output,"%d\n",y);
            fflush(output);

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

            int *ising_lattice_node_segment = NULL;

            if(( ising_lattice_node_segment = calloc(chunksize, sizeof(*ising_lattice_node_segment) )) == NULL)
            {
                printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD,MPI_error);
                exit(EXIT_FAILURE);
            }

            int *ising_lattice_node_boundaries= NULL;

            if(( ising_lattice_node_boundaries = calloc(x*z*2, sizeof(*ising_lattice_node_boundaries) )) == NULL)
            {
                printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD,MPI_error);
                exit(EXIT_FAILURE);
            }

            //
            // Assign random spins to the lattice, ie 1 or -1
            //
            randLattice(ising_lattice,x,y,z);
          /* for( int write = 0; write < x*y*z; write++)
            {
                ising_lattice[write] =write;

            }*/

            //
            //  Allocate random number generators
            //
            gsl_rng *rndarray[gNumOfthreads];

            for(int i = 0; i<gNumOfthreads; i++)
            {
                rndarray[i] = gsl_rng_alloc(gsl_rng_taus);
                gsl_rng_set(rndarray[i], i);
            }



            for(int task = 1; task<gNumOfNodes; task++)
            {
                MPI_Send(&ising_lattice[task*chunksize],chunksize,MPI_INT,task,messageTag[1],MPI_COMM_WORLD);
            //    printf("Sent %d elements to rank %d\n",chunksize,task);
            //    fflush(stdout);

            }
            for( int write = 0; write < chunksize; write++)
            {
                ising_lattice_node_segment[write] =ising_lattice[write];

            }
            //memcpy(ising_lattice_node_segment,ising_lattice,chunksize);
            /*****************MEMORY MANAGEMENT END****************************/

        /*   {
                fprintf(output,"%d ",ising_lattice[write]);

            }
            fprintf(output,"\n");
            fflush(output);*/
            while(temperature < tempmax)
            {
                double beta=1/temperature;
                int coreSpinFlip =numOfSpinFlips/(gNumOfNodes*gNumOfthreads);

                for(int loop = 0; loop < coreSpinFlip; loop ++)
                {
                    //
                    // update boundaries
                    //
                    if(taskid%2)
                    {
                        //
                        // Send own top boundary
                        //
                        int send ;
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);

                        MPI_Send(&ising_lattice_node_segment[0],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);
                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);
                    }
                    else
                    {

                        int send ;

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Send own top boundary
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);

                        MPI_Send(&ising_lattice_node_segment[0],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                    }

                    #pragma omp parallel
                    {
                        int *ising_lattice_core_boundaries = NULL;
                        if(( ising_lattice_core_boundaries = calloc((y/gNumOfNodes) * z * 2, sizeof(*ising_lattice_core_boundaries) )) == NULL)
                        {
                            printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                            fflush(stdout);
                            MPI_Abort(MPI_COMM_WORLD,MPI_error);
                            exit(EXIT_FAILURE);
                        }

                        for(int i = 0; i < y/gNumOfNodes*z; i++)
                        {
                            int left, right;
                            //
                            // Get left core boundary
                            //
                            (omp_get_thread_num() == 0) ? (left = gNumOfthreads -1 ) : (left = omp_get_thread_num()-1);
                            ising_lattice_core_boundaries[i] = ising_lattice_node_segment[i * x + x/gNumOfthreads*(left+1)-1 ];
                            //
                            // Get right core boundary
                            //
                            (omp_get_thread_num() == gNumOfthreads-1) ? (right = 0) : (right = omp_get_thread_num()+1);
                            ising_lattice_core_boundaries[i+y/gNumOfNodes] =ising_lattice_node_segment[i * x + x/gNumOfthreads*right];
                        }



                        //
                        // Calculate spin flip
                        //
                        energy_comparison(x, y, rndarray[omp_get_thread_num()], ising_lattice_node_segment, ising_lattice_core_boundaries,ising_lattice_node_boundaries, beta,  J, chunksize);

                        free(ising_lattice_core_boundaries);
                    }

                }

                for(int task = 1; task<gNumOfNodes; task++)
                {
                    int check =1;
                    MPI_Recv(&ising_lattice[task*chunksize],chunksize,MPI_INT,task,messageTag[2],MPI_COMM_WORLD,&status);
                    //MPI_Send(&check,1,MPI_INT,task,messageTag[0],MPI_COMM_WORLD);
                //    printf("Reviced %d elements from rank %d\n",chunksize,task);
                //    fflush(stdout);

                }
                for( int write = 0; write < chunksize; write++)
                {
                    ising_lattice[write] =ising_lattice_node_segment[write];

                }

                double magnetisation = 0;
                double magnetisation2 = 0;
                double energy = 0;
                double energy2 =0;
                for( int write = 0; write < x*y*z; write++)
                {
                    magnetisation +=ising_lattice[write];
                //    printf("magnetisation %f, loop %d\n",magnetisation, write );
                    magnetisation2 += (ising_lattice[write]*ising_lattice[write]);

                    int icolumn = write%(x); // the remainder is the column number
                    int irow = (write - icolumn)/(x) ;// number of rows
                    int temp=0;
                    int temp2=0;

                    (irow == 0) ? (temp += J*ising_lattice[irow * x + icolumn] * ising_lattice[x*y*z -x +icolumn ]*(-1) )
                                : (temp += J*ising_lattice[irow * x + icolumn] * ising_lattice[(irow-1)*x + icolumn ]*(-1));
                                temp2 +=temp;
                //printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ], lattice[(irow-1)*x + icolumn + columnOffset]);

                    //          printf("%d : \n",J*lattice[irow * x + icolumn + columnOffset ], nodeBoundaries[columnOffset + icolumn ]);
                    // Calculate energy of cell below
                    //
                    (irow == (x-1)) ? (temp += J*ising_lattice[irow * x + icolumn  ] * ising_lattice[icolumn]*(-1))
                                    : (temp += J*ising_lattice[irow * x + icolumn  ] * ising_lattice[(irow+1)*x + icolumn ]*(-1));
                                    temp2 +=temp;
                //printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ], lattice[(irow+1)*x + icolumn + columnOffset]);
                    //
                    // Calculate energy of cell left
                    //

                    (icolumn == 0) ? (temp += J*ising_lattice[irow * x + icolumn   ] * ising_lattice[irow*x + x -1]*(-1))
                                    : (temp += J*ising_lattice[irow * x + icolumn  ] * ising_lattice[irow*x + icolumn -1]*(-1));
                //printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ],lattice[irow*x + icolumn + columnOffset -1]);
                    //
                    // Calculate energy of cell right
                    //
                    temp2 +=temp;
                    (icolumn == (y-1)) ? (temp += J*ising_lattice[irow * x + icolumn   ] * ising_lattice[irow*x]*(-1))
                                    : (temp += J*ising_lattice[irow * x + icolumn  ] * ising_lattice[irow*x + icolumn +1]*(-1));
                                    //printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ], lattice[irow*x + icolumn + columnOffset +1]);
                                    temp2 +=temp;
                                    energy +=0.5*(double)temp2 ;
                                    energy2+= pow(0.5*(double)temp2,2);

                }
                magnetisation/=(double)(x*y*z);
                magnetisation2/=(double)(x*y*z);
                energy2/=(double)(x*y*z);
                energy/=(double)(x*y*z);
                double Susceptability = (magnetisation2 - magnetisation*magnetisation)*beta;
                double heat_capacity = (energy2 - energy*energy)*beta*beta/(x*y*z);
                fprintf(output,"%lf\t%lf\t%lf\t%lf\n",temperature,energy,magnetisation,Susceptability,heat_capacity);
                fflush(output);

                for(int task = 1; task<gNumOfNodes; task++)
                {
                    int check =1;
                    MPI_Send(&check,1,MPI_INT,task,messageTag[0],MPI_COMM_WORLD);
                //    printf("Reviced %d elements from rank %d\n",chunksize,task);
                //    fflush(stdout);

                }
                temperature+=tempdelta;
            }
            fclose(output);
            free(ising_lattice_node_boundaries);
            free(ising_lattice_node_segment);
            free(ising_lattice);

        }


        /*************************MASTER TASK END*********************************/





        /*************************OTHER TASK*************************************/
        else
        {
            /*****************MEMORY MANAGEMENT****************************/
            int *ising_lattice_node_segment= NULL;

            if(( ising_lattice_node_segment = calloc(chunksize, sizeof(*ising_lattice_node_segment) )) == NULL)
            {
                printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD,MPI_error);
                exit(EXIT_FAILURE);
            }

            int *ising_lattice_node_boundaries= NULL;

            if(( ising_lattice_node_boundaries = calloc(x*z*2, sizeof(*ising_lattice_node_boundaries) )) == NULL)
            {
                printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD,MPI_error);
                exit(EXIT_FAILURE);
            }

            gsl_rng *rndarray[gNumOfthreads];

            for(int i = 0; i<gNumOfthreads; i++)
            {
                rndarray[i] = gsl_rng_alloc(gsl_rng_taus);
                gsl_rng_set(rndarray[i], i + gNumOfthreads*taskid);
            }

            MPI_Recv(&ising_lattice_node_segment[0],chunksize,MPI_INT,MASTER,messageTag[1],MPI_COMM_WORLD,&status);
        //    printf("Rank %d recieved %d elements from task %d\n",taskid,chunksize,MASTER);
            //fflush(stdout);

            /*****************MEMORY MANAGEMENT END****************************/
            while(temperature < tempmax)
            {
                double beta=1/temperature;
                int coreSpinFlip =numOfSpinFlips/(gNumOfNodes*gNumOfthreads);

                for(int loop = 0; loop < coreSpinFlip; loop ++)
                {
                        //printf("loop %d\n",loop);
                //        fflush(stdout);
                    //
                    // update boundaries
                    //
                    if(taskid%2)
                    {
                        //
                        // Send own top boundary
                        //
                        int send ;
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);

                        MPI_Send(&ising_lattice_node_segment[0],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);
                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);
                    }
                    else
                    {

                        int send ;

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x*z,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Send own top boundary
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);

                        MPI_Send(&ising_lattice_node_segment[0],x*z,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                    }

                    #pragma omp parallel
                    {
                        int *ising_lattice_core_boundaries = NULL;
                        if(( ising_lattice_core_boundaries = calloc((y/gNumOfNodes) *z* 2, sizeof(*ising_lattice_core_boundaries) )) == NULL)
                        {
                            printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                            fflush(stdout);
                            MPI_Abort(MPI_COMM_WORLD,MPI_error);
                            exit(EXIT_FAILURE);
                        }

                        for(int i = 0; i < y/gNumOfNodes *z; i++)
                        {
                            int left, right;
                            //
                            // Get left core boundary
                            //
                            (omp_get_thread_num() == 0) ? (left = gNumOfthreads -1 ) : (left = omp_get_thread_num()-1);
                            ising_lattice_core_boundaries[i] = ising_lattice_node_segment[i * x + x/gNumOfthreads*(left+1)-1 ];
                            //
                            // Get right core boundary
                            //
                            (omp_get_thread_num() == gNumOfthreads-1) ? (right = 0) : (right = omp_get_thread_num()+1);
                            ising_lattice_core_boundaries[i+y/gNumOfNodes] =ising_lattice_node_segment[i * x + x/gNumOfthreads*right];
                        }

                        //
                        // Calculate spin flip
                        //

                        energy_comparison(x, y, rndarray[omp_get_thread_num()], ising_lattice_node_segment, ising_lattice_core_boundaries,ising_lattice_node_boundaries, beta,  J, chunksize);
                        free(ising_lattice_core_boundaries);

                    }


                }


                int check;
                MPI_Send(&ising_lattice_node_segment[0],chunksize,MPI_INT,MASTER,messageTag[2],MPI_COMM_WORLD);
                //printf("Sent %d elements from rank %d\n",chunksize,taskid);
                MPI_Recv(&check,1,MPI_INT,MASTER,messageTag[0],MPI_COMM_WORLD,&status);
                temperature+=tempdelta;
            }


            free(ising_lattice_node_boundaries);
            free(ising_lattice_node_segment);
        }
        /*************************OTHER TASK END*************************************/




    MPI_Finalize();
    return 0;
}
