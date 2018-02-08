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
#include <stdbool.h>

#include "ising_functions.h"


#define MASTER  0

int gNumOfNodes ;
int gNumOfthreads = 16;

int main(int argc,char **argv)
{
    int numOfSpinFlips = 1E6;
    int x_Max = 32; // x - dimension lattice sites
    int y_Max = 32; // y - dimension lattice sites
    int z_Max = 32; // z - dimension latties sites
    double temperature =0.001; // Kelvin
    double tempmax =6;
    double tempdelta =0.1;
    int J = 1; // anti/ferromagnetic
    int MPI_error = 0;
    //
    // Handle input arguments
    //
    bool isFirstLoop = true ;

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
    if((x_Max*y_Max*z_Max)%gNumOfNodes != 0 || (((x_Max*y_Max*z_Max)%gNumOfNodes))%gNumOfthreads!=0)
    {
        printf("Exit: Grid size not evenly divisable by number of threads times nodes\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, MPI_error);
        return 0;
    }

    int taskid;

    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    //
    // Determine number of cells per core
    //

    int chunksize = x_Max*y_Max*z_Max  / gNumOfNodes ;

    omp_set_num_threads(gNumOfthreads);
    /*****************INITIALISE END****************************/



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
            printf("Started\n" );
            fflush(stdout);
            /*****************INITIALISE FILE****************************/
            FILE *output = NULL ;

            output=fopen("output.txt","a");

            if(output == NULL)
            {
                MPI_Abort(MPI_COMM_WORLD,MPI_error);
                exit(EXIT_FAILURE);
            }

            fprintf(output,"%d\n",x_Max);
            fprintf(output,"%d\n",y_Max);
            fflush(output);

            /*****************INITIALISE FILE END****************************/


            /*****************MEMORY MANAGEMENT****************************/

            //
            // Allocate lattice memory
            //
            int *ising_lattice = NULL;

            if(( ising_lattice = calloc(x_Max * y_Max * z_Max, sizeof(*ising_lattice) )) == NULL)
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

            if(( ising_lattice_node_boundaries = calloc(x_Max*z_Max*2, sizeof(*ising_lattice_node_boundaries) )) == NULL)
            {
                printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD,MPI_error);
                exit(EXIT_FAILURE);
            }

            //
            // Assign random spins to the lattice, ie 1 or -1
            //
            randLattice(ising_lattice,x_Max,y_Max,z_Max);

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

            }
            for( int write = 0; write < chunksize; write++)
            {
                ising_lattice_node_segment[write] =ising_lattice[write];

            }
            /*****************MEMORY MANAGEMENT END****************************/

            while(temperature < tempmax)
            {
                double beta=1/temperature;
                int coreSpinFlip =numOfSpinFlips/(gNumOfNodes*gNumOfthreads);

                if(isFirstLoop == true)
                {
                    coreSpinFlip = 1E7/(gNumOfNodes*gNumOfthreads);
                    isFirstLoop = false ;
                }
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

                        MPI_Send(&ising_lattice_node_segment[0],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);
                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x_Max*z_Max],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x_Max*z_Max],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);
                    }
                    else
                    {

                        int send ;

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x_Max*z_Max],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x_Max*z_Max],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Send own top boundary
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);

                        MPI_Send(&ising_lattice_node_segment[0],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                    }

                    #pragma omp parallel
                    {
                        int *ising_lattice_core_boundaries = NULL;
                        if(( ising_lattice_core_boundaries = calloc((y_Max/gNumOfNodes) * z_Max * 2, sizeof(*ising_lattice_core_boundaries) )) == NULL)
                        {
                            printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                            fflush(stdout);
                            MPI_Abort(MPI_COMM_WORLD,MPI_error);
                            exit(EXIT_FAILURE);
                        }

                        for(int i = 0; i < y_Max/gNumOfNodes; i++)
                        {
                            for(int j = 0; j < z_Max; j++)
                            {
                                int left, right, offset;

                                //
                                // Get left core boundary
                                //
                                (omp_get_thread_num() == 0) ? (left = gNumOfthreads -1 ) : (left = omp_get_thread_num()-1);
                                offset = x_Max/gNumOfthreads*(left+1);
                                ising_lattice_core_boundaries[j*y_Max/gNumOfNodes +i] = ising_lattice_node_segment[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  offset-1,  i, j)];

                                //
                                // Get right core boundary
                                //
                                (omp_get_thread_num() == gNumOfthreads-1) ? (right = 0) : (right = omp_get_thread_num()+1);

                                offset = x_Max/gNumOfthreads*(right);
                                ising_lattice_core_boundaries[j*y_Max/gNumOfNodes +i +y_Max/gNumOfNodes*z_Max]=ising_lattice_node_segment[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  offset,  i, j)];

                            }
                        }

                        //
                        // Calculate spin flip
                        //
                        energy_comparison(x_Max, y_Max, z_Max, rndarray[omp_get_thread_num()], ising_lattice_node_segment, ising_lattice_core_boundaries,ising_lattice_node_boundaries, beta,  J, chunksize);

                        free(ising_lattice_core_boundaries);
                    }

                }

                //
                // Reconstruct the lattice, ready for file write
                //
                for(int task = 1; task<gNumOfNodes; task++)
                {
                    MPI_Recv(&ising_lattice[task*chunksize],chunksize,MPI_INT,task,messageTag[2],MPI_COMM_WORLD,&status);


                }

                for( int write = 0; write < chunksize; write++)
                {
                    ising_lattice[write] =ising_lattice_node_segment[write];

                }

                double magnetisation = 0;
                double magnetisation2 = 0;
                double energy = 0;
                double energy2 =0;


                for( int write = 0; write <x_Max*y_Max*z_Max; write++)
                {
                    magnetisation +=ising_lattice[write];

                    magnetisation2 += (ising_lattice[write]*ising_lattice[write]);


                    int x,y,z;
                    coordinates_from_linear_index(write,x_Max,y_Max,&x, &y, &z);


                    int temp=0;
                    int temp2=0;
                    //
                    // Calculate energy of cell above
                    //
                    (y == 0) ? (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  x ,  y_Max -1 ,  z)]*(-1) )
                                : (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  x ,  y-1 ,  z)]*(-1));
                                temp2 +=temp;

                    //
                    // Calculate energy of cell below
                    //
                    (y == (y_Max-1)) ? (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  x , 0 ,  z)]*(-1))
                                    : (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  x ,  y+1 ,  z)]*(-1));
                                    temp2 +=temp;

                    //
                    // Calculate energy of cell left
                    //

                    (x == 0) ? (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  x_Max -1 ,  y ,  z)]*(-1))
                                    : (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  x-1 ,  y ,  z)]*(-1));
                    //
                    // Calculate energy of cell right
                    //
                    temp2 +=temp;

                    (x == (x_Max-1)) ? (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  0 ,  y ,  z)]*(-1))
                                    : (temp = J*ising_lattice[write] * ising_lattice[linear_index_from_coordinates(x_Max, y_Max,  x+1 ,  y ,  z)]*(-1));
                    temp2 +=temp;

                    if(z_Max>1)
                    {

                        //
                        // Calculate energy of cell infront
                        //
                        (z == 0) ? (energy += J*ising_lattice[write  ] * ising_lattice[linear_index_from_coordinates(x_Max,x,  x ,  y,  z_Max -1)]*(-1))
                                            : (energy += J*ising_lattice[write ] * ising_lattice[linear_index_from_coordinates(x_Max,y_Max,  x  ,  y,  z-1)]*(-1));
                                            temp2 +=temp;
                        //
                        // Calculate energy of cell behind
                        //
                        (z == (z_Max-1)) ? (energy += J*ising_lattice[write ] * ising_lattice[linear_index_from_coordinates(x_Max,y_Max,  x,  y,  0)]*(-1))
                                        : (energy += J*ising_lattice[write ] * ising_lattice[linear_index_from_coordinates(x_Max,y_Max,  x ,  y,  z+1)]*(-1));
                                        temp2 +=temp;

                    }

                    energy +=0.5*(double)temp2 ;
                    energy2+= pow(0.5*(double)temp2,2);

                }
                magnetisation/=(double)(x_Max*y_Max*z_Max);
                magnetisation2/=(double)(x_Max*y_Max*z_Max);
                energy2/=(double)(x_Max*y_Max*z_Max);
                energy/=(double)(x_Max*y_Max*z_Max);
                double Susceptability = (magnetisation2 - magnetisation*magnetisation)*beta;
                double heat_capacity = (energy2 - energy*energy)*beta*beta;//(x_Max*y_Max*z_Max);
                fprintf(output,"%lf\t%lf\t%lf\t%lf\t%lf\n",temperature,energy,fabs(magnetisation),Susceptability,heat_capacity);
                fflush(output);
            /*    for(int write = 0; write<x*y*z; write++)
                {
                    fprintf(output,"%d ",ising_lattice[write]);
                //    printf("Reviced %d elements from rank %d\n",chunksize,task);
                //    fflush(stdout);

            }
                fprintf(output,"\n");*/

                for(int task = 1; task<gNumOfNodes; task++)
                {
                    int check =1;
                    MPI_Send(&check,1,MPI_INT,task,messageTag[0],MPI_COMM_WORLD);

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

            if(( ising_lattice_node_boundaries = calloc(x_Max*z_Max*2, sizeof(*ising_lattice_node_boundaries) )) == NULL)
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

            /*****************MEMORY MANAGEMENT END****************************/
            while(temperature < tempmax)
            {
                double beta=1/temperature;
                int coreSpinFlip =numOfSpinFlips/(gNumOfNodes*gNumOfthreads);

                if(isFirstLoop == true)
                {
                    coreSpinFlip = 1E7/(gNumOfNodes*gNumOfthreads);
                    isFirstLoop = false ;
                }
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

                        MPI_Send(&ising_lattice_node_segment[0],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);
                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x_Max*z_Max],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x_Max*z_Max],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);
                    }
                    else
                    {

                        int send ;

                        //
                        // Recieve neighbours top boundary (our bottom)
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);
                        MPI_Recv(&ising_lattice_node_boundaries[x_Max*z_Max],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Recieve neighbours bottom boundary (our top)
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);
                        MPI_Recv(&ising_lattice_node_boundaries[0],x_Max*z_Max,MPI_INT,send,messageTag[3],MPI_COMM_WORLD,&status);

                        //
                        // Send own bottom boundary
                        //
                        (taskid == gNumOfNodes-1) ? (send = 0) : (send = taskid+1);

                        MPI_Send(&ising_lattice_node_segment[chunksize - x_Max*z_Max],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                        //
                        // Send own top boundary
                        //
                        (taskid == 0) ? (send = gNumOfNodes -1) : (send = taskid-1);

                        MPI_Send(&ising_lattice_node_segment[0],x_Max*z_Max,MPI_INT,send ,messageTag[3],MPI_COMM_WORLD);

                    }

                    #pragma omp parallel
                    {
                        int *ising_lattice_core_boundaries = NULL;
                        if(( ising_lattice_core_boundaries = calloc((y_Max/gNumOfNodes) * z_Max * 2, sizeof(*ising_lattice_core_boundaries) )) == NULL)
                        {
                            printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
                            fflush(stdout);
                            MPI_Abort(MPI_COMM_WORLD,MPI_error);
                            exit(EXIT_FAILURE);
                        }

                        for(int i = 0; i < y_Max/gNumOfNodes; i++)
                        {
                            for(int j = 0; j < z_Max; j++)
                            {
                                int left, right, offset;

                                //
                                // Get left core boundary
                                //
                                (omp_get_thread_num() == 0) ? (left = gNumOfthreads -1 ) : (left = omp_get_thread_num()-1);
                                offset = x_Max/gNumOfthreads*(left+1);
                                                                                                    //return x+x_Max*y+x_Max*y_Max*z;
                                ising_lattice_core_boundaries[j*y_Max/gNumOfNodes +i] = ising_lattice_node_segment[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  offset-1,  i, j)];

                                //
                                // Get right core boundary
                                //
                                (omp_get_thread_num() == gNumOfthreads-1) ? (right = 0) : (right = omp_get_thread_num()+1);

                                offset = x_Max/gNumOfthreads*(right);
                                ising_lattice_core_boundaries[j*y_Max/gNumOfNodes +i +y_Max/gNumOfNodes*z_Max]=ising_lattice_node_segment[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  offset,  i, j)];
                            }
                        }

                        //
                        // Calculate spin flip
                        //

                        energy_comparison(x_Max, y_Max, z_Max, rndarray[omp_get_thread_num()], ising_lattice_node_segment, ising_lattice_core_boundaries,ising_lattice_node_boundaries, beta,  J, chunksize);
                        free(ising_lattice_core_boundaries);

                    }


                }

                //
                // Send chunk for reconstruction. Wait for check bit.
                //

                int check;
                MPI_Send(&ising_lattice_node_segment[0],chunksize,MPI_INT,MASTER,messageTag[2],MPI_COMM_WORLD);
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
