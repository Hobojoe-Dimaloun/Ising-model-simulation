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

//Screen dimension constants

static void randLattice(int *lattice, int x, int y, int z);


const double gBoltzmann = 1.38064852E-23;
const double gPi = 3.14159;

int main(int argc,char **argv)
{
    int numOfNodes;
    int debug = 0;  // Debug
    int x = 4; // x - dimension lattice sites
    int y = 4; // y - dimension lattice sites
    int z = 1; // z - dimension latties sites
    double temperature = 0.0001; // Kelvin
    int J = 1;
    //
    // Handle input arguments
    //
    int numOfThreads;
    for(int i = 0; i < argc; i++)
    {
        switch (i)
        {
            case 0: numOfThreads = omp_get_max_threads(); break;
            case 1: numOfThreads =  atoi(argv[i]); break;
            case 2: numOfNodes = *argv[i]; break;
            case 3: debug = 1; break;
            default: break;
        }
    }

    omp_set_num_threads(numOfThreads);
    if(debug == 1)
    {
        printf("Warning! Debug mode entered. Press any key to continue...\n");
        getchar();
    }

    printf("Ising model begin.\n");

    int *ising_lattice = NULL;

    //
    // Allocate lattice memory
    //

    if(( ising_lattice = calloc(x * y * z, sizeof(*ising_lattice) )) == NULL)
    {
        printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
    }

    //
    // Assign random spins to the lattice, ie 1 or -1
    //

    long int loop = 0;

    randLattice(ising_lattice,x,y,z);

    gsl_rng *rndarray[numOfThreads];
    //printf("%d\n", numOfThreads);
    for(int i = 0; i<numOfThreads; i++)
    {
        rndarray[i] = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rndarray[i],i);
    }

    //
    //
    //
    bool quit = false;

    double time = omp_get_wtime();

    FILE *output=fopen("output.txt","w");

    fprintf(output,"%d\n",x);
    fprintf(output,"%d\n",y);

    if(output == NULL)
    {
        exit(EXIT_FAILURE);
    }

	//While application is running
	while( !quit && (loop/1000 <= 100) )
	{

        printf("loop %ld\n",loop);
        #pragma omp parallel
        {
            gsl_rng *r = rndarray[omp_get_thread_num()];
            #pragma omp single
            {
    		//Handle events on queue

                if(loop%1000 == 0)
                {
                    int sum = 0;
                    for(int i = 0; i < x*y; i++)
                    {
                        fprintf(output,"%d ",ising_lattice[i]);
                        sum += ising_lattice[i];
                    }
                    fprintf(output,"\n");

                    fprintf(output,"%d\n",sum);

                }
            }
            #pragma omp barrier

            int location = (double)gsl_rng_uniform(r)* x * y * (omp_get_thread_num()+1) / (double)numOfThreads;
            int icolumn = location%y; // the remainder is the column number
            int irow = (location - icolumn)/x  ;// number of rows
            int temp;


            double dE = 0, E1=0, E2=0;
            //
            // Calculate energy of cell above
            //
            (irow == 0) ? (temp = x-1) : (temp = irow-1);

            E1+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell below
            //
            (irow == (x-1)) ? (temp = 0) : (temp = irow+1);

            E1+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell left
            //

            (icolumn == 0) ? (temp = y-1) : (temp = icolumn-1);
            E1+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            //
            // Calculate energy of cell right
            //
            (icolumn == (y-1)) ? (temp = 0) : (temp = icolumn+1);
            E1+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            //
            // Flip the spin
            //

            ising_lattice[location] *=-1 ;

            (irow == 0) ? (temp = x-1) : (temp = irow-1);

            E2+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell below
            //
            (irow == (x-1)) ? (temp = 0) : (temp = irow+1);

            E2+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell left
            //

            (icolumn == 0) ? (temp = y-1) : (temp = icolumn-1);
            E2+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            //
            // Calculate energy of cell right
            //
            (icolumn == (y-1)) ? (temp = 0) : (temp = icolumn+1);
            E2+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            dE= E2 - E1;

            if( dE > 0)
            {

                double beta =1/temperature ;

                double prob = exp(beta *dE*(-1)) ;
            //    printf(" prob = %g\n",prob);

                if( prob < gsl_rng_uniform(r))
                {
                //    printf(" t = %g\n",dE);
                    ising_lattice[location]*= -1;

                }
            }

        }

        loop ++;
    }
    time = omp_get_wtime() - time;
    printf("Time = %lf\n",time);

    return 0;
}

static void randLattice(int *lattice, int x, int y, int z)
{
    for(int i = 0; i < x * y * z; i ++)
    {
        if((rand()/(double)RAND_MAX) <0.5)
        {
            lattice[i] = 1;
        }
        else
        {
            lattice[i] = -1;
        }
    }
}
