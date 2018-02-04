
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <omp.h>

#include "ising_functions.h"

//
// Randomises lattice
//

extern int gNumOfNodes;
extern int gNumOfthreads;

void randLattice(int *lattice, int x, int y, int z)
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

//
// Perform energy comparison
//

void energy_comparison(int x, int y, gsl_rng *r, int *ising_lattice, int * ising_lattice_core_boundaries,int *ising_lattice_node_boundaries, double beta, int J,int chunksize)
{
    int location = (double)gsl_rng_uniform(r) * (x /gNumOfthreads) * (y /gNumOfNodes);
    int columnOffset = (x/gNumOfthreads) * omp_get_thread_num();

    int icolumn = location%(x /gNumOfthreads); // the remainder is the column number
    int irow = (location - icolumn)/(x /gNumOfthreads);// number of rows
    int E1=0, E2=0;
    //
    // Calculate spin energy
    //
    E1 = energy_calculation(ising_lattice,ising_lattice_core_boundaries,ising_lattice_node_boundaries, location, x,  y, J);
//    printf("E1 %d\n",E1 );
    //
    // Flip spin
    //
    ising_lattice[irow * x + icolumn + columnOffset ] *=-1 ;
    //
    // Calcualte new spin-flipped energy
    //
    E2= energy_calculation(ising_lattice,ising_lattice_core_boundaries, ising_lattice_node_boundaries, location,x,  y, J);
    //printf("E2 %d\n",E2 );

    //
    // Keep spin flip with probability determined in spin_flip_check()
    //
    ising_lattice[irow * x + icolumn + columnOffset ]*=spin_flip_check(E1, E2, beta, r);
}


//
// Calculates the energy due to the surrounding cells
//

int energy_calculation(int *lattice, int *coreBoundaries, int *nodeBoundaries, int location, int x, int y, int J)
{

    int icolumn = location%(x /gNumOfthreads); // the remainder is the column number
    int irow = (location - icolumn)/(x /gNumOfthreads) ;// number of rows
    int temp;
    int energy = 0;
    int columnOffset = x/gNumOfthreads * omp_get_thread_num();
    //
    // Calculate energy of cell above
    //
    (irow == 0) ? (energy += J*lattice[irow * x + icolumn + columnOffset  ] * nodeBoundaries[columnOffset + icolumn ]*(-1) )
                : (energy += J*lattice[irow * x + icolumn + columnOffset ] * lattice[(irow-1)*x + icolumn + columnOffset]*(-1));

//printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ], lattice[(irow-1)*x + icolumn + columnOffset]);

    //          printf("%d : \n",J*lattice[irow * x + icolumn + columnOffset ], nodeBoundaries[columnOffset + icolumn ]);
    // Calculate energy of cell below
    //
    (irow == (x-1)) ? (energy += J*lattice[irow * x + icolumn + columnOffset  ] * nodeBoundaries[columnOffset + icolumn + x]*(-1))
                    : (energy += J*lattice[irow * x + icolumn + columnOffset ] * lattice[(irow+1)*x + icolumn + columnOffset]*(-1));

//printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ], lattice[(irow+1)*x + icolumn + columnOffset]);
    //
    // Calculate energy of cell left
    //

    (icolumn == 0) ? (energy += J*lattice[irow * x + icolumn + columnOffset  ] * coreBoundaries[icolumn]*(-1))
                    : (energy += J*lattice[irow * x + icolumn + columnOffset ] * lattice[irow*x + icolumn + columnOffset -1]*(-1));
//printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ],lattice[irow*x + icolumn + columnOffset -1]);
    //
    // Calculate energy of cell right
    //

    (icolumn == (y-1)) ? (energy += J*lattice[irow * x + icolumn + columnOffset  ] * coreBoundaries[icolumn + y/gNumOfNodes]*(-1))
                    : (energy += J*lattice[irow * x + icolumn + columnOffset ] * lattice[irow*x + icolumn + columnOffset +1]*(-1));
                    //printf("%d : %d\n",lattice[irow * x + icolumn + columnOffset ], lattice[irow*x + icolumn + columnOffset +1]);

    return energy;
}

int spin_flip_check(int E1, int E2, double beta, gsl_rng *r)
{

    int dE = E2 - E1;
    if( dE > 0)
    {
        double prob = exp(beta * (double)dE * (-1));
    //    printf(" prob = %g\n",prob);

        if( prob < gsl_rng_uniform(r))
        {
            return -1;
        }
        else
        {
            return 1;
        }
    }
    else
    {
        return 1;
    }
}
