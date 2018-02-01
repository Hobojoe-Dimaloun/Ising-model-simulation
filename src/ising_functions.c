
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <omp.h>

#include "ising_functions.h"

//
// Randomises lattice
//

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

void energy_comparison(int x, int y, gsl_rng *r, int *ising_lattice, int beta, int J,int numOfThreads, int numOfNodes,int chunksize)
{
    int location = (double)gsl_rng_uniform(r)*chunksize*y* (omp_get_thread_num()+1) / (double)numOfThreads;

    int E1=0, E2=0;
    //
    // Calculate spin energy
    //
    E1 = energy_calculation(ising_lattice, location, x/numOfNodes,  y, J);
    //
    // Flip spin
    //
    ising_lattice[location] *=-1 ;
    //
    // Calcualte new spin-flipped energy
    //
    E2= energy_calculation(ising_lattice, location, x/numOfNodes,  y, J);
    //
    // Keep spin flip with probability determined in spin_flip_check()
    //
    ising_lattice[location]*=spin_flip_check(E1, E2, beta, r);
}


//
// Calculates the energy due to the surrounding cells
//

int energy_calculation(int *lattice, int location, int x, int y, int J)
{

    int icolumn = location%y; // the remainder is the column number
    int irow = (location - icolumn)/x  ;// number of rows
    int temp;
    int energy = 0;
    //
    // Calculate energy of cell above
    //
    (irow == 0) ? (temp = x-1) : (temp = irow-1);

    energy+= J*lattice[location] * lattice[temp*x + icolumn]*(-1);

    //
    // Calculate energy of cell below
    //
    (irow == (x-1)) ? (temp = 0) : (temp = irow+1);

    energy+= J*lattice[location] * lattice[temp*x + icolumn]*(-1);

    //
    // Calculate energy of cell left
    //

    (icolumn == 0) ? (temp = y-1) : (temp = icolumn-1);
    energy+= J*lattice[location] * lattice[irow*x + temp]*(-1);

    //
    // Calculate energy of cell right
    //
    (icolumn == (y-1)) ? (temp = 0) : (temp = icolumn+1);
    energy+= J*lattice[location] * lattice[irow*x + temp]*(-1);

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
