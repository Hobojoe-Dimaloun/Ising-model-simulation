
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>


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
