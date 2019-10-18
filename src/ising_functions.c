
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

void energy_comparison(int x_Max, int y_Max, int z_Max, gsl_rng *r, int *ising_lattice, int * ising_lattice_core_boundaries,int *ising_lattice_node_boundaries, double beta, int J,int chunksize)
{
    int location = (double)gsl_rng_uniform(r) * (x_Max /gNumOfthreads) * (y_Max /gNumOfNodes) *z_Max;

    int E1=0, E2=0;

    int columnOffset = x_Max/gNumOfthreads * omp_get_thread_num();

    int x,y,z;
    coordinates_from_linear_index(location,x_Max/gNumOfthreads, y_Max/gNumOfNodes, &x, &y, &z);
    int offsetLocation = linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset,  y,  z);
    //
    // Calculate spin energy
    //
    E1 = energy_calculation(ising_lattice,ising_lattice_core_boundaries,ising_lattice_node_boundaries, location, x_Max,  y_Max, z_Max, J);
    //
    // Flip spin
    //
    ising_lattice[offsetLocation ] *=-1 ;
    //
    // Calcualte new spin-flipped energy
    //
    E2= energy_calculation(ising_lattice,ising_lattice_core_boundaries, ising_lattice_node_boundaries, location, x_Max,  y_Max, z_Max, J);
    //printf("E2 %d\n",E2 );

    //
    // Keep spin flip with probability determined in spin_flip_check()
    //
    ising_lattice[offsetLocation]*=spin_flip_check(E1, E2, beta, r);
}


//
// Calculates the energy due to the surrounding cells
//

int energy_calculation(int *lattice, int *coreBoundaries, int *nodeBoundaries, int location, int x_Max, int y_Max, int z_Max, int J)
{
    int energy = 0;
    int columnOffset = x_Max/gNumOfthreads * omp_get_thread_num();

    int x,y,z;
    coordinates_from_linear_index(location,x_Max/gNumOfthreads, y_Max/gNumOfNodes, &x, &y, &z);
    int offsetLocation = linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset,  y,  z);
    //
    // Calculate energy of cell above
    //
    (y == 0) ? (energy += J*lattice[offsetLocation  ] * nodeBoundaries[columnOffset + x +z*x_Max]*(-1) )
                : (energy += J*lattice[offsetLocation ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset,  y - 1,  z)]*(-1));

    //
    // Calculate energy of cell below
    //
    (y == (y_Max/gNumOfNodes-1)) ? (energy += J*lattice[offsetLocation  ] * nodeBoundaries[(columnOffset + x + z*x_Max)+x_Max*z_Max]*(-1))
                    : (energy += J*lattice[offsetLocation ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset,  y+1,  z)]*(-1));

    //
    // Calculate energy of cell left
    //

    (x == 0) ? (energy += J*lattice[offsetLocation  ] * coreBoundaries[y + z*y_Max/gNumOfNodes ]*(-1))
                    : (energy += J*lattice[offsetLocation ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset -1,  y,  z)]*(-1));
    //
    // Calculate energy of cell right
    //

    (x == (x_Max-1)) ? (energy += J*lattice[offsetLocation ] * coreBoundaries[y + z*y_Max/gNumOfNodes +(y_Max/gNumOfNodes)*z_Max]*(-1))
                    : (energy += J*lattice[offsetLocation ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset +1,  y,  z)]*(-1));
    if( z_Max > 1)
    {
        //
        // Calculate energy of cell front
        //
        (z == 0) ? (energy += J*lattice[offsetLocation  ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset ,  y,  z_Max -1)]*(-1))
                            : (energy += J*lattice[offsetLocation ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset ,  y,  z-1)]*(-1));
        //
        // Calculate energy of cell behind
        //
        (z == (z_Max-1)) ? (energy += J*lattice[offsetLocation ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset ,  y,  0)]*(-1))
                        : (energy += J*lattice[offsetLocation ] * lattice[linear_index_from_coordinates(x_Max,y_Max/gNumOfNodes,  x + columnOffset ,  y,  z+1)]*(-1));
    }
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

int linear_index_from_coordinates(int x_Max,int y_Max, int x, int y, int z)
{
    return x+x_Max*y+x_Max*y_Max*z;
}

void coordinates_from_linear_index(int location,int x_Max,int y_Max, int *x, int *y, int *z)
{
     *x = location%x_Max;
     location=(location - *x)/x_Max;
     *y = location%y_Max;
     location=(location-*y)/y_Max;
     *z = location;
}
