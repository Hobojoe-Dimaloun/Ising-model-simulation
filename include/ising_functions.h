#ifndef _ISING_FUNCTIONS_H
#define _ISING_FUNCTIONS_H

#include <stdlib.h>
#include <gsl/gsl_rng.h>

void energy_comparison(int x_Max, int y_Max, int z_Max, gsl_rng *r, int *ising_lattice, int *ising_lattice_core_boundaries,int *ising_lattice_node_boundaries, double beta, int J,int chunksize);

void randLattice(int *lattice, int x, int y, int z);

int energy_calculation(int *lattice, int *coreBoundaries, int *nodeBoundaries, int location, int x_Max, int y_Max, int z_Max, int J);

int spin_flip_check(int E1, int E2, double beta, gsl_rng *r);

int linear_index_from_coordinates(int x_Max,int y_Max, int x, int y, int z);

void coordinates_from_linear_index(int location,int x_Max,int y_Max, int *x, int *y, int *z);

#endif // _ISING_FUNCTIONS_H
