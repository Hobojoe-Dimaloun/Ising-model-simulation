#ifndef _ISING_FUNCTIONS_H
#define _ISING_FUNCTIONS_H

#include <stdlib.h>
#include <gsl/gsl_rng.h>

void randLattice(int *lattice, int x, int y, int z);

int energy_calculation(int *lattice, int location, int x, int y, int J);

int spin_flip_check(int E1, int E2, double beta, gsl_rng *r);

#endif // _ISING_FUNCTIONS_H
