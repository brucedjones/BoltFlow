#ifndef SOLVER_H
#define SOLVER_H

#include "data_types.h"

// CUDA KERNEL PROTOTYPES
void iterate_kernel (Lattice *lattice, Domain *domain_arrays, bool store_macros, int t, int x, int y, int z);

#endif