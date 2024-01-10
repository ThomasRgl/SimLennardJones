#pragma once 

#include "system.h"
#include <stdbool.h>


// simulation
int start_simulation( microscopic_system_t * sys, const char * trajectory_filename, 
                        const size_t max_ite, const double dt, const int N_sym,
                        bool verbose, bool log );

// lsj
double ljs_potential( microscopic_system_t * sys, size_t N_sym);
int ljs_forces( microscopic_system_t * sys, size_t N_sym);
double ljs_sum_forces( microscopic_system_t * sys);

// vv
double velocity_verlet_update( microscopic_system_t * sys, int n_sym, double dt);
