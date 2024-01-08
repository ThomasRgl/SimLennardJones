#pragma once 

#include "system.h"

double ljs_potential( microscopic_system_t * sys, size_t N_sym);

int ljs_forces( microscopic_system_t * sys, size_t N_sym);
double ljs_sum_forces( microscopic_system_t * sys);
