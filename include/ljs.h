#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "system.h"

double ljs_potential( microscopic_system_t * sys);

int ljs_forces( microscopic_system_t * sys);
double ljs_sum_forces( microscopic_system_t * sys);
