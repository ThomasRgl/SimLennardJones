#pragma once

#include "system.h"


int save_trajectory( FILE * fp, microscopic_system_t * sys, int ite);
size_t parse_data( const char * filename, double ** x, double ** y, double ** z );


