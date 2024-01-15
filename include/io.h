#pragma once

#include "system.h"


int save_trajectory( FILE * fp, microscopic_system_t * sys, int ite);
size_t parse_data( const char * filename, double ** x, double ** y, double ** z );

int init_log( FILE * fp_log );
int add_log( FILE * fp_log, microscopic_system_t * sys );

