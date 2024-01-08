#include <stdio.h>


#pragma once 

typedef struct {
    double x;
    double y;
    double z;
} dim3_t;

typedef struct {
    // particule_t particules;
    double * px;
    double * py;
    double * pz;

    double * fx;
    double * fy;
    double * fz;

    size_t N_particules_local;
    size_t N_particules_total;

} microscopic_system_t;


#define MAX_LINES 10000

// particules
double distance( microscopic_system_t * sys, const size_t i,
                        const size_t j, const dim3_t offset);

double squared_distance( microscopic_system_t * sys, const size_t i, 
                        const size_t j, const dim3_t offset);

//systems
int print_system( microscopic_system_t sys );
int allocate_system( microscopic_system_t * sys, size_t indice );
int free_system( microscopic_system_t sys );
int parse_system_data( microscopic_system_t * sys, char * filename );

