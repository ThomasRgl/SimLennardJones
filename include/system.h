#pragma once 
#include <stdio.h>

typedef struct {
    double x;
    double y;
    double z;
} dim3_t;

typedef struct {
    // particule_t particules;
    double * x;
    double * y;
    double * z;

    double * vx;
    double * vy;
    double * vz;

    double * fx;
    double * fy;
    double * fz;

    double * px;
    double * py;
    double * pz;

    double dimx;
    double dimy;
    double dimz;

    double cinetic_energy;
    double potential_energy;
    double temperature;
    double sum_forces;

    double mass;

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
int create_system( microscopic_system_t * sys, const char * filename );

