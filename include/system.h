#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// typedef struct {
//     double x;
//     double y;
//     double z;
//     double buffer;
// } particule_t;

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


int print_system( microscopic_system_t sys );
int allocate_system( microscopic_system_t * sys, size_t indice );
int free_system( microscopic_system_t sys );
int parse_system_data( microscopic_system_t * sys, char * filename );
