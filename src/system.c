#include <iso646.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "system.h"
#include "io.h"


double distance( microscopic_system_t * sys, const size_t i, 
                        const size_t j, const dim3_t offset){

    double dist_x = ( sys->x[i] - ( sys->x[j] + offset.x ) );
    double dist_y = ( sys->y[i] - ( sys->y[j] + offset.y ) );
    double dist_z = ( sys->z[i] - ( sys->z[j] + offset.z ) );

    return sqrt( dist_x * dist_x + dist_y * dist_y + dist_z * dist_z );
}

double squared_distance( microscopic_system_t * sys, const size_t i, 
                        const size_t j, const dim3_t offset){
    double dist_x = ( sys->x[i] - ( sys->x[j] + offset.x ) );
    double dist_y = ( sys->y[i] - ( sys->y[j] + offset.y ) );
    double dist_z = ( sys->z[i] - ( sys->z[j] + offset.z ) );

    // if(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z < 0.05){
    //     printf(" %lu %lu  = %f\n", i, j, sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z));
    // }
    return dist_x * dist_x + dist_y * dist_y + dist_z * dist_z ;
}



int print_system( microscopic_system_t sys ){

    printf("----------------------------------------\n");
    printf("x\ty\tz\n");
    for (size_t i = 0; i < sys.N_particules_local; i++) {
        printf("%lu\t%lf\t%lf\t%lf\n",i, sys.x[i], sys.y[i], sys.z[i]);
    }
    
    printf("\n\nvx\tvy\tvz\n");
    for (size_t i = 0; i < sys.N_particules_local; i++) {
        printf("%lu\t%lf\t%lf\t%lf\n",i, sys.vx[i], sys.vy[i], sys.vz[i]);
    }

    printf("\n\nfx\tfy\tfz\n");
    for (size_t i = 0; i < sys.N_particules_local; i++) {
        printf("%lu\t%lf\t%lf\t%lf\n",i, sys.fx[i], sys.fy[i], sys.fz[i]);
    }
    printf("----------------------------------------\n");
    return 0;
}

int allocate_system( microscopic_system_t * sys, size_t indice ){
 
    sys->N_particules_total = indice;
    sys->N_particules_local = indice;

    // aligned alloc is chouining
    sys->x = aligned_alloc(32, indice * sizeof(double) );
    sys->y = aligned_alloc(32, indice * sizeof(double) );
    sys->z = aligned_alloc(32, indice * sizeof(double) );

    sys->fx = aligned_alloc(32, indice * sizeof(double) );
    sys->fy = aligned_alloc(32, indice * sizeof(double) );
    sys->fz = aligned_alloc(32, indice * sizeof(double) );

    sys->vx = aligned_alloc(32, indice * sizeof(double) );
    sys->vy = aligned_alloc(32, indice * sizeof(double) );
    sys->vz = aligned_alloc(32, indice * sizeof(double) );

    sys->px = aligned_alloc(32, indice * sizeof(double) );
    sys->py = aligned_alloc(32, indice * sizeof(double) );
    sys->pz = aligned_alloc(32, indice * sizeof(double) );
 
    memset(sys->fx, 0, indice * sizeof(double) );
    memset(sys->fy, 0, indice * sizeof(double) );
    memset(sys->fz, 0, indice * sizeof(double) );
   
    memset(sys->vx, 0, indice * sizeof(double) );
    memset(sys->vy, 0, indice * sizeof(double) );
    memset(sys->vz, 0, indice * sizeof(double) );

    memset(sys->px, 0, indice * sizeof(double) );
    memset(sys->py, 0, indice * sizeof(double) );
    memset(sys->pz, 0, indice * sizeof(double) );

      
    return 0;
}

int free_system( microscopic_system_t sys ){
    free(sys.x);
    free(sys.y);
    free(sys.z);

    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    free(sys.vx);
    free(sys.vy);
    free(sys.vz);

    free(sys.px);
    free(sys.py);
    free(sys.pz);


    return 0;
}


int create_system( microscopic_system_t * sys, const char * filename ){
    
    double * x; 
    double * y;
    double * z;


    size_t indice = parse_data(filename, &x, &y, &z);
    allocate_system(sys, indice);
    
    for (size_t i = 0; i < indice; i++) {
        sys->x[i] = x[i];
        sys->y[i] = y[i];
        sys->z[i] = z[i];
    }

    sys->dimx = 30;
    sys->dimy = 30;
    sys->dimz = 30;

    

    free(x);
    free(y);
    free(z);
    return 0;
}
