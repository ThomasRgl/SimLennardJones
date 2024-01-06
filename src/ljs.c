#include "ljs.h"
#include <stdio.h>
#include <string.h>
#include <math.h>




double ljs_potential( microscopic_system_t * sys){
    const double r = 3.0;
    const double epsilon = 0.2;

    double tmp = 0;
    size_t N = sys->N_particules_total;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            double r_2 = squared_distance(sys, i, j);
            double rr_2 = ( r * r ) / r_2;
            double rr_6 =  rr_2 * rr_2 * rr_2;
            double rr_12 = rr_6 * rr_6;
            tmp += epsilon * ( rr_12 - 2 * rr_6 );
        }
    }
    
    double U = 4 * tmp;
    return U;
}

int ljs_forces( microscopic_system_t * sys){

    const double r = 3.0;
    const double epsilon = 0.2;

    size_t N = sys->N_particules_total;

    memset(sys->fx, 0, N * sizeof(double) );
    memset(sys->fy, 0, N * sizeof(double) );
    memset(sys->fz, 0, N * sizeof(double) );

    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {

            double r_2 = squared_distance(sys, i, j);
            double rr_2 = ( r * r ) / r_2;
            double rr_6 =  rr_2 * rr_2 * rr_2;
            double rr_12 = rr_6 * rr_6;
            
            double fx = 48 * epsilon * rr_2 * (rr_12 - rr_6) * (sys->px[i] - sys->px[j]);
            double fy = 48 * epsilon * rr_2 * (rr_12 - rr_6) * (sys->py[i] - sys->py[j]);
            double fz = 48 * epsilon * rr_2 * (rr_12 - rr_6) * (sys->pz[i] - sys->pz[j]);

            sys->fx[i] -= fx;
            sys->fy[i] -= fy;
            sys->fz[i] -= fz;

            sys->fx[j] += fx;
            sys->fy[j] += fy;
            sys->fz[j] += fz;
        }
    }

    return 0;
}

double ljs_sum_forces( microscopic_system_t * sys){

    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;

    size_t N = sys->N_particules_total;

    for (size_t i = 0; i < N; i++) {
            sum_x += sys->fx[i] ;
            sum_y += sys->fy[i] ;
            sum_z += sys->fz[i] ;
    }

    return sum_x + sum_y + sum_z;
}

