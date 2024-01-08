#include <stdio.h>
#include <string.h>

#include "ljs.h"
#include "system.h"

const double r = 3.0;
const double epsilon = 0.2;
const double L = 30;
const double Rcut = 10;



double ljs_potential( microscopic_system_t * sys, size_t N_sym){

    size_t N = sys->N_particules_total;

    // main box
    dim3_t sym_offset = {0,0,0};

    double U = 0;
    double u_tmp = 0;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            double r_2 = squared_distance(sys, i, j, sym_offset);
            if( r_2 < Rcut * Rcut ) {
                double rr_2 = ( r * r ) / r_2;
                double rr_6 =  rr_2 * rr_2 * rr_2;
                double rr_12 = rr_6 * rr_6;
                u_tmp +=  ( rr_12 - 2 * rr_6 );
            }
        }
    }
    U += 4 * u_tmp * epsilon ;

    u_tmp = 0;
    for (int l = 1; l < N_sym; l++) {
        dim3_t sym_offset = {0,0,0};
        sym_offset.x = ((l + 1) % 3 - 1) * L;
        sym_offset.y = ((l / 3 + 1) % 3 - 1) * L;
        sym_offset.z = ((l / 9 + 1) % 3 - 1) * L;

        // printf("%d %d %d\n", vecteurs[i][0], vecteurs[i][1], vecteurs[i][2]);
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                double r_2 = squared_distance(sys, i, j, sym_offset);
                if( r_2 < Rcut * Rcut ) {
                    double rr_2 = ( r * r ) / r_2;
                    double rr_6 =  rr_2 * rr_2 * rr_2;
                    double rr_12 = rr_6 * rr_6;

                    u_tmp += ( rr_12 - 2 * rr_6 );
                }
            }
        }
    }
    
    U += 4 * u_tmp * epsilon;
    return U;
}

int ljs_forces( microscopic_system_t * sys, size_t N_sym){


    size_t N = sys->N_particules_total;

    memset(sys->fx, 0, N * sizeof(double) );
    memset(sys->fy, 0, N * sizeof(double) );
    memset(sys->fz, 0, N * sizeof(double) );

    dim3_t sym_offset = {0,0,0};
    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            double r_2 = squared_distance(sys, i, j, sym_offset);
            if( r_2 < Rcut * Rcut ) {
                double rr_2 = ( r * r ) / r_2;
                double rr_6 =  rr_2 * rr_2 * rr_2;
                double rr_12 = rr_6 * rr_6;
                
                double fx = rr_2 * (rr_12 - rr_6) * (sys->px[i] - ( sys->px[j] + sym_offset.x ) );
                double fy = rr_2 * (rr_12 - rr_6) * (sys->py[i] - ( sys->py[j] + sym_offset.y ) );
                double fz = rr_2 * (rr_12 - rr_6) * (sys->pz[i] - ( sys->pz[j] + sym_offset.z ) );
                // printf(" -- { %f } \n", - rr_2 * (rr_12 - rr_6) * (sys->px[i] - ( sys->px[j] + sym_offset.x ) ) );

                sys->fx[i] += fx;
                sys->fy[i] += fy;
                sys->fz[i] += fz;

                sys->fx[j] -= fx;
                sys->fy[j] -= fy;
                sys->fz[j] -= fz;
            }
        }
    }

    for (int l = 1; l < N_sym; l++) {
        dim3_t sym_offset = {0,0,0};
        sym_offset.x = ((l + 1) % 3 - 1) * L ;
        sym_offset.y = ((l / 3 + 1) % 3 - 1) * L ;
        sym_offset.z = ((l / 9 + 1) % 3 - 1) * L;

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {

                double r_2 = squared_distance(sys, i, j, sym_offset);
                if( r_2 <= Rcut * Rcut ) {
                    double rr_2 = ( r * r ) / r_2;
                    double rr_6 =  rr_2 * rr_2 * rr_2;
                    double rr_12 = rr_6 * rr_6;
                    
                    double fx = rr_2 * (rr_12 - rr_6) * (sys->px[i] - ( sys->px[j] + sym_offset.x ) );
                    double fy = rr_2 * (rr_12 - rr_6) * (sys->py[i] - ( sys->py[j] + sym_offset.y ) );
                    double fz = rr_2 * (rr_12 - rr_6) * (sys->pz[i] - ( sys->pz[j] + sym_offset.z ) );

                    sys->fx[i] += fx;
                    sys->fy[i] += fy;
                    sys->fz[i] += fz;

                    // sys->fx[j] += fx;
                    // sys->fy[j] += fy;
                    // sys->fz[j] += fz;
                }
            }
        }
    }

    for (size_t i = 0; i < N; i++) {
        sys->fx[i] *= - 48 * epsilon;
        sys->fy[i] *= - 48 * epsilon;
        sys->fz[i] *= - 48 * epsilon;
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

    printf(" ljs sum forces : { %g, %g, %g} \n", sum_x, sum_y, sum_z );

    return sum_x + sum_y + sum_z;
}



// double ljs_potential( microscopic_system_t * sys, size_t N_sym){
//
//     double u_tmp = 0;
//     size_t N = sys->N_particules_total;
//
//     for (int l = 0; l < N_sym; l++) {
//         dim3_t sym_offset = {0,0,0};
//         sym_offset.x = ((l + 1) % 3 - 1) * L;
//         sym_offset.y = ((l / 3 + 1) % 3 - 1) * L;
//         sym_offset.z = ((l / 9 + 1) % 3 - 1) * L;
//
//         // printf("%d %d %d\n", vecteurs[i][0], vecteurs[i][1], vecteurs[i][2]);
//         for (size_t i = 0; i < N; i++) {
//             for (size_t j = i + 1; j < N; j++) {
//                 double r_2 = squared_distance(sys, i, j, sym_offset);
//                 if( r_2 < Rcut * Rcut ) {
//                     double rr_2 = ( r * r ) / r_2;
//                     double rr_6 =  rr_2 * rr_2 * rr_2;
//                     double rr_12 = rr_6 * rr_6;
//                     u_tmp += epsilon * ( rr_12 - 2 * rr_6 );
//                 }
//             }
//         }
//     }
//     
//     double U = 4 * u_tmp;
//     return U;
// }


// int ljs_forces( microscopic_system_t * sys, size_t N_sym){
//
//
//     size_t N = sys->N_particules_total;
//
//     memset(sys->fx, 0, N * sizeof(double) );
//     memset(sys->fy, 0, N * sizeof(double) );
//     memset(sys->fz, 0, N * sizeof(double) );
//
//     for (int l = 0; l < N_sym; l++) {
//         dim3_t sym_offset = {0,0,0};
//         sym_offset.x = ((l + 1) % 3 - 1) * L ;
//         sym_offset.y = ((l / 3 + 1) % 3 - 1) * L ;
//         sym_offset.z = ((l / 9 + 1) % 3 - 1) * L;
//         // printf(" { %f %f %f } \n", sym_offset.z, sym_offset.y, sym_offset.z );
//
//         for (size_t i = 0; i < N; i++) {
//             for (size_t j = i + 1; j < N; j++) {
//
//                 double r_2 = squared_distance(sys, i, j, sym_offset);
//                 if( r_2 < Rcut * Rcut ) {
//                     double rr_2 = ( r * r ) / r_2;
//                     double rr_6 =  rr_2 * rr_2 * rr_2;
//                     double rr_12 = rr_6 * rr_6;
//                     
//                     double fx = 48 * epsilon * rr_2 * (rr_12 - rr_6) * (sys->px[i] - ( sys->px[j] + sym_offset.x ) );
//                     double fy = 48 * epsilon * rr_2 * (rr_12 - rr_6) * (sys->py[i] - ( sys->py[j] + sym_offset.y ) );
//                     double fz = 48 * epsilon * rr_2 * (rr_12 - rr_6) * (sys->pz[i] - ( sys->pz[j] + sym_offset.z ) );
//
//                     sys->fx[i] -= fx;
//                     sys->fy[i] -= fy;
//                     sys->fz[i] -= fz;
//
//                     sys->fx[j] += fx;
//                     sys->fy[j] += fy;
//                     sys->fz[j] += fz;
//                 }
//             }
//         }
//     }
//
//     return 0;
// }
