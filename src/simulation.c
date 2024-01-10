#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "simulation.h"
#include "system.h"
#include "io.h"

const double r = 3.0;
const double epsilon = 0.2;
const double L = 30;
const double Rcut = 10;



int start_simulation( microscopic_system_t * sys, const char * trajectory_filename, 
                        const size_t max_ite, const double dt, const int n_sym,
                        bool verbose, bool log ){

    FILE * trajectory_file = fopen(trajectory_filename, "w+");

    save_trajectory( trajectory_file, sys, 0);
    for (int i = 1; i <= max_ite; i++) {
    
        // printf("%d - (%d) ljs potential : %f \n", i, n_sym, ljs_potential(&sys, n_sym) );
        // ljs_forces( &sys, n_sym );
        // double sum = ljs_sum_forces( &sys );
        // printf("%d - (%d) ljs sum forces : %g \n", i, n_sym, sum );
        // print_system(sys);
        
        printf("%d ...\n", i);
        velocity_verlet_update( sys, n_sym, dt);
        save_trajectory( trajectory_file, sys, 2);
        
        // print_system(sys);
        // printf("(%d) ljs potential : %f \n", n_sym, ljs_potential(&sys, n_sym) );
        // ljs_forces( &sys, n_sym );
        // double sum = ljs_sum_forces( &sys );
        // printf("(%d) ljs sum forces : %g \n", n_sym, sum );

    }

    fclose(trajectory_file);



    return 0;
}

double ljs_potential( microscopic_system_t * sys, size_t N_sym){

    size_t N = sys->N_particules_local;

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


    size_t N = sys->N_particules_local;

    memset(sys->fx, 0, N * sizeof(double) );
    memset(sys->fy, 0, N * sizeof(double) );
    memset(sys->fz, 0, N * sizeof(double) );

    dim3_t sym_offset = {0,0,0};
    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            double r_2 = squared_distance(sys, i, j, sym_offset);
            // if(i == 0){printf("dist 26 to %lu, %f\n", j, sqrt(r_2) );}
            if(sqrt(r_2) < 1){printf("dist %lu to %lu, %f\n", i, j, sqrt(r_2) );}
            if( r_2 < Rcut * Rcut ) {
                double rr_2 = ( r * r ) / r_2;
                
                // if(i == 0){printf("rr_2, %f\n", rr_2 );}
                double rr_6 =  rr_2 * rr_2 * rr_2;
                // if(i == 0){printf("rr_6, %f\n", rr_6 );}
                double rr_12 = rr_6 * rr_6;
                // if(i == 0){printf("rr_12, %f\n", rr_12 );}
                
                double fx = rr_2 * (rr_12 - rr_6) * (sys->x[i] - ( sys->x[j] + sym_offset.x ) );
                double fy = rr_2 * (rr_12 - rr_6) * (sys->y[i] - ( sys->y[j] + sym_offset.y ) );
                double fz = rr_2 * (rr_12 - rr_6) * (sys->z[i] - ( sys->z[j] + sym_offset.z ) );
                // if(i == 0){printf("fx, %f\n", fx );}
                // if(i == 0){printf("fy, %f\n", fy );}
                // if(i == 0){printf("fz, %f\n", fz );}
                // printf(" -- { %f } \n", - rr_2 * (rr_12 - rr_6) * (sys->x[i] - ( sys->x[j] + sym_offset.x ) ) );

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
                    
                    double fx = rr_2 * (rr_12 - rr_6) * (sys->x[i] - ( sys->x[j] + sym_offset.x ) );
                    double fy = rr_2 * (rr_12 - rr_6) * (sys->y[i] - ( sys->y[j] + sym_offset.y ) );
                    double fz = rr_2 * (rr_12 - rr_6) * (sys->z[i] - ( sys->z[j] + sym_offset.z ) );

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

    size_t N = sys->N_particules_local;

    for (size_t i = 0; i < N; i++) {
            sum_x += sys->fx[i] ;
            sum_y += sys->fy[i] ;
            sum_z += sys->fz[i] ;
    }

    printf(" ljs sum forces : { %g, %g, %g} \n", sum_x, sum_y, sum_z );

    return sum_x + sum_y + sum_z;
}


double velocity_verlet_update( microscopic_system_t * sys, int n_sym, double dt){
    
    size_t N = sys->N_particules_local;
    // double dt = 0.001f;
    const double mass = 1;
    const double inv_mass = 1 / mass;
    const double conversion_force =0.0001*4.186;


    for (int i = 0; i < N ; i++) {
        sys->x[i] -= ( sys->vx[i] * dt ) + ( conversion_force * sys->fx[i] * inv_mass * dt * dt * 0.5 );
        sys->y[i] -= ( sys->vy[i] * dt ) + ( conversion_force * sys->fy[i] * inv_mass * dt * dt * 0.5 );
        sys->z[i] -= ( sys->vz[i] * dt ) + ( conversion_force * sys->fz[i] * inv_mass * dt * dt * 0.5 );
    }

    for (int i = 0; i < N ; i++) {
        sys->vx[i] += ( conversion_force * sys->fx[i] * inv_mass * dt * 0.5 );
        sys->vy[i] += ( conversion_force * sys->fy[i] * inv_mass * dt * 0.5 );
        sys->vz[i] += ( conversion_force * sys->fz[i] * inv_mass * dt * 0.5 );
    }

    ljs_forces( sys, n_sym );

    for (int i = 0; i < N ; i++) {
        sys->vx[i] += ( conversion_force * sys->fx[i] * inv_mass * dt * 0.5 );
        sys->vy[i] += ( conversion_force * sys->fy[i] * inv_mass * dt * 0.5 );
        sys->vz[i] += ( conversion_force * sys->fz[i] * inv_mass * dt * 0.5 );
    }

    return 0;
}

