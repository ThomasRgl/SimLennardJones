#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

#include "simulation.h"
#include "system.h"
#include "io.h"

const double r_star = 3.0;
const double epsilon = 0.2;
const double L = 33;
const double Rcut = 10;
const double T0 = 300;
const double CONSTANTE_R = 0.00199;
const double conversion_force =0.0001*4.186;
const double gamma_ = 0.01;


int init_simulation( microscopic_system_t * sys, simulation_config_t config ){

    initialize_moments( sys );

    calibrate_moments( sys );
    correct_moments( sys );
    calibrate_moments( sys );

    return 0;

}



int start_simulation( microscopic_system_t * sys, const char * trajectory_filename, 
                        simulation_config_t config ){

    const int save_freq = config.save_freq;
    const int verbose = config.verbose;
    const bool log = config.log;
    const bool step_by_step = config.step_by_step;

    const int max_ite = config.max_ite;
    const double dt = config.dt;
    const int n_sym = config.n_sym;

    const int mass = sys->mass;

    FILE * fp_traj = fopen(trajectory_filename, "w+");
    FILE * fp_log = fopen("../logs/log.csv", "w+");

    init_log( fp_log );

    // initial position

    // si
    // correct_position(sys);
    ljs( sys, n_sym );
    for (int i = 0; i <= max_ite; i++) {
    
        printf("%d ...\n", i);

         if( i % save_freq == 0){
            save_trajectory( fp_traj, sys, i);
        }
        
        compute_cinetic_energy(sys, mass);
        compute_temperature(sys);
        ljs_sum_forces(sys);

        // ljs_sum_speed(sys);

        // warning(sys);

        if( verbose >= 1 ){
            printf("(%d) ljs potential : %f \n", n_sym, sys->potential_energy );
            printf("(%d) ljs sum forces : %g \n", n_sym, sys->sum_forces );
            printf("(%d) cinetic energy : %g \n", n_sym, sys->cinetic_energy );
            printf("(%d) temperature : %g \n", n_sym, sys->temperature );
        }

        if( verbose >= 3 )
            print_system(*sys);

        add_log(fp_log, sys); 

        if( step_by_step )
            getchar();

        // update system
        velocity_verlet_update( sys, n_sym, dt, mass);
        compute_cinetic_energy(sys, mass);
        // if( i < 400 )
            Berendsen_thermostat(sys);
        correct_position(sys);


    }

    fclose(fp_traj);
    fclose(fp_log);



    return 0;
}




// LJS

double ljs( microscopic_system_t * sys, size_t N_sym){

    size_t N = sys->N_particules_local;

    memset(sys->fx, 0, N * sizeof(double) );
    memset(sys->fy, 0, N * sizeof(double) );
    memset(sys->fz, 0, N * sizeof(double) );

    double U = 0;

    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {

            // printf("Central box : %3lu %3lu \n", i, j );

            double distance_x = sys->x[i] - sys->x[j] + 0.f;
            double distance_y = sys->y[i] - sys->y[j] + 0.f;
            double distance_z = sys->z[i] - sys->z[j] + 0.f;
            double r_2 = distance_x * distance_x + distance_y * distance_y + distance_z * distance_z;

            if( r_2 < Rcut * Rcut ) { //TODO
                double rr_2 = ( r_star * r_star ) / r_2;
                
                double rr_6 =  rr_2 * rr_2 * rr_2;
                double rr_8 =  rr_6 * rr_2;
                double rr_12 = rr_6 * rr_6;
                double rr_14 = rr_12 * rr_2;
                
                double fx = (rr_14 - rr_8) * distance_x;
                double fy = (rr_14 - rr_8) * distance_y;
                double fz = (rr_14 - rr_8) * distance_z;

                sys->fx[i] += fx;
                sys->fy[i] += fy;
                sys->fz[i] += fz;

                sys->fx[j] -= fx;
                sys->fy[j] -= fy;
                sys->fz[j] -= fz;

                // if( i == 0 && j == 1){
                //     printf("U 0 - 1 = %f\n", 4 * epsilon * ( rr_12 - 2 * rr_6 ));
                //     printf("r2 %f\n", r_2);
                //     printf("r_star %f\n", r_star);
                //     printf("rr_2 %f\n", rr_2);
                //     printf("rr_6 %g\n", rr_6);
                //     printf("rr_12 %g\n", rr_12);
                // }        

                U += ( rr_12 - 2 * rr_6 );
            }

            
        }
        // if( i == 0 ){
        //     printf("U 0 = %f\n", 4 * epsilon * U);
        // }
        
    }

    #pragma omp parallel for
    for (int l = 1; l < N_sym; l++) {

        double offset_x = ((l + 1) % 3 - 1) * L ;
        double offset_y = ((l / 3 + 1) % 3 - 1) * L ;
        double offset_z = ((l / 9 + 1) % 3 - 1) * L;

        // printf(" ========================= \n" );
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {

                // printf("Periodic : %3lu %3lu \n", i, j );

                double distance_x = sys->x[i] - sys->x[j] + offset_x;
                double distance_y = sys->y[i] - sys->y[j] + offset_y;
                double distance_z = sys->z[i] - sys->z[j] + offset_z;
                double r_2 = distance_x * distance_x + distance_y * distance_y + distance_z * distance_z;

                // if( r_2 < 1.5)
                //     r_2 = 1.5; ??

                // double r_2 = squared_distance(sys, i, j, sym_offset);
                if( r_2 <= Rcut * Rcut ) {
                    double rr_2 = ( r_star * r_star ) / r_2;
 
                    double rr_6 =  rr_2 * rr_2 * rr_2;
                    double rr_8 =  rr_6 * rr_2;
                    double rr_12 = rr_6 * rr_6;
                    double rr_14 = rr_12 * rr_2;
                    
                    double fx = (rr_14 - rr_8) * distance_x;
                    double fy = (rr_14 - rr_8) * distance_y;
                    double fz = (rr_14 - rr_8) * distance_z;
                    
                    #pragma omp critical
                    {
                    sys->fx[i] += fx;
                    sys->fy[i] += fy;
                    sys->fz[i] += fz;

                    U += ( rr_12 - 2 * rr_6 );
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < N; i++) {
        sys->fx[i] *= - 48 * epsilon;
        sys->fy[i] *= - 48 * epsilon;
        sys->fz[i] *= - 48 * epsilon;
    }

    U *= 4 * epsilon;

    sys->potential_energy = U;

    return U;
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

    // printf(" ljs sum forces : { %g, %g, %g} \n", sum_x, sum_y, sum_z );

    sys->sum_forces = sum_x + sum_y + sum_z ;
    return sum_x + sum_y + sum_z;
}

// 
int warning( microscopic_system_t * sys ){

    size_t N = sys->N_particules_local;

    // printf("=================WARNING=================\n");

    #pragma omp parallel for
    for (int l = 1; l < 27; l++) 
    // for (int l = 1; l < 1; l++) 
    {

        double offset_x = ((l + 1) % 3 - 1) * L ;
        double offset_y = ((l / 3 + 1) % 3 - 1) * L ;
        double offset_z = ((l / 9 + 1) % 3 - 1) * L;

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                if( i == j )
                    continue;
                double distance_x = sys->x[i] - sys->x[j] + offset_x;
                double distance_y = sys->y[i] - sys->y[j] + offset_y;
                double distance_z = sys->z[i] - sys->z[j] + offset_z;
                double r_2 = distance_x * distance_x + distance_y * distance_y + distance_z * distance_z;

                // double r_2 = squared_distance(sys, i, j, sym_offset);

                if(sqrt(r_2) < 2){printf("/!\\ dist %lu to %lu, %f\n", i, j, sqrt(r_2) );}
            }
        }
    }

    // printf("=========================================\n\n");

    return 0;

}

int correct_position( microscopic_system_t * sys ){

    size_t N = sys->N_particules_local;
    double pos = 0;
    for (size_t i = 0; i < N; i++) {
        // Correction pour la coordonnée x
        double old_x = sys->x[i];
        if (sys->x[i] > L/2) {
            sys->x[i] = fmod(sys->x[i] - L/2, L) - L/2;
        } else if (sys->x[i] < -L/2) {
            sys->x[i] = fmod(sys->x[i] + L/2, L) + L/2;
        }

        // Correction pour la coordonnée y
        double old_y = sys->y[i];
        if (sys->y[i] > L/2) {
            sys->y[i] = fmod(sys->y[i] - L/2, L) - L/2;
        } else if (sys->y[i] < -L/2) {
            sys->y[i] = fmod(sys->y[i] + L/2, L) + L/2;
        }

        // Correction pour la coordonnée z
        double old_z = sys->z[i];
        if (sys->z[i] > L/2) {
            sys->z[i] = fmod(sys->z[i] - L/2, L) - L/2;
        } else if (sys->z[i] < -L/2) {
            sys->z[i] = fmod(sys->z[i] + L/2, L) + L/2;
        }
    
    }
    
    return 0;

}


int compute_temperature( microscopic_system_t * sys ){
    const double ndl = ( 3 * sys->N_particules_local - 3 );
    sys->temperature = ( 1 / ( ndl * CONSTANTE_R ) ) * sys->cinetic_energy;

    return 0;
}

int compute_cinetic_energy( microscopic_system_t * sys, double mass ){
    
    double cinetic_energy = 0;
    for (int i = 0; i < sys->N_particules_local ; i++) {
        
        cinetic_energy += ( sys->px[i] * sys->px[i] + 
                            sys->py[i] * sys->py[i] + 
                            sys->pz[i] * sys->pz[i]) / mass;
    }

    cinetic_energy = cinetic_energy * (1 / ( 2 * conversion_force ));
    sys->cinetic_energy = cinetic_energy;
     
    return 0;
}

double velocity_verlet_update( microscopic_system_t * sys, int n_sym, double dt, double mass){
    
    size_t N = sys->N_particules_local;
    // double dt = 0.001f;
    const double inv_mass = 1 / mass;

    for (int i = 0; i < N ; i++) {
        sys->vx[i] = sys->px[i] * mass;
        sys->vy[i] = sys->py[i] * mass;
        sys->vz[i] = sys->pz[i] * mass;
    }

    for (int i = 0; i < N ; i++) {
        sys->x[i] -= ( sys->vx[i] * dt ) + ( sys->fx[i] * inv_mass * dt * dt * 0.5 );
        sys->y[i] -= ( sys->vy[i] * dt ) + ( sys->fy[i] * inv_mass * dt * dt * 0.5 );
        sys->z[i] -= ( sys->vz[i] * dt ) + ( sys->fz[i] * inv_mass * dt * dt * 0.5 );

    }

    for (int i = 0; i < N ; i++) {
        sys->vx[i] += ( sys->fx[i] * inv_mass * dt * 0.5 );
        sys->vy[i] += ( sys->fy[i] * inv_mass * dt * 0.5 );
        sys->vz[i] += ( sys->fz[i] * inv_mass * dt * 0.5 );
    }

    ljs( sys, n_sym );

    for (int i = 0; i < N ; i++) {
        sys->vx[i] += ( sys->fx[i] * inv_mass * dt * 0.5 );
        sys->vy[i] += ( sys->fy[i] * inv_mass * dt * 0.5 );
        sys->vz[i] += ( sys->fz[i] * inv_mass * dt * 0.5 );
    }
    
    for (int i = 0; i < N ; i++) {
        sys->px[i] = sys->vx[i] / mass;
        sys->py[i] = sys->vy[i] / mass;
        sys->pz[i] = sys->vz[i] / mass;
    }

    return 0;
}

int initialize_moments( microscopic_system_t * sys ){

    size_t N = sys->N_particules_local;

    // init
    for (int i = 0; i < N ; i++) {
        double s_x = (double)rand() / RAND_MAX;
        double s_y = (double)rand() / RAND_MAX;
        double s_z = (double)rand() / RAND_MAX;

        double c_x = (double)rand() / RAND_MAX;
        double c_y = (double)rand() / RAND_MAX;
        double c_z = (double)rand() / RAND_MAX;

        sys->px[i] = copysign(1.0, 0.5 - s_x) * c_x;
        sys->py[i] = copysign(1.0, 0.5 - s_y) * c_y;
        sys->pz[i] = copysign(1.0, 0.5 - s_z) * c_z;

    }

    return 0;
}

int calibrate_moments( microscopic_system_t * sys ){

    size_t N = sys->N_particules_local;

    //recalibration
    const double ndl = ( 3 * N - 3 );
    compute_cinetic_energy(sys, sys->mass);
    double cinetic_energy = sys->cinetic_energy;
    const double RAPPORT = ndl * CONSTANTE_R * T0 / cinetic_energy;

    for (int i = 0; i < N ; i++) {
        sys->px[i] *= RAPPORT;
        sys->py[i] *= RAPPORT;
        sys->pz[i] *= RAPPORT;
    }

    for (int i = 0; i < N ; i++) {
        sys->vx[i] = sys->px[i] * sys->mass;
        sys->vy[i] = sys->py[i] * sys->mass;
        sys->vz[i] = sys->pz[i] * sys->mass;
    }

    return 0;
}

int correct_moments( microscopic_system_t * sys ){

    size_t N = sys->N_particules_local;

    // init
    double Px = 0;
    double Py = 0;
    double Pz = 0;

    for (int i = 0; i < N ; i++) {
        Px += sys->px[i];
        Py += sys->py[i];
        Pz += sys->pz[i];
    }

    Px = Px / N;
    Py = Py / N;
    Pz = Pz / N;

    for (int i = 0; i < N ; i++) {
        sys->px[i] -= Px;
        sys->py[i] -= Py;
        sys->pz[i] -= Pz;
    }

    return 0;
}


int Berendsen_thermostat( microscopic_system_t * sys ){
    
    compute_temperature(sys);

    const int N = sys->N_particules_local;
    const double temp = sys->temperature;

    for (int i = 0; i < N ; i++) {
        sys->px[i] += ( T0/temp - 1 ) * gamma_ * sys->px[i];
        sys->py[i] += ( T0/temp - 1 ) * gamma_ * sys->py[i];
        sys->pz[i] += ( T0/temp - 1 ) * gamma_ * sys->pz[i];
    }

    return 0;
}
