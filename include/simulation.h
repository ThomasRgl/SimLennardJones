#pragma once 

#include "system.h"
#include <stdbool.h>
#include <stdlib.h>

typedef struct{
    //
    bool log; 
    int verbose; 
    bool step_by_step; 
    int save_freq; 

    //
    // int n_local_particules;
    int max_ite;
    double dt;
    // double mass;
    int n_sym;

} simulation_config_t;


// simulation
int init_simulation( microscopic_system_t * sys, simulation_config_t config );
int start_simulation( microscopic_system_t * sys, const char * trajectory_filename, 
                        simulation_config_t config );

// lsj
double ljs( microscopic_system_t * sys, size_t N_sym);
double ljs_sum_forces( microscopic_system_t * sys);

// vv
int warning( microscopic_system_t * sys );
int correct_position( microscopic_system_t * sys );
double velocity_verlet_update( microscopic_system_t * sys, int n_sym, 
                              double dt, double mass);

// moments
int compute_cinetic_energy( microscopic_system_t * sys, double mass );
int initialize_moments( microscopic_system_t * sys );
int calibrate_moments( microscopic_system_t * sys );
int correct_moments( microscopic_system_t * sys );

// temperature
int Berendsen_thermostat( microscopic_system_t * sys );
int compute_temperature( microscopic_system_t * sys );
