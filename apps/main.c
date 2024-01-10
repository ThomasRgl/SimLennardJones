#include "system.h"
#include "simulation.h"
#include "argtable3.h"

#include <stdio.h>


const int default_iteration = 100;
const int default_dt = 1;
const double default_mass = 18;
const double default_periodicity = 27;


int main( int argc, char * argv[] ) {

    struct arg_file* arg_datafile = arg_file1(NULL, NULL, "DATA FILE", "A file containing particules positions");
    struct arg_int* arg_n = arg_int0("n", "local_particule", "INT", "Number of particule to simulate. Defaults to the number of particule in the datafile");
    struct arg_int* arg_i = arg_int0("i", "iteration", "INT", "Number of iteration to run in the simulation. Defaults to 100");
    struct arg_dbl* arg_t = arg_dbl0("t", "dt", "FLOAT", "delta time. Default to 1 [fs]");
    struct arg_dbl* arg_m = arg_dbl0("m", "mass", "FLOAT", "mass of one particule. Default to 18 ");
    struct arg_int* arg_p = arg_int0("p", "periodicity", "INT", "number of virtual boxes. Default to 27 ");
    struct arg_lit* verbose = arg_lit0("v", "verbose", "If present, additional information will be printed to stdout");
    struct arg_lit* help = arg_lit0("h", "help", "Display this help message");
    struct arg_lit* log = arg_lit0("l", "log", "print simulation information in a logfile");
    struct arg_end* end = arg_end(20);
 

    void* argtable[] = { arg_n, arg_i, arg_t, arg_m, arg_p, arg_datafile,
                        verbose, help, log, end};

    const int nb_opts = sizeof argtable / sizeof argtable[0];

    if (arg_nullcheck(argtable) != 0) {
        fprintf(stderr, "Argument parsing failed\n");
        arg_freetable(argtable, nb_opts);
        return 1;
    }

    arg_i->ival[0] = default_iteration;
    arg_t->dval[0] = default_dt;
    arg_m->dval[0] = default_mass;
    arg_p->ival[0] = default_periodicity;

    int nb_errors = arg_parse(argc, argv, argtable);

    if (help->count > 0) {
        printf("Usage: %s ", argv[0]);
        arg_print_syntax(stdout, argtable, "\n\n");
        arg_print_glossary(stdout, argtable, "  %-30s %s\n");
        arg_freetable(argtable, nb_opts);
        return 0;
    }

    if (nb_errors > 0) {
        arg_print_errors(stderr, end, argv[0]);
        arg_freetable(argtable, nb_opts);
        return 1;
    }


    microscopic_system_t sys;
    create_system(&sys, arg_datafile->filename[0] );

    if (arg_n->count ) {

        if( arg_n->ival[0] < 0){
            fprintf(stderr, "local_particule must be positive \n");
            arg_freetable(argtable, nb_opts);
            free_system(sys);
            return 1;
        }

        if( arg_n->ival[0] > sys.N_particules_total){
            fprintf(stderr, "local_particule must be inferior to the total number of particules in the datafile \n");
            arg_freetable(argtable, nb_opts);
            free_system(sys);
            return 1;
        }
        

        sys.N_particules_local = arg_n->ival[0];
        
    } 

    start_simulation( &sys, "../trajectories/traj.pdb", arg_i->ival[0], 
                     arg_t->dval[0], arg_p->ival[0], verbose->count, log->count);

    

    arg_freetable(argtable, nb_opts);

    free_system(sys);

    return 0;
}
