#include "system.h"


double distance( microscopic_system_t * sys, const size_t i, const size_t j){
    double dist_x = ( sys->px[i] - sys->px[j] );
    double dist_y = ( sys->py[i] - sys->py[j] );
    double dist_z = ( sys->pz[i] - sys->pz[j] );

    return sqrt( dist_x * dist_x + dist_y * dist_y + dist_z * dist_z );
}

double squared_distance( microscopic_system_t * sys, const size_t i, const size_t j){
    double dist_x = ( sys->px[i] - sys->px[j] );
    double dist_y = ( sys->py[i] - sys->py[j] );
    double dist_z = ( sys->pz[i] - sys->pz[j] );

    return dist_x * dist_x + dist_y * dist_y + dist_z * dist_z ;
}



int print_system( microscopic_system_t sys ){

    printf("x\ty\tz\n");
    for (size_t i = 0; i < sys.N_particules_total; i++) {
        printf("%lu\t%lf\t%lf\t%lf\n",i, sys.px[i], sys.py[i], sys.pz[i]);
    }
    return 0;
}

int allocate_system( microscopic_system_t * sys, size_t indice ){
 
    sys->N_particules_total = indice;
    sys->N_particules_local = indice;

    sys->px = aligned_alloc(32, indice * sizeof(double) );
    sys->py = aligned_alloc(32, indice * sizeof(double) );
    sys->pz = aligned_alloc(32, indice * sizeof(double) );

    sys->fx = aligned_alloc(32, indice * sizeof(double) );
    sys->fy = aligned_alloc(32, indice * sizeof(double) );
    sys->fz = aligned_alloc(32, indice * sizeof(double) );

       
    return 0;
}

int free_system( microscopic_system_t sys ){
    free(sys.px);
    free(sys.py);
    free(sys.pz);

    free(sys.fx);
    free(sys.fy);
    free(sys.fz);


    return 0;
}


int parse_system_data( microscopic_system_t * sys, char * filename ){
    FILE * fp;
    fp = fopen(filename, "r");

    if (fp == NULL) {
        perror("error when opening particule system file");
        return 1;
    }

    size_t indice = 0;
    double * x = malloc( MAX_LINES * sizeof(double) );
    double * y = malloc( MAX_LINES * sizeof(double) );
    double * z = malloc( MAX_LINES * sizeof(double) );

    // Ignorer la première ligne
    fscanf(fp, "%*d %*d");

    // Lire les données et stocker dans les tableaux
    while (fscanf(fp, "%*d %lf %lf %lf", &x[indice], &y[indice], &z[indice]) == 3) {
        indice++;
    }

    // Fermer le fp
    fclose(fp);

    allocate_system(sys, indice);
    
    for (size_t i = 0; i < indice; i++) {
        sys->px[i] = x[i];
        sys->py[i] = y[i];
        sys->pz[i] = z[i];
    }

    

    free(x);
    free(y);
    free(z);
    return 0;
}
