#include "io.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

//  1 -  6        Record name   "ATOM  "
//  7 - 11        Integer       serial       Atom  serial number.
// 13 - 16        Atom          name         Atom name.
// 17             Character     altLoc       Alternate location indicator.
// 18 - 20        Residue name  resName      Residue name.
// 22             Character     chainID      Chain identifier.
// 23 - 26        Integer       resSeq       Residue sequence number.
// 27             AChar         iCode        Code for insertion of residues.
// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
//         
int save_trajectory( FILE * fp, microscopic_system_t * sys, int ite){
    
    //Crystallographic

    const double alpha = 90; 
    const double beta = 90; 
    const double gamma = 90; 
    const char * sgroup = " P             ";
    const int zvalue = 1;

    fprintf(fp, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%d\n", sys->dimx, sys->dimy, sys->dimz, alpha, beta, gamma, sgroup, zvalue);

    // model
    fprintf(fp, "MODEL     %d\n", ite);

    for ( int i = 0; i < sys->N_particules_local; i++ ) {
        fprintf(fp, "ATOM  "); // 1 - 6
        fprintf(fp, "%5d", i+1); // 7 - 11  Atom  serial number
        fprintf(fp, " ");// 12
        fprintf(fp, "  C "); // 13 - 16 atome name
        fprintf(fp, " "); // 17 alternate location indicatior
        fprintf(fp, "   "); // 18-20 Residue name
        fprintf(fp, " "); // 21 chain id
        fprintf(fp, " "); // 22
        fprintf(fp, "    "); // 23 - 26 resSeq
        fprintf(fp, "0"); // 27 iCode 
        fprintf(fp, "   "); // 28-30 
        fprintf(fp, "%8.3f", sys->x[i] ); // 31-38 x
        fprintf(fp, "%8.3f", sys->y[i] ); // 39-46y 
        fprintf(fp, "%8.3f", sys->z[i] ); // 47-54 z
        fprintf(fp,"      "); // 55-60 occupancy
        fprintf(fp,"      "); // 61-66 temp factor
        fprintf(fp,"          "); // 67-76
        fprintf(fp,"MR"); // 77-78
        fprintf(fp,"ES"); // 79-80
        fprintf(fp, "\n");
    }

    // end
    fprintf(fp, "TER\n");
    fprintf(fp, "ENDMDL\n");

    return 0;
}

size_t parse_data( const char * filename, double ** x, double ** y, double ** z ){

    FILE * fp;
    fp = fopen(filename, "r");

    if (fp == NULL) {
        perror("error when opening particule system file");
        return 1;
    }

    size_t indice = 0;
    *x = malloc( MAX_LINES * sizeof(double) );
    *y = malloc( MAX_LINES * sizeof(double) );
    *z = malloc( MAX_LINES * sizeof(double) );

    // Ignorer la première ligne
    fscanf(fp, "%*d %*d");

    // Lire les données et stocker dans les tableaux
    while (fscanf(fp, "%*d %lf %lf %lf", &(*x)[indice], &(*y)[indice], &(*z)[indice]) == 3) {
        indice++;
    }

    // Fermer le fp
    fclose(fp);

    return indice;

}
