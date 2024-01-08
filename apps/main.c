#include "system.h"
#include "ljs.h"
#include <stdio.h>


int main() {
    int n_sym = 27;

    microscopic_system_t sys;
    parse_system_data(&sys, "../data/particule.xyz");
    // print_system(sys);
    printf("(%d) ljs potential : %f \n", n_sym, ljs_potential(&sys, n_sym) );

    ljs_forces( &sys, n_sym );
    double sum = ljs_sum_forces( &sys );
    printf("(%d) ljs sum forces : %g \n", n_sym, sum );
    
    free_system(sys);

    return 0;
}
