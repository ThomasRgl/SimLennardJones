#include "system.h"
#include "ljs.h"
#include <stdio.h>



int main() {

    microscopic_system_t sys;
    parse_system_data(&sys, "../data/particule.xyz");
    print_system(sys);
    printf("ljs potential : %f \n", ljs_potential(&sys) );

    ljs_forces( &sys );
    double sum = ljs_sum_forces( &sys );
    printf("ljs sum forces : %.30f \n", sum );
    
    free_system(sys);

    return 0;
}
