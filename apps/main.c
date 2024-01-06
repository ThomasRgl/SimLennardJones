#include "system.h"



int main() {

    microscopic_system_t sys;
    parse_system_data(&sys, "../data/particule.xyz");
    print_system(sys);
    free_system(sys);

    return 0;
}
