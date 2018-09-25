#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"
#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "gsl/gsl_sf_erf.h"

float REDSHIFT;

int main(int argc, char ** argv){
    
    omp_set_num_threads(1);
    
    char filename[500];
    char cmnd[1000];
    char dummy_string[500];
    FILE *F;
    
    INDIVIDUAL_ID = atof(argv[1]);
    INDIVIDUAL_ID_2 = atof(argv[2]);
    
    // Redshift for which Ts.c is evolved down to, i.e. z'
    REDSHIFT = atof(argv[3]);
    
    float z_prime,prev_z_prime;
    
    z_prime = REDSHIFT*1.0001; //higher for rounding
    while (z_prime < Z_HEAT_MAX)
        z_prime = ((1.+z_prime)*ZPRIME_STEP_FACTOR - 1.);
    prev_z_prime = Z_HEAT_MAX;
    z_prime = ((1.+z_prime)/ ZPRIME_STEP_FACTOR - 1.);
    
    while (z_prime > REDSHIFT){
    
        sprintf(cmnd, "./perturb_field %1.6f %1.6f %1.6f", z_prime,INDIVIDUAL_ID,INDIVIDUAL_ID_2);
        
        system(cmnd);
        
        prev_z_prime = z_prime;
        z_prime = ((1.+prev_z_prime) / ZPRIME_STEP_FACTOR - 1.);

    }

    return 0;
}
