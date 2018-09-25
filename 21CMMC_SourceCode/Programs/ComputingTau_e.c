#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"

int main(int argc, char ** argv){

    char filename[500];
    char dummy_string[500];
    FILE *F;
    
    int i,N_REDSHIFTS;
    
    float taue;
    
    INDIVIDUAL_ID = atof(argv[1]);
    INDIVIDUAL_ID_2 = atof(argv[2]);
    
    double *PARAM_COSMOLOGY_VALS = calloc(TOTAL_COSMOLOGY_FILEPARAMS,sizeof(double));
    
    /////////////////   Read in the cosmological parameter data     /////////////////
    
    sprintf(filename,"WalkerCosmology_%1.6lf_%1.6lf.txt",INDIVIDUAL_ID,INDIVIDUAL_ID_2);
    F = fopen(filename,"rt");
    
    for(i=0;i<TOTAL_COSMOLOGY_FILEPARAMS;i++) {
        fscanf(F,"%s\t%lf\n",&dummy_string,&PARAM_COSMOLOGY_VALS[i]);
    }
    fclose(F);
    
    // Assign these values. Hard-coded, so order is important
    RANDOM_SEED = (unsigned long long)PARAM_COSMOLOGY_VALS[0];
    SIGMA8 = (float)PARAM_COSMOLOGY_VALS[1];
    hlittle = (float)PARAM_COSMOLOGY_VALS[2];
    OMm = (float)PARAM_COSMOLOGY_VALS[3];
    OMl = (float)PARAM_COSMOLOGY_VALS[4];
    OMb = (float)PARAM_COSMOLOGY_VALS[5];
    POWER_INDEX = (float)PARAM_COSMOLOGY_VALS[6]; //power law on the spectral index, ns
    
    
    
    
    
    // Minus 3 for argv[0]=0, Random ID (argv[1]) and Zeta (argv[2])
    // Use Zeta just incase the Random ID ends up the same (shouldn't happen)
    N_REDSHIFTS = (argc-3)/2;
    
    float *Redshifts = calloc(N_REDSHIFTS,sizeof(float));
    float *xH = calloc(N_REDSHIFTS,sizeof(float));
    
    for(i=0;i<N_REDSHIFTS;i++) {
        Redshifts[i] = (float)(atof(argv[2*i+3]));
        xH[i] = (float)(atof(argv[2*(i+2)]));
    }
    
    taue = tau_e(0, Redshifts[N_REDSHIFTS-1], Redshifts, xH, N_REDSHIFTS);

//    printf("Tau = %lf\n",taue);
    
    sprintf(filename, "Tau_e_%s_%s.txt",argv[1],argv[2]);
    F=fopen(filename, "wt");
    fprintf(F, "%lf\n",taue);
    fclose(F);
    
    free(Redshifts);
    free(xH);
    
    return 0;
}
