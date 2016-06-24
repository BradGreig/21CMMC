#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"

int main(int argc, char ** argv){

    char filename[500];
    FILE *F;
    
    int i,N_REDSHIFTS;
    
    float taue;
    
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
