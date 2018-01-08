#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"
#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "gsl/gsl_sf_erf.h"
#include "filter.c"

float REDSHIFT;

#define TOTAL_COSMOLOGY_FILEPARAMS (int)7

int main(int argc, char ** argv){
    
    omp_set_num_threads(1);
    
    char filename[500];
    char dummy_string[500];
    FILE *F;
    fftwf_plan plan;
    
    // Redshift for which Ts.c is evolved down to, i.e. z'
    REDSHIFT = atof(argv[1]);
    R_BUBBLE_MAX = atof(argv[2]);
    
    INDIVIDUAL_ID = atof(argv[3]);
    INDIVIDUAL_ID_2 = atof(argv[4]);
    
    int i,temp_int,temp_int2;
    
    double *PARAM_COSMOLOGY_VALS = calloc(TOTAL_COSMOLOGY_FILEPARAMS,sizeof(double));
    
    sprintf(filename,"WalkerCosmology_%1.6lf_%1.6lf.txt",INDIVIDUAL_ID,INDIVIDUAL_ID_2);
    F = fopen(filename,"rt");
    
    for(i=0;i<TOTAL_COSMOLOGY_FILEPARAMS;i++) {
        fscanf(F,"%s\t%lf\n",&dummy_string,&PARAM_COSMOLOGY_VALS[i]);
    }
    fclose(F);
    
    RANDOM_SEED = (unsigned long long)PARAM_COSMOLOGY_VALS[0];
    SIGMA8 = (float)PARAM_COSMOLOGY_VALS[1];
    hlittle = (float)PARAM_COSMOLOGY_VALS[2];
    OMm = (float)PARAM_COSMOLOGY_VALS[3];
    OMl = (float)PARAM_COSMOLOGY_VALS[4];
    OMb = (float)PARAM_COSMOLOGY_VALS[5];
    POWER_INDEX = (float)PARAM_COSMOLOGY_VALS[6];
    
    init_ps();
    
    float z_prime,prev_z_prime;
    int counter;
    
    int j,k,LAST_FILTER_STEP,first_step_R,N_RSTEPS;
    float cell_length_factor,stored_R;
    
    int n_x, n_y, n_z;
    float k_x, k_y, k_z, k_mag,kR;
    unsigned long long ct;
    
    dR = (BOX_LEN / (double) HII_DIM) * CMperMPC; // size of cell (in comoving cm)
    
    deltax_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    
    float **stored_smoothed_density;
    
    R=fmax(R_BUBBLE_MIN, (cell_length_factor*BOX_LEN/(float)HII_DIM));
    while (R < fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN)) {
        R*= DELTA_R_HII_FACTOR;
        if(R >= fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN)) {
            stored_R = R/DELTA_R_HII_FACTOR;
        }
        counter += 1;
    }
    
    N_RSTEPS = counter;
    
    stored_smoothed_density = (float **)calloc(N_RSTEPS,sizeof(float *));
    for(i=0;i<N_RSTEPS;i++) {
        stored_smoothed_density[i] = (float *)calloc(HII_TOT_FFT_NUM_PIXELS,sizeof(float));
    }
    
    // Determine the number of redshifts within the Ts.c calculation to set N_USER_REDSHIFT for the light-cone version of the computation.
        
    z_prime = REDSHIFT*1.0001; //higher for rounding
    while (z_prime < Z_HEAT_MAX)
        z_prime = ((1.+z_prime)*ZPRIME_STEP_FACTOR - 1.);
    prev_z_prime = Z_HEAT_MAX;
    z_prime = ((1.+z_prime)/ ZPRIME_STEP_FACTOR - 1.);
        
    while (z_prime > REDSHIFT){
        
        sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", z_prime, HII_DIM, BOX_LEN);
	printf("filename = %s\n",filename);
        F = fopen(filename, "rb");
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    if (fread((float *)deltax_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
                        printf("Read error occured while reading deltax box.\n");
                    }
                }
            }
        }
        fclose(F);
        
        for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
            deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
        }
        
	printf("Do FFT\n");
        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
//        fftwf_cleanup();
        
	printf("Complete FFT\n");

        cell_length_factor = L_FACTOR;
        
        counter = 0;
        
        
        R=fmax(R_BUBBLE_MIN, (cell_length_factor*BOX_LEN/(float)HII_DIM));
        while (R < fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN)) {
            R*= DELTA_R_HII_FACTOR;
            if(R >= fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN)) {
                stored_R = R/DELTA_R_HII_FACTOR;
            }
            counter += 1;
        }
        
        R=fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN);
        
        LAST_FILTER_STEP = 0;
        
        first_step_R = 1;
        
	printf("Start loop\n");
        while (!LAST_FILTER_STEP){//(R > (cell_length_factor*BOX_LEN/(HII_DIM+0.0))){
            
            if ( ((R/DELTA_R_HII_FACTOR - cell_length_factor*BOX_LEN/(float)HII_DIM) <= FRACT_FLOAT_ERR) || ((R/DELTA_R_HII_FACTOR - R_BUBBLE_MIN) <= FRACT_FLOAT_ERR) ) {
                LAST_FILTER_STEP = 1;
            }
	    printf("R = %1.6e counter = %d\n",R,counter-1);
            
            memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
            HII_filter(deltax_filtered, HII_FILTER, R);
	    printf("Before FFT\n");            
            plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
//            fftwf_cleanup();
            
	    printf("Fill array\n");

            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){
                        stored_smoothed_density[counter-1][HII_R_FFT_INDEX(i,j,k)] = *((float *)deltax_filtered + HII_R_FFT_INDEX(i,j,k));
                    }
                }
            }
            
            if(first_step_R) {
                R = stored_R;
                first_step_R = 0;
            }
            else {
                R /= DELTA_R_HII_FACTOR;
            }
            counter -= 1;
        }
  
	printf("Loop complete, write file\n");      
        sprintf(filename, "delta_x_filtered_RMAX%1.6f_z%1.6f",R_BUBBLE_MAX,z_prime);
        F = fopen(filename, "wb");
        for(i=0;i<N_RSTEPS;i++) {
            mod_fwrite(stored_smoothed_density[i], sizeof(float)*(HII_TOT_FFT_NUM_PIXELS), 1, F);
        }
        fclose(F);
        
        prev_z_prime = z_prime;
        z_prime = ((1.+prev_z_prime) / ZPRIME_STEP_FACTOR - 1.);
    }

    fftwf_free(deltax_unfiltered);
    fftwf_free(deltax_filtered);
    
    free(PARAM_COSMOLOGY_VALS);
    
    for(i=0;i<N_RSTEPS;i++) {
        free(stored_smoothed_density[i]);
    }
    free(stored_smoothed_density);

    
    return 0;
}
