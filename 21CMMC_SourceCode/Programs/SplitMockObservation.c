#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "../Parameter_files/Variables.h"

// This file splits the 21cm cubic boxes used for the light-cone into smaller boxes (necessary for MCMC sampling) and computes the 21cm PS for each box
// This is useful if the light-cone box is to be sampled within 21CMMC

// This can be called as follows:
//
// ./SplitMockObservation filename BOX_LENGTH_LARGE HII_DIM_LARGE BOX_LENGTH_SMALL HII_DIM_SMALL

// 1) filename: this is the location of the 21cm brightness temperature box output by drive_21cmMC_streamlined.c (the flag PRINT_LIGHTCONE_21cmBoxes in Variables.h needs to be set to 1)

//    NOTE: Remember to set PRINT_LIGHTCONE_21cmBoxes back to 0 prior to commencing a full MCMC run!).

// 2) BOX_LENGTH_LARGE: Value of BOX_LEN used for creating the larger box (i.e. for the mock observation) from (INIT_PARAMS.H)
// 3) HII_DIM_LARGE: Value of HII_DIM used for creating the larger box (i.e. for the mock observation) from (INIT_PARAMS.H)
// 4) BOX_LENGTH_SMALL: Value of BOX_LEN that is to be used for sampling the smaller box within 21CMMC (i.e. the smaller box) from (INIT_PARAMS.H)
// 5) HII_DIM_SMALL: Value of HII_DIM used for creating the larger box (i.e. the smaller box) from (INIT_PARAMS.H)


// NOTE: This script only works provided that the ratio HII_DIM_LARGE/HII_DIM_SMALL is an integer multiple.


// This will output X 21cm PS files (X being the ratio of HII_DIM_LARGE/HII_DIM_SMALL) corresponding to the smaller boxes.


// NOTE: The output filenaming convention is incremented from 0 to X, therefore, these files need to be manually renamed for subsequent splitting of the 21cm brightness cubes (otherwise
//       the files will be overwritted). This is not ideal, but I couldn't think of an easier approach to do this (i.e. at least it keeps things quite general).


/* INDEXING MACROS */
#define LOS_RED_R_INDEX(x,y,z)((unsigned long long)((z)+DIM_MOCK_OBS_CUBIC*((y)+DIM_MOCK_OBS_CUBIC*(x)))) // for 3D real array with no padding

#define HII_LOS_RED_C_INDEX(x,y,z)((unsigned long long)((z)+(DIM_MOCK_OBS_MID+1llu)*((y)+DIM_MOCK_OBS_CUBIC*(x))))// for 3D complex array
#define HII_LOS_RED_R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(DIM_MOCK_OBS_MID+1llu)*((y)+DIM_MOCK_OBS_CUBIC*(x)))) // for 3D real array with the FFT padding
#define HII_LOS_RED_R_INDEX(x,y,z)((unsigned long long)((z)+DIM_MOCK_OBS*((y)+DIM_MOCK_OBS_CUBIC*(x)))) // for 3D real array with no padding

void init_21cmPS_arrays();

int main(int argc, char ** argv){
    
    FILE *F;
    
    char filename[1000];
    int i,j,k,ii,jj;
    
    // Provide the filename of the light-cone box to split
    strcpy(filename, argv[1]);
    
    printf("filename = %s\n",filename);
    
    // Provide the length (in Mpc) of the larger box (e.g. the mock observation)
    CUBIC_BOX_LENGTH = atof(argv[2]);
    
    // Provide HII_DIM of the larger box (e.g. the mock observation)
    DIM_MOCK_OBS_CUBIC = atof(argv[3]);
    DIM_MOCK_OBS_CUBIC_MID = DIM_MOCK_OBS_CUBIC/2;
    
    DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS = (unsigned long long)(DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC);
    
    // Provide the length (in Mpc) of the smaller box (e.g. the size of the box to be sampled within 21CMMC)
    RED_BOX_LENGTH = atof(argv[4]);
    
    // Provide HII_DIM of the smaller box (e.g. the size of the box to be sampled within 21CMMC)
    DIM_MOCK_OBS = atof(argv[5]);
    DIM_MOCK_OBS_MID = DIM_MOCK_OBS/2;
    
    DIM_MOCK_OBS_TOT_NUM_PIXELS = (unsigned long long)(DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS);
    DIM_MOCK_OBS_TOT_FFT_NUM_PIXELS = ((unsigned long long)(DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC*2llu*(DIM_MOCK_OBS_MID+1llu)));
    DIM_MOCK_OBS_KSPACE_NUM_PIXELS = ((unsigned long long)(DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC*(DIM_MOCK_OBS_MID+1llu)));

    printf("DIM_MOCK_OBS_CUBIC = %d DIM_MOCK_OBS = %d\n",DIM_MOCK_OBS_CUBIC,DIM_MOCK_OBS);
    printf("DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS = %d DIM_MOCK_OBS_TOT_NUM_PIXELS = %d\n",DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS,DIM_MOCK_OBS_TOT_NUM_PIXELS);
    printf("DIM_MOCK_OBS_TOT_FFT_NUM_PIXELS = %d DIM_MOCK_OBS_KSPACE_NUM_PIXELS = %d\n",DIM_MOCK_OBS_TOT_FFT_NUM_PIXELS,DIM_MOCK_OBS_KSPACE_NUM_PIXELS);
    
    printf("RATIO = %d %f\n",DIM_MOCK_OBS_CUBIC/DIM_MOCK_OBS,(float)DIM_MOCK_OBS_CUBIC/(float)DIM_MOCK_OBS);
    
    if((DIM_MOCK_OBS_CUBIC % DIM_MOCK_OBS)!=0) {
        printf("The L.O.S dimensions of the smaller (asymmetric) box is not an integer multiple of the larger cubic box containing the mock observation!\n");
        printf("This piece of code only works when the ratio is an integer multiple");
        return -1;
    }
    
    float *delta_T_LCBox,*delta_T_LCBox_Reduced;

    delta_T_LCBox = (float *) malloc(sizeof(float)*DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS);
    delta_T_LCBox_Reduced = (float *) malloc(sizeof(float)*DIM_MOCK_OBS_TOT_NUM_PIXELS);
    
    F = fopen(filename, "rb");
    for (i=0; i<DIM_MOCK_OBS_CUBIC; i++){
        for (j=0; j<DIM_MOCK_OBS_CUBIC; j++){
            for (k=0; k<DIM_MOCK_OBS_CUBIC; k++){
                mod_fread(delta_T_LCBox, sizeof(float)*DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS, 1, F);
            }
        }
    }
    fclose(F);
    
    double ave;
    int ratio, counter;
    int n_x, n_y, n_z;
    float k_x, k_y, k_z, k_mag;
    unsigned long long ct;

    ratio = DIM_MOCK_OBS_CUBIC/DIM_MOCK_OBS;
    
    counter = 0;

    fftwf_complex *deldel_T_asymmetric;
    fftwf_plan plan;
    
    deldel_T_asymmetric = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*DIM_MOCK_OBS_KSPACE_NUM_PIXELS);
    
    init_21cmPS_arrays();
    
    for(jj=0;jj<ratio;jj++) {
        
        for (ct=0; ct<NUM_BINS; ct++){
            p_box[ct] = k_ave[ct] = 0;
            in_bin_ct[ct] = 0;
        }
        
        ave = 0.0;
        for (i=0; i<DIM_MOCK_OBS_CUBIC; i++){
            for (j=0; j<DIM_MOCK_OBS_CUBIC; j++){
                for (k=0; k<DIM_MOCK_OBS; k++){
                    delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)] = delta_T_LCBox[LOS_RED_R_INDEX(i,j,k+counter*DIM_MOCK_OBS)];
                
                    ave += delta_T_LCBox[LOS_RED_R_INDEX(i,j,k+counter*DIM_MOCK_OBS)];
                }
            }
        }
        ave /= DIM_MOCK_OBS_TOT_NUM_PIXELS;
    
        for (i=0; i<DIM_MOCK_OBS_CUBIC; i++){
            for (j=0; j<DIM_MOCK_OBS_CUBIC; j++){
                for (k=0; k<DIM_MOCK_OBS; k++){
                    *((float *)deldel_T_asymmetric + HII_LOS_RED_R_FFT_INDEX(i,j,k)) = (delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)]/ave - 1)*(CUBIC_BOX_LENGTH*CUBIC_BOX_LENGTH*RED_BOX_LENGTH)/(DIM_MOCK_OBS_TOT_NUM_PIXELS+0.0);
                    if (DIMENSIONAL_T_POWER_SPEC){
                        *((float *)deldel_T_asymmetric + HII_LOS_RED_R_FFT_INDEX(i,j,k)) *= ave;
                    }
                    // Note: we include the V/N factor for the scaling after the fft
                }
            }
        }
    
        plan = fftwf_plan_dft_r2c_3d(DIM_MOCK_OBS_CUBIC, DIM_MOCK_OBS_CUBIC, DIM_MOCK_OBS, (float *)deldel_T_asymmetric, (fftwf_complex *)deldel_T_asymmetric, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
//      fftwf_cleanup();
    
        for (n_x=0; n_x<DIM_MOCK_OBS_CUBIC; n_x++){
            if (n_x>DIM_MOCK_OBS_CUBIC_MID)
                k_x =(n_x-(int)DIM_MOCK_OBS_CUBIC) * (TWOPI/CUBIC_BOX_LENGTH);  // wrap around for FFT convention
            else
                k_x = n_x * (TWOPI/CUBIC_BOX_LENGTH);
        
            for (n_y=0; n_y<DIM_MOCK_OBS_CUBIC; n_y++){
                
                // avoid the k(k_x = 0, k_y = 0, k_z) modes
                if(n_x != 0 && n_y != 0) {

                    if (n_y>DIM_MOCK_OBS_CUBIC_MID)
                        k_y =(n_y-(int)DIM_MOCK_OBS_CUBIC) * (TWOPI/CUBIC_BOX_LENGTH);
                    else
                        k_y = n_y * (TWOPI/CUBIC_BOX_LENGTH);
            
                    for (n_z=0; n_z<=DIM_MOCK_OBS_MID; n_z++){
                        k_z = n_z * (TWOPI/RED_BOX_LENGTH);
                
                        k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
                
                        // now go through the k bins and update
                        ct = 0;
                        k_floor = 0;
                        k_ceil = k_first_bin_ceil;
                        while (k_ceil < k_max){
                            // check if we fal in this bin
                            if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                                in_bin_ct[ct]++;
                                p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T_asymmetric[HII_LOS_RED_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*(CUBIC_BOX_LENGTH*CUBIC_BOX_LENGTH*RED_BOX_LENGTH));
                                // note the 1/VOLUME factor, which turns this into a power density in k-space
                        
                                k_ave[ct] += k_mag;
                                break;
                            }
                    
                            ct++;
                            k_floor=k_ceil;
                            k_ceil*=k_factor;
                        }
                    }
                }
            }
        }
        
        sprintf(filename, "delTps_estimate_carvedup_section_%d_%i_%.0fMpc_lighttravel.txt",jj, HII_DIM, BOX_LEN);
        F=fopen(filename, "wt");
        for (ct=1; ct<NUM_BINS; ct++){
            if (in_bin_ct[ct]>0) {
                fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
            }
        }
        fclose(F);
        
        for (i=0; i<DIM_MOCK_OBS_CUBIC; i++){
            for (j=0; j<DIM_MOCK_OBS_CUBIC; j++){
                for (k=0; k<(DIM_MOCK_OBS_CUBIC_MID+1); k++){
                    *((float *)deldel_T_asymmetric + HII_LOS_RED_R_FFT_INDEX(i,j,k)) = 0.0;
                }
            }
        }
        
        counter += 1;
    
        /*
        for (ct=1; ct<NUM_BINS; ct++){
            if (in_bin_ct[ct]>0)
                printf("%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
        }*/
        
        
        
        
    }
    
    return 0;
}

void init_21cmPS_arrays() {
    
    //    k_factor = 1.25;
    k_factor = 1.35;
    //    k_factor = 1.5;
    //    k_factor = 2.0;
    k_first_bin_ceil = (TWOPI/CUBIC_BOX_LENGTH);
    k_max = (TWOPI/CUBIC_BOX_LENGTH)*DIM_MOCK_OBS_CUBIC;
    // initialize arrays
    // ghetto counting (lookup how to do logs of arbitrary bases in c...)
    NUM_BINS = 0;
    k_floor = 0;
    k_ceil = k_first_bin_ceil;
    while (k_ceil < k_max){
        NUM_BINS++;
        k_floor=k_ceil;
        k_ceil*=k_factor;
    }
    
    p_box = calloc(NUM_BINS,sizeof(double));
    k_ave = calloc(NUM_BINS,sizeof(double));
    in_bin_ct = (unsigned long long *)calloc(NUM_BINS,sizeof(unsigned long long));
    
}


