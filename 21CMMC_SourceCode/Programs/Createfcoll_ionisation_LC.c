#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"
#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "gsl/gsl_sf_erf.h"

/* Throughout this and other 21cmMC drivers, the use of Deltac is not for checking against
 the z=0 collapsed fraction, but rather, it is used as a threshold for the validity of the
 collapse fraction expression. This expression is only valid up to Deltac
 */

float REDSHIFT;

void init_21cmMC_HII_arrays();
void destroy_21cmMC_HII_arrays();

int main(int argc, char ** argv){
    
    omp_set_num_threads(1);
    
    char dummy_string[500];
    
    char filename[500];
    FILE *F;
    fftwf_plan plan;
    
    // Other parameters used in the code
    int i,j,k, x,y,z, N_min_cell, LAST_FILTER_STEP, short_completely_ionised,skip_deallocate,first_step_R;
    unsigned long long ct;
    
    float MAX_RANGE_R, LOG_TVIR_MIN, LOG_TVIR_MAX, EFF_FACTOR_PL_INDEX, MIN_ALLOWED_ALPHA, MAX_ALLOWED_ALPHA;
    
    float growth_factor,MFEEDBACK, R, pixel_mass, cell_length_factor, ave_N_min_cell, M_MIN, nf;
    float f_coll_crit, erfc_denom, erfc_denom_cell, res_xH, Splined_Fcoll, sqrtarg, xHI_from_xrays, curr_dens, stored_R, massofscaleR;
    
    double global_xH, global_step_xH, ST_over_PS, mean_f_coll_st, f_coll, f_coll_temp;
    
    const gsl_rng_type * T;
    gsl_rng * r;
    
    skip_deallocate = 0;
    
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
    POWER_INDEX = (float)PARAM_COSMOLOGY_VALS[6];
    
    // Redshift for which Ts.c is evolved down to, i.e. z'
    REDSHIFT = atof(argv[3]);
    
    //R_MFP for evaluating f_coll
    R_BUBBLE_MAX = atof(argv[4]);
    
    //Tvir,min for evaluating f_coll
    ION_Tvir_MIN = atof(argv[5]);
    ION_Tvir_MIN = pow(10.,ION_Tvir_MIN);

    //Power-law index for evaluating f_coll
    EFF_FACTOR_PL_INDEX = atof(argv[6]);
    
    LAST_FILTER_STEP = atof(argv[7]);
    
    erfc_arg_min = -15.0;
    erfc_arg_max = 15.0;
    
    ERFC_NUM_POINTS = 10000;
    
    ERFC_VALS = calloc(ERFC_NUM_POINTS,sizeof(double));
    ERFC_VALS_DIFF = calloc(ERFC_NUM_POINTS,sizeof(double));
    
    ArgBinWidth = (erfc_arg_max - erfc_arg_min)/((double)ERFC_NUM_POINTS - 1.);
    InvArgBinWidth = 1./ArgBinWidth;
    
    for(i=0;i<ERFC_NUM_POINTS;i++) {
        
        erfc_arg_val = erfc_arg_min + ArgBinWidth*(double)i;
        
        ERFC_VALS[i] = splined_erfc(erfc_arg_val);
    }
    
    for(i=0;i<(ERFC_NUM_POINTS-1);i++) {
        ERFC_VALS_DIFF[i] = ERFC_VALS[i+1] - ERFC_VALS[i];
    }
    
    init_ps();
    
    /////////////////////////////////   BEGIN INITIALIZATION   //////////////////////////////////
    
    // perform a very rudimentary check to see if we are underresolved and not using the linear approx
    if ((BOX_LEN > DIM) && !EVOLVE_DENSITY_LINEARLY){
        printf("perturb_field.c: WARNING: Resolution is likely too low for accurate evolved density fields\n It Is recommended that you either increase the resolution (DIM/Box_LEN) or set the EVOLVE_DENSITY_LINEARLY flag to 1\n");
    }
    
    init_21cmMC_HII_arrays();
    
    // initialize power spectrum
    growth_factor = dicke(REDSHIFT);
    
    sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
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
    
    memcpy(deltax_unfiltered_original, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    
    i=0;
    
    pixel_mass = RtoM(L_FACTOR*BOX_LEN/(float)HII_DIM);
    cell_length_factor = L_FACTOR;
    
    // this parameter choice is sensitive to noise on the cell size, at least for the typical
    // cell sizes in RT simulations.  it probably doesn't matter for larger cell sizes.
    if (USE_HALO_FIELD && (FIND_BUBBLE_ALGORITHM==2) && ((BOX_LEN/(float)HII_DIM) < 1)){ // fairly arbitrary length based on 2 runs i did
        cell_length_factor = 1;
    }
    
    //set the minimum source mass
    if (ION_Tvir_MIN > 0){ // use the virial temperature for Mmin
        if (ION_Tvir_MIN < 9.99999e3) // neutral IGM
            M_MIN = TtoM(REDSHIFT, ION_Tvir_MIN, 1.22);
        else // ionized IGM
            M_MIN = TtoM(REDSHIFT, ION_Tvir_MIN, 0.6);
    }
    else if (ION_Tvir_MIN < 0){ // use the mass
        M_MIN = ION_M_MIN;
    }
    // check for WDM
    
    if (P_CUTOFF && ( M_MIN < M_J_WDM())){
        printf( "The default Jeans mass of %e Msun is smaller than the scale supressed by the effective pressure of WDM.\n", M_MIN);
        M_MIN = M_J_WDM();
        printf( "Setting a new effective Jeans mass from WDM pressure supression of %e Msun\n", M_MIN);
    }
    
    // lets check if we are going to bother with computing the inhmogeneous field at all...
    
    MFEEDBACK = M_MIN;
    
    if(EFF_FACTOR_PL_INDEX != 0.) {
        mean_f_coll_st = FgtrM_st_PL(REDSHIFT,M_MIN,MFEEDBACK,EFF_FACTOR_PL_INDEX);
    }
    else {
        mean_f_coll_st = FgtrM_st(REDSHIFT, M_MIN);
    }
    
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
    //  real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
        
    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
        deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
    }
    
    // loop through the filter radii (in Mpc)
    erfc_denom_cell=1; //dummy value
    //        R=fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN);
    
    initialiseSplinedSigmaM(M_MIN,1e16);
    
    R = R_BUBBLE_MAX;
    massofscaleR = RtoM(R);
    
    
    memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    HII_filter(deltax_filtered, HII_FILTER, R);
        
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
            
    // Check if this is the last filtering scale.  If so, we don't need deltax_unfiltered anymore.
    // We will re-read it to get the real-space field, which we will use to se the residual
    // neutral fraction
    ST_over_PS = 0;
    f_coll = 0;
    printf("Here?\n");
    if (LAST_FILTER_STEP){
        
        memcpy(deltax_unfiltered, deltax_unfiltered_original, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
                
        sqrtarg = 2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(massofscaleR), 2) );
        if (sqrtarg < 0){  // our filtering scale has become too small
            // Set to zero to avoid a seg. fault in the generation of the table.
            // However, the excursion set code itself will be able to deal with this anyway
            f_coll = 0.0;
        }
        else {
            erfc_denom_cell = sqrt(sqrtarg);
            erfc_denom_cell = 1./(growth_factor * sqrt(sqrtarg));
                
            if(EFF_FACTOR_PL_INDEX!=0.) {
                initialiseGL_Fcoll(NGLlow,NGLhigh,M_MIN,massofscaleR);
                initialiseFcoll_spline(REDSHIFT,M_MIN,massofscaleR,massofscaleR,MFEEDBACK,EFF_FACTOR_PL_INDEX);
            }
                
            // renormalize the collapse fraction so that the mean matches ST,
            // since we are using the evolved (non-linear) density field
            for (x=0; x<HII_DIM; x++){
                for (y=0; y<HII_DIM; y++){
                    for (z=0; z<HII_DIM; z++){
                        
                        curr_dens = *((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z));
                            
                        if(EFF_FACTOR_PL_INDEX!=0.) {
                            if(curr_dens < Deltac) {
                                FcollSpline(curr_dens,&(Splined_Fcoll));
                                    
                                f_coll += Splined_Fcoll;
                            }
                            else {
                                f_coll += 1.;
                            }
                        }
                        else {
                            
                            erfc_arg_val = (Deltac - curr_dens)*erfc_denom_cell;
                                
                            if( erfc_arg_val < erfc_arg_min || erfc_arg_val > erfc_arg_max ) {
                                f_coll_temp = splined_erfc((Deltac - curr_dens)*erfc_denom_cell);
                            }
                            else {
                                erfc_arg_val_index = (int)floor(( erfc_arg_val - erfc_arg_min )*InvArgBinWidth);
                                f_coll_temp = ERFC_VALS[erfc_arg_val_index] + (erfc_arg_val - (erfc_arg_min + ArgBinWidth*(double)erfc_arg_val_index))*ERFC_VALS_DIFF[erfc_arg_val_index]*InvArgBinWidth;
                            }
                            f_coll += f_coll_temp;
                        }
                    }
                }
            }
            f_coll /= (double) HII_TOT_NUM_PIXELS;
        } // end if last filter step conditional statement
    }
    else {         // not the last filter step, and we operating on the density field
    
        erfc_denom = 2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(massofscaleR), 2) );
        if (erfc_denom < 0){  // our filtering scale has become too small
            // Set to zero to avoid a seg. fault in the generation of the table.
            // However, the excursion set code itself will be able to deal with this anyway
            f_coll = 0.0;
        }
        else {
            erfc_denom = sqrt(erfc_denom);
            erfc_denom = 1./(erfc_denom * growth_factor);
                
            if(EFF_FACTOR_PL_INDEX!=0.) {
                initialiseGL_Fcoll(NGLlow,NGLhigh,M_MIN,massofscaleR);
                initialiseFcoll_spline(REDSHIFT,M_MIN,massofscaleR,massofscaleR,MFEEDBACK,EFF_FACTOR_PL_INDEX);
            }
                
            // renormalize the collapse fraction so that the mean matches ST,
            // since we are using the evolved (non-linear) density field
            for (x=HII_DIM; x--;){
                for (y=HII_DIM; y--;){
                    for (z=HII_DIM; z--;){
                    
                        curr_dens = *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));
                        if(EFF_FACTOR_PL_INDEX!=0.) {
                            if(curr_dens < Deltac) {
                                FcollSpline(curr_dens,&(Splined_Fcoll));
                                f_coll += Splined_Fcoll;
                            }
                            else {
                                f_coll += 1.0;
                            }
                        }
                        else {
                            erfc_arg_val = (Deltac - curr_dens)*erfc_denom;
                            
                            if( erfc_arg_val < erfc_arg_min || erfc_arg_val > erfc_arg_max ) {
                                f_coll_temp = splined_erfc((Deltac - curr_dens)*erfc_denom);
                            }
                            else {
                                erfc_arg_val_index = (int)floor(( erfc_arg_val - erfc_arg_min )*InvArgBinWidth);
                                f_coll_temp = ERFC_VALS[erfc_arg_val_index] + (erfc_arg_val - (erfc_arg_min + ArgBinWidth*(double)erfc_arg_val_index))*ERFC_VALS_DIFF[erfc_arg_val_index]*InvArgBinWidth;
                            }
                            f_coll += f_coll_temp;
                        }
                    }
                }
            }
            f_coll /= (double) HII_TOT_NUM_PIXELS;
        }
    }
    sprintf(filename, "Box_fcoll_z%1.6f_%1.6f_%1.6f_%1.6f_%d.txt",atof(argv[3]), R_BUBBLE_MAX, atof(argv[5]), EFF_FACTOR_PL_INDEX, LAST_FILTER_STEP);
    F=fopen(filename, "wt");
    fprintf(F, "%1.16lf\n",f_coll);
    fclose(F);
    
    destroy_21cmMC_HII_arrays(skip_deallocate);
    
    free(ERFC_VALS);
    free(ERFC_VALS_DIFF);
    
    free(redshifts);
    
    return 0;
}

/**** Arrays declared and used *****/

void init_21cmMC_HII_arrays() {
    
    Overdense_spline_GL_low = calloc(Nlow,sizeof(float));
    Fcoll_spline_GL_low = calloc(Nlow,sizeof(float));
    second_derivs_low_GL = calloc(Nlow,sizeof(float));
    Overdense_spline_GL_high = calloc(Nhigh,sizeof(float));
    Fcoll_spline_GL_high = calloc(Nhigh,sizeof(float));
    second_derivs_high_GL = calloc(Nhigh,sizeof(float));
    
    deltax_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_unfiltered_original = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    xe_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    xe_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    
    deltax = (float *) calloc(HII_TOT_FFT_NUM_PIXELS,sizeof(float));
    Fcoll = (float *) calloc(HII_TOT_NUM_PIXELS,sizeof(float));
    xH = (float *)calloc(HII_TOT_NUM_PIXELS,sizeof(float));
    delta_T = (float *) calloc(HII_TOT_NUM_PIXELS,sizeof(float));
    v = (float *) calloc(HII_TOT_FFT_NUM_PIXELS,sizeof(float));
    
    xi_low = calloc((NGLlow+1),sizeof(float));
    wi_low = calloc((NGLlow+1),sizeof(float));
    
    xi_high = calloc((NGLhigh+1),sizeof(float));
    wi_high = calloc((NGLhigh+1),sizeof(float));
    
//    k_factor = 1.25;
    k_factor = 1.35;
//    k_factor = 1.5;
//    k_factor = 2.0;
    k_first_bin_ceil = DELTA_K;
    k_max = DELTA_K*HII_DIM;
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

void destroy_21cmMC_HII_arrays(int skip_deallocate) {
    
    fftwf_free(deltax_unfiltered);
    fftwf_free(deltax_unfiltered_original);
    fftwf_free(deltax_filtered);
    fftwf_free(deldel_T);
    fftwf_free(xe_unfiltered);
    fftwf_free(xe_filtered);
    
    free(xH);
    free(deltax);
    free(Fcoll);
    free(delta_T);
    free(v);
    free(p_box);
    free(k_ave);
    free(in_bin_ct);
    
    free(Overdense_spline_GL_low);
    free(Fcoll_spline_GL_low);
    free(second_derivs_low_GL);
    free(Overdense_spline_GL_high);
    free(Fcoll_spline_GL_high);
    free(second_derivs_high_GL);

    free(xi_low);
    free(wi_low);
    
    free(xi_high);
    free(wi_high);
    
    if(skip_deallocate!=1) {
        free(Mass_Spline);
        free(Sigma_Spline);
        free(dSigmadm_Spline);
        free(second_derivs_sigma);
        free(second_derivs_dsigma);
    }
}