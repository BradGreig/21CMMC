#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"
#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "gsl/gsl_sf_erf.h"

float REDSHIFT;

void init_21cmMC_FcollTable_arrays();

void destroy_21cmMC_FcollTable_arrays();

int main(int argc, char ** argv){

    char filename[500];
    char dummy_string[500];
    FILE *F;
    fftwf_plan plan;
    
    omp_set_num_threads((int)atof(argv[4]));
    
    // Redshift for which Ts.c is evolved down to, i.e. z'
    REDSHIFT = atof(argv[1]);
    
    X_RAY_Tvir_LOWERBOUND = atof(argv[2]);
    X_RAY_Tvir_UPPERBOUND = atof(argv[3]);
    
    LOG10_X_RAY_Tvir_LOWERBOUND = pow(10.,atof(argv[2]));
    LOG10_X_RAY_Tvir_UPPERBOUND = pow(10.,atof(argv[3]));
        
    init_21cmMC_FcollTable_arrays();
    
    init_ps();
    
    init_heat();
    
    clock_t start_section, end_section, start_subsection, end_subsection;
    
    unsigned long long ct;
    
    int R_ct,i,j,i_Tvir,ii,k, zpp_gridpoint1_int, zpp_gridpoint2_int,zpp_evolve_gridpoint1_int, zpp_evolve_gridpoint2_int;
    float growth_factor_z, inverse_growth_factor_z, R, R_factor, zp, mu_for_Ts, filling_factor_of_HI_zp, dzp, prev_zp, zpp, prev_zpp, prev_R, Tk_BC, xe_BC;
    
    float determine_zpp_max, determine_zpp_min, min_density, max_density, zpp_grid;
    float zpp_gridpoint1, zpp_gridpoint2,zpp_evolve_gridpoint1, zpp_evolve_gridpoint2;
    float grad1, grad2, grad3, grad4, OffsetValue, DensityValueLow, delNL0_bw_val, delNL0_ibw_val, log10delNL0_diff_val;

    double dens_grad, grid_dens;
    
    float zpp_for_evolve,dzpp_for_evolve,zpp_bin_width;
    double curr_delNL0, inverse_val,prefactor_1,prefactor_2,dfcoll_dz_val, density_eval1, density_eval2, grid_sigmaTmin, fcoll_R;
    
    int counter = 0;
    
    float M_MIN_WDM =  M_J_WDM();
    
    growth_factor_z = dicke(REDSHIFT);
    inverse_growth_factor_z = 1./growth_factor_z;
    
    /////////////// Create the z=0 non-linear density fields smoothed on scale R to be used in computing fcoll //////////////
    R = L_FACTOR*BOX_LEN/(float)HII_DIM;
    R_factor = pow(R_XLy_MAX/R, 1/(float)NUM_FILTER_STEPS_FOR_Ts);
    //      R_factor = pow(E, log(HII_DIM)/(float)NUM_FILTER_STEPS_FOR_Ts);
    
    ///////////////////  Read in density box at z-prime  ///////////////
    
    // allocate memory for the nonlinear density field and open file
    sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc",REDSHIFT, HII_DIM, BOX_LEN);
    F = fopen(filename, "rb");
    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                fread((float *)unfiltered_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F);
            }
        }
    }
    fclose(F);
    
    ////////////////// Transform unfiltered box to k-space to prepare for filtering /////////////////
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)unfiltered_box, (fftwf_complex *)unfiltered_box, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    
    // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
        unfiltered_box[ct] /= (float)HII_TOT_NUM_PIXELS;
    }
    
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){

        R_values[R_ct] = R;
        sigma_atR[R_ct] = sigma_z0(RtoM(R));
        
        // copy over unfiltered box
        memcpy(box, unfiltered_box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
        
        if (R_ct > 0){ // don't filter on cell size
            HII_filter(box, HEAT_FILTER, R);
        }
        
        // now fft back to real space
        plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)box, (float *)box, FFTW_ESTIMATE);
        fftwf_execute(plan);
        
        min_density = 0.0;
        max_density = 0.0;
        
        // copy over the values
        for (i=HII_DIM; i--;){
            for (j=HII_DIM; j--;){
                for (k=HII_DIM; k--;){
                    curr_delNL0 = *((float *) box + HII_R_FFT_INDEX(i,j,k));
                    
                    if (curr_delNL0 < -1){ // correct for alliasing in the filtering step
                        curr_delNL0 = -1+FRACT_FLOAT_ERR;
                    }
                    // and linearly extrapolate to z=0
                    curr_delNL0 *= inverse_growth_factor_z;
                    
                    delNL0[R_ct][HII_R_INDEX(i,j,k)] = curr_delNL0;
                    
                    if(curr_delNL0 < min_density) {
                        min_density = curr_delNL0;
                    }
                    if(curr_delNL0 > max_density) {
                        max_density = curr_delNL0;
                    }
                    
                }
            }
        }
        
        if(min_density < 0.0) {
            delNL0_LL[R_ct] = min_density*1.001;
            delNL0_Offset[R_ct] = 1.e-6 - (delNL0_LL[R_ct]);
        }
        else {
            delNL0_LL[R_ct] = min_density*0.999;
            delNL0_Offset[R_ct] = 1.e-6 + (delNL0_LL[R_ct]);
        }
        if(max_density < 0.0) {
            delNL0_UL[R_ct] = max_density*0.999;
        }
        else {
            delNL0_UL[R_ct] = max_density*1.001;
        }
        
        R *= R_factor;
        
    } //end for loop through the filter scales R
    
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    
    for(ii=0;ii<NUM_FILTER_STEPS_FOR_Ts;ii++) {
        delNL0_bw[ii] = ( log10(delNL0_UL[ii] + delNL0_Offset[ii]) - log10(delNL0_LL[ii] + delNL0_Offset[ii]) )/((float)dens_Ninterp - 1.);
    }
    
    log10delNL0_diff_UL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    
    for (R_ct=NUM_FILTER_STEPS_FOR_Ts; R_ct--;){
        log10delNL0_diff_UL[R_ct] = log10( delNL0_UL[R_ct] + delNL0_Offset[R_ct] );
    }
    
    // Determine the grid point locations for solving the interpolation tables
    for (R_ct=NUM_FILTER_STEPS_FOR_Ts; R_ct--;){
        delNL0_ibw_val = 1./delNL0_bw[R_ct];
        log10delNL0_diff_val = log10( delNL0_LL[R_ct] + delNL0_Offset[R_ct] );
        log10delNL0_diff[R_ct] = log10delNL0_diff_val;
        OffsetValue = delNL0_Offset[R_ct];
        
        for (box_ct=HII_TOT_NUM_PIXELS; box_ct--;){
            dens_grid_int_vals[R_ct][box_ct] = (short)floor( ( log10(delNL0[R_ct][box_ct] + OffsetValue ) - log10delNL0_diff_val )*delNL0_ibw_val);
        }
    }
    
    // Evaluating the interpolated density field points (for using the interpolation tables for fcoll and dfcoll_dz
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
        OffsetValue = delNL0_Offset[R_ct];
        DensityValueLow = delNL0_LL[R_ct];
        delNL0_bw_val = delNL0_bw[R_ct];
        
        for(i=0;i<dens_Ninterp;i++) {
            density_gridpoints[R_ct][i] = pow(10.,( log10( DensityValueLow + OffsetValue) + delNL0_bw_val*((float)i) )) - OffsetValue;
        }
    }
    
    // main trapezoidal integral over z' (see eq. ? in Mesinger et al. 2009)
    zp = REDSHIFT*1.0001; //higher for rounding
    while (zp < Z_HEAT_MAX)
        zp = ((1+zp)*ZPRIME_STEP_FACTOR - 1);
    prev_zp = Z_HEAT_MAX;
    zp = ((1+zp)/ ZPRIME_STEP_FACTOR - 1);
    dzp = zp - prev_zp;
    
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
        if (R_ct==0){
            prev_zpp = zp;
            prev_R = 0;
        }
        else{
            prev_zpp = zpp_edge[R_ct-1];
            prev_R = R_values[R_ct-1];
        }
        zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
        zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
    }
    
    determine_zpp_min = REDSHIFT*0.999;
    determine_zpp_max = zpp*1.001;
    
    zpp_bin_width = (determine_zpp_max - determine_zpp_min)/((float)zpp_interp_points-1.0);
    
    while (zp > REDSHIFT){
        
        counter += 1.;
        
        // let's initialize an array of redshifts (z'') corresponding to the
        // far edge of the dz'' filtering shells
        // and the corresponding minimum halo scale, sigma_Tmin,
        // as well as an array of the frequency integrals
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
            if (R_ct==0){
                prev_zpp = zp;
            }
            else{
                prev_zpp = zpp_edge[R_ct-1];
            }
        } // end loop over R_ct filter steps
        
        prev_zp = zp;
        zp = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
        dzp = zp - prev_zp;
        
    } // end main integral loop over z'
    
    int TOT_NUM_ZPRIME = counter;
    
    double ***Fcoll_R_Table = (double ***)calloc(TOT_NUM_ZPRIME,sizeof(double **));
    for(i=0;i<TOT_NUM_ZPRIME;i++) {
        Fcoll_R_Table[i] = (double **)calloc(X_RAY_Tvir_POINTS,sizeof(double *));
        for(j=0;j<X_RAY_Tvir_POINTS;j++) {
            Fcoll_R_Table[i][j] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
        }
    }
    
    for(i_Tvir=0;i_Tvir<X_RAY_Tvir_POINTS;i_Tvir++) {
        
        X_RAY_Tvir_MIN = X_RAY_Tvir_LOWERBOUND + (X_RAY_Tvir_UPPERBOUND - X_RAY_Tvir_LOWERBOUND)*(double)i_Tvir/( (double)X_RAY_Tvir_POINTS - 1. );
        X_RAY_Tvir_MIN = pow(10.,X_RAY_Tvir_MIN);
    
        //NOTE: This assumes ION_Tvir_MIN =X_RAY_Tvir_MIN;
        ION_Tvir_MIN = X_RAY_Tvir_MIN;
        
        if (X_RAY_Tvir_MIN < 9.99999e3) // neutral IGM
            mu_for_Ts = 1.22;
        else // ionized IGM
            mu_for_Ts = 0.6;
        
        // main trapezoidal integral over z' (see eq. ? in Mesinger et al. 2009)
        zp = REDSHIFT*1.0001; //higher for rounding
        while (zp < Z_HEAT_MAX)
            zp = ((1+zp)*ZPRIME_STEP_FACTOR - 1);
        prev_zp = Z_HEAT_MAX;
        zp = ((1+zp)/ ZPRIME_STEP_FACTOR - 1);
        dzp = zp - prev_zp;
        
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
            if (R_ct==0){
                prev_zpp = zp;
                prev_R = 0;
            }
            else{
                prev_zpp = zpp_edge[R_ct-1];
                prev_R = R_values[R_ct-1];
            }
            zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
            zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
        }
        
        ////////////////////////////    Create and fill interpolation tables to be used by Ts.c   /////////////////////////////
        
        init_FcollTable(determine_zpp_min,determine_zpp_max);
        
        for(i=0;i<zpp_interp_points;i++) {
            
            zpp_grid = determine_zpp_min + (determine_zpp_max - determine_zpp_min)*(float)i/((float)zpp_interp_points-1.0);
            
            grid_sigmaTmin = sigma_z0(FMAX(TtoM(zpp_grid, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM));
            
            Sigma_Tmin_grid[i] = grid_sigmaTmin;
            
            for(ii=0;ii<NUM_FILTER_STEPS_FOR_Ts;ii++) {
                for(j=0;j<dens_Ninterp;j++) {
                    
                    grid_dens = log10delNL0_diff[ii] + ( log10delNL0_diff_UL[ii] - log10delNL0_diff[ii] )*(float)j/((float)dens_Ninterp - 1.);
                    grid_dens = pow(10,grid_dens) - delNL0_Offset[ii];
                    
                    fcoll_R_grid[ii][i][j] = sigmaparam_FgtrM_bias(zpp_grid, grid_sigmaTmin, grid_dens, sigma_atR[ii]);
                }
            }
        }
        
        counter = 0;
        
        while (zp > REDSHIFT){
            
            // let's initialize an array of redshifts (z'') corresponding to the
            // far edge of the dz'' filtering shells
            // and the corresponding minimum halo scale, sigma_Tmin,
            // as well as an array of the frequency integrals
            
            for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
                if (R_ct==0){
                    prev_zpp = zp;
                    prev_R = 0;
                }
                else{
                    prev_zpp = zpp_edge[R_ct-1];
                    prev_R = R_values[R_ct-1];
                }
                
                zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
                zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
                
                zpp_gridpoint1_int = (int)floor((zpp - determine_zpp_min)/zpp_bin_width);
                zpp_gridpoint2_int = zpp_gridpoint1_int + 1;
                
                zpp_gridpoint1 = determine_zpp_min + zpp_bin_width*(float)zpp_gridpoint1_int;
                zpp_gridpoint2 = determine_zpp_min + zpp_bin_width*(float)zpp_gridpoint2_int;
                
                grad1 = ( zpp_gridpoint2 - zpp )/( zpp_gridpoint2 - zpp_gridpoint1 );
                grad2 = ( zpp - zpp_gridpoint1 )/( zpp_gridpoint2 - zpp_gridpoint1 );
                
                sigma_Tmin[R_ct] = Sigma_Tmin_grid[zpp_gridpoint1_int] + grad2*( Sigma_Tmin_grid[zpp_gridpoint2_int] - Sigma_Tmin_grid[zpp_gridpoint1_int] );
                
                // let's now normalize the total collapse fraction so that the mean is the
                // Sheth-Torman collapse fraction
                
                for(i=0;i<(dens_Ninterp-1);i++) {
                    dens_grad = 1./( density_gridpoints[R_ct][i+1] - density_gridpoints[R_ct][i] );

                    fcoll_interp1[R_ct][i] = ( ( fcoll_R_grid[R_ct][zpp_gridpoint1_int][i] )*grad1 + ( fcoll_R_grid[R_ct][zpp_gridpoint2_int][i] )*grad2 )*dens_grad;
                    fcoll_interp2[R_ct][i] = ( ( fcoll_R_grid[R_ct][zpp_gridpoint1_int][i+1] )*grad1 + ( fcoll_R_grid[R_ct][zpp_gridpoint2_int][i+1] )*grad2 )*dens_grad;
                }
            } // end loop over R_ct filter steps

            for (R_ct=0;R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
                
                fcoll_R = 0.0;
                
#pragma omp parallel shared (box_ct,fcoll_R_array,fcoll_interp1,fcoll_interp2,dens_grid_int_vals,density_gridpoints,delNL0_rev)
                {
#pragma omp for schedule(static,1) reduction(+:fcoll_R)
                    for (box_ct=0;box_ct<HII_TOT_NUM_PIXELS; box_ct++){
                        fcoll_R += ( fcoll_interp1[R_ct][dens_grid_int_vals[R_ct][box_ct]]*( density_gridpoints[R_ct][dens_grid_int_vals[R_ct][box_ct] + 1] - delNL0[R_ct][box_ct] ) + fcoll_interp2[R_ct][dens_grid_int_vals[R_ct][box_ct]]*( delNL0[R_ct][box_ct] - density_gridpoints[R_ct][dens_grid_int_vals[R_ct][box_ct]] ) );
                    }
                }
                Fcoll_R_Table[counter][i_Tvir][R_ct] = fcoll_R/(double)HII_TOT_NUM_PIXELS;
            }
            
            prev_zp = zp;
            zp = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
            dzp = zp - prev_zp;
            
            counter += 1;
        } // end main integral loop over z'
    }
    
    sprintf(filename, "FcollTvirTable_Numzp_ZPRIME_FACTOR%0.2f_logTvirmin%0.6f_logTvirmax%0.6f_XRAY_POINTS%d_z_end%06.6f_%0.2fMpc_%d.dat",ZPRIME_STEP_FACTOR,X_RAY_Tvir_LOWERBOUND,X_RAY_Tvir_UPPERBOUND,X_RAY_Tvir_POINTS,REDSHIFT,BOX_LEN,HII_DIM);
    F = fopen(filename, "wb");
    fwrite(&counter, sizeof(int),1,F);
    fclose(F);

    sprintf(filename, "FcollTvirTable_ZPRIME_FACTOR%0.2f_logTvirmin%0.6f_logTvirmax%0.6f_XRAY_POINTS%d_z_end%06.6f_%0.2fMpc_%d.dat",ZPRIME_STEP_FACTOR,X_RAY_Tvir_LOWERBOUND,X_RAY_Tvir_UPPERBOUND,X_RAY_Tvir_POINTS,REDSHIFT,BOX_LEN,HII_DIM);
    F = fopen(filename, "wb");
    for(i=0;i<counter;i++) {
        for(i_Tvir=0;i_Tvir<X_RAY_Tvir_POINTS;i_Tvir++) {
            fwrite(Fcoll_R_Table[i][i_Tvir], sizeof(double),NUM_FILTER_STEPS_FOR_Ts,F);
        }
    }
    fclose(F);
    
    for(i=0;i<TOT_NUM_ZPRIME;i++) {
        for(j=0;j<X_RAY_Tvir_POINTS;j++) {
            free(Fcoll_R_Table[i][j]);
        }
        free(Fcoll_R_Table[i]);
    }
    free(Fcoll_R_Table);
    
    destroy_21cmMC_FcollTable_arrays();
    
    destruct_heat();
    
    return 0;
}


/**** Arrays declared and used *****/

void init_21cmMC_FcollTable_arrays() {

    int i,j;
    
    box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    unfiltered_box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    
    fcoll_R_grid = (double ***)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double **));
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        fcoll_R_grid[i] = (double **)calloc(zpp_interp_points,sizeof(double *));
        for(j=0;j<zpp_interp_points;j++) {
            fcoll_R_grid[i][j] = (double *)calloc(dens_Ninterp,sizeof(double));
        }
    }
    
    Sigma_Tmin_grid = (double *)calloc(zpp_interp_points,sizeof(double));
    
    density_gridpoints = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double *));
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        density_gridpoints[i] = (double *)calloc(dens_Ninterp,sizeof(double));
    }
    
    dens_grid_int_vals = (short **)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(short *));
    delNL0 = (float **)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float *));
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        dens_grid_int_vals[i] = (short *)calloc((float)HII_TOT_NUM_PIXELS,sizeof(short));
        delNL0[i] = (float *)calloc((float)HII_TOT_NUM_PIXELS,sizeof(float));
    }
    
    fcoll_interp1 = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double *));
    fcoll_interp2 = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double *));
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        fcoll_interp1[i] = (double *)calloc(dens_Ninterp,sizeof(double));
        fcoll_interp2[i] = (double *)calloc(dens_Ninterp,sizeof(double));
    }

    zpp_edge = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    sigma_atR = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    sigma_Tmin = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    
    delNL0_bw = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    R_values = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    delNL0_Offset = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    delNL0_LL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    delNL0_UL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    SingleVal_float = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    delNL0_ibw = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    log10delNL0_diff = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
}

void destroy_21cmMC_FcollTable_arrays() {

    int i,j;
    
    free(box);
    free(unfiltered_box);
    
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        for(j=0;j<zpp_interp_points;j++) {
            free(fcoll_R_grid[i][j]);
        }
        free(fcoll_R_grid[i]);
    }
    free(fcoll_R_grid);
    
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        free(density_gridpoints[i]);
    }
    free(density_gridpoints);
    
    free(ST_over_PS_arg_grid);
    free(Sigma_Tmin_grid);
    
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        free(fcoll_interp1[i]);
        free(fcoll_interp2[i]);
        free(dens_grid_int_vals[i]);
        free(delNL0[i]);
    }
    free(fcoll_interp1);
    free(fcoll_interp2);
    free(dens_grid_int_vals);
    free(delNL0);
    
    free(delNL0_bw);
    free(R_values);
    free(delNL0_Offset);
    free(delNL0_LL);
    free(delNL0_UL);
    free(SingleVal_int);
    free(delNL0_ibw);
    free(log10delNL0_diff);

    free(zpp_edge);
    free(sigma_atR);
    free(sigma_Tmin);
    
    free_FcollTable();
}
