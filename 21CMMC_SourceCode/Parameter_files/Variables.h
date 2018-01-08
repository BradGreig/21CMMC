// Globals for drive_21cmMC_streamlined.c

/* Maximum allowed value for the kinetic temperature. Useful to set to avoid some spurious behaviour
 when the code is run with redshift poor resolution and very high X-ray heating efficiency */
#define MAX_TK (float) 5e4

/* The minimum file length containing the relevant astrophysical parameters, plus related information requisite for computing the 21cm signal. This number should remain fixed,
 unless additional parameters are to be added and suitably varied. This number ** does not ** include the redshift sampling (i.e. if the co-eval box option is set. That comes
 through at the command line, where number of redshifts is provided.   */
#define TOTAL_AVAILABLE_PARAMS (int) (16)

/* A text file containing the cosmological parameters that can be varied within 21CMMC. This number is fixed, and includes the standard 6 parameter Lambda-CDM cosmological parameters,
 and also a random initial seed for generating the initial conditions (if this option is set). Note that OMEGA_L is assumed to be 1 - OMEGA_M, therefore it is unnecessary to vary OMEGA_L */
#define TOTAL_COSMOLOGY_FILEPARAMS (int)7

/* Allows the output of the global averaged quantities in the computation of the IGM spin temperature. Used for de-bugging purposes only */
#define OUTPUT_AVE  0

/* Whether or not to print the text files containing the neutral fraction or 21cm PS files */
#define PRINT_FILES (int)1

/* Whether or not to store the cubic 21cm boxes (brightness temperature) from the co-eval boxes */
#define PRINT_COEVAL_21cmBoxes (int)0

/* Whether or not to store the cubic 21cm boxes (brightness temperature) from the light-cone boxes */
#define PRINT_LIGHTCONE_21cmBoxes (int)0

/* Pre-determined number of parameters to generate the file strings for reading and using the f_coll ionisation table for reducing the computation of the light-cone boxes */
#define TOTAL_AVAILABLE_PARAMS_FCOLL_TABLE  (int)9

/* If one wants to compute the reduced light-cone only (i.e. not sample the full co-eval boxes at each output redshift from Ts.c) and the sub-cell redshift space distortions are added,
 then padding must be added to this reduced box to account for the fact that cells can move in/out of the reduced box. The default choice of 12 Mpc was based on my choice of 8 cells for
 the 1.5 Mpc voxels used in my work. This should not need to be changed */
#define LC_BOX_PADDING_IN_MPC (float)12


// Declaration of various required variables


fftwf_complex *box, *unfiltered_box, *deldel_T, *deldel_T_LC, *deltax_unfiltered, *deltax_unfiltered_original, *deltax_filtered, *xe_unfiltered, *xe_filtered;

fftwf_complex *HIRES_box, *HIRES_box_saved;

float *LOWRES_density, *LOWRES_vx, *LOWRES_vy, *LOWRES_vz, *LOWRES_vx_2LPT, *LOWRES_vy_2LPT, *LOWRES_vz_2LPT, *LOWRES_density_REDSHIFT, *LOWRES_velocity_REDSHIFT, *HIRES_density;

float *xH, *deltax, *Fcoll, *delta_T, *v, *vel_gradient, *zpp_growth, *inverse_diff, *Tk_box, *x_e_box, *Ts;
//float *delNL0_bw,*zpp_for_evolve_list,*R_values,*delNL0_Offset,*delNL0_LL,*delNL0_UL,*SingleVal_float,*delNL0_ibw,*log10delNL0_diff,*log10delNL0_diff_UL;
float *zpp_for_evolve_list,*R_values,*SingleVal_float;
float *delNL0_bw,*delNL0_Offset,*delNL0_LL,*delNL0_UL,*delNL0_ibw,*log10delNL0_diff,*log10delNL0_diff_UL;

float *box_z1, *box_z2, *box_interpolate, *box_interpolate_remainder;
double *redshifts_LC,*slice_redshifts;
int *start_index_LC, *end_index_LC,*full_index_LC;

float *x_pos_offset, *x_pos, *delta_T_RSD_LOS;

double *p_box, *k_ave, *fcoll_R_array, *Sigma_Tmin_grid, *ST_over_PS_arg_grid, *dstarlya_dt_prefactor, *zpp_edge, *sigma_atR, *sigma_Tmin, *ST_over_PS, *sum_lyn, *ERFC_VALS, *ERFC_VALS_DIFF;
double *aveTb, *aveNF, *redshifts;

unsigned long long *in_bin_ct;

double ***fcoll_R_grid, ***dfcoll_dz_grid;

double **density_gridpoints,**grid_dens, **freq_int_heat_tbl,**freq_int_ion_tbl,**freq_int_lya_tbl,**freq_int_heat_tbl_diff,**freq_int_ion_tbl_diff,**freq_int_lya_tbl_diff;

//float **fcoll_interp1, **fcoll_interp2, **dfcoll_interp1, **dfcoll_interp2, **Ts_z, **x_e_z, **delNL0_rev, **delNL0;
float *Ts_z, *x_e_z;
float **delNL0_rev, **delNL0;
double **fcoll_interp1, **fcoll_interp2, **dfcoll_interp1, **dfcoll_interp2;

double *Ionisation_fcoll_table, *Ionisation_fcoll_table_final;

short **dens_grid_int_vals, *SingleVal_int;

int NUM_BINS, CREATE_FFT_DATA_FROM_FILE,READ_FFT_DATA_FROM_FILE,SHORTEN_FCOLL,N_USER_REDSHIFT, WALKER_FILE_LENGTH, USE_LIGHTCONE,CALC_PS,USE_TS_FLUCT,INHOMO_RECO,STORE_DATA,ERFC_NUM_POINTS, erfc_arg_val_index;

float R_MFP_MIN, R_MFP_BINWIDTH, TVIR_BINWIDTH, PL_BINWIDTH, R_MFP_VAL_1, R_MFP_VAL_2, TVIR_VAL_1, TVIR_VAL_2, ZETA_PL_VAL_1, ZETA_PL_VAL_2;
int R_MFP_INT_1, R_MFP_INT_2, TVIR_INT_1, TVIR_INT_2, ZETA_PL_INT_1, ZETA_PL_INT_2;

float R_MFP_UB, TVIR_LB_FCOLL, TVIR_UB_FCOLL, ZETA_PL_LB, ZETA_PL_UB;
int R_MFP_STEPS, R_INTER_STEPS, TVIR_STEPS, PL_STEPS, SIZE_FIRST, SIZE_INTERMEDIATE, SIZE_FINAL, USE_FCOLL_IONISATION_TABLE, LC_BOX_PADDING, SUBCELL_RSD, GenerateNewICs, N_RSD_STEPS, LOS_direction;

int LOS_direction, slice_ct, total_slice_ct, num_boxes_interp,N_USER_REDSHIFT_LC,total_num_boxes,remainder_LC, Original_LOS_direction, Default_LOS_direction, Stored_LOS_direction_state_1, Stored_LOS_direction_state_2, N_RSTEPS_TOT;
double z1_LC, z2_LC, z_LC, dR, final_z;
double t_z1_LC,t_z2_LC,t_z_slice;
float start_z, end_z;

float k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor;
double INDIVIDUAL_ID, INDIVIDUAL_ID_2;

double EFF_FACTOR_PL_INDEX, HII_EFF_FACTOR, R_BUBBLE_MAX, ION_Tvir_MIN, L_X, NU_X_THRESH, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, NU_X_BAND_MAX, NU_X_MAX;
double BinWidth_pH,inv_BinWidth_pH,BinWidth_elec,inv_BinWidth_elec,BinWidth_10,inv_BinWidth_10, erfc_arg_min, erfc_arg_max, erfc_arg_val, ArgBinWidth, InvArgBinWidth;
double X_RAY_Tvir_LOWERBOUND, X_RAY_Tvir_UPPERBOUND,LOG10_X_RAY_Tvir_LOWERBOUND, LOG10_X_RAY_Tvir_UPPERBOUND, INCLUDE_ZETA_PL;

double F_STAR;
float t_STAR;

unsigned long long RANDOM_SEED;
float SIGMA8, hlittle, OMm, OMl, OMb, POWER_INDEX, INHOMO_RECO_R_BUBBLE_MAX;

float RED_BOX_LENGTH,CUBIC_BOX_LENGTH;

unsigned long long DIM_MOCK_OBS_CUBIC,DIM_MOCK_OBS,DIM_MOCK_OBS_MID,DIM_MOCK_OBS_CUBIC_MID,DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS,DIM_MOCK_OBS_TOT_NUM_PIXELS, DIM_MOCK_OBS_TOT_FFT_NUM_PIXELS,DIM_MOCK_OBS_KSPACE_NUM_PIXELS, DIM_MOCK_OBS_CUBIC_TOT_FFT_NUM_PIXELS, DIM_MOCK_OBS_CUBIC_KSPACE_NUM_PIXELS;

// Recombinations

float *z_re, *Gamma12;
fftwf_complex *N_rec_unfiltered, *N_rec_filtered;




