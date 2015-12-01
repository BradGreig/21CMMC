
fftwf_complex *box, *k_box, *box_x, *updated, *save_updated, *deldel_T, *deltax_unfiltered, *deltax_unfiltered_original, *deltax_filtered, *deltax_unfiltered_t;
fftwf_complex *xe_unfiltered, *xe_filtered;
fftwf_complex *M_coll_unfiltered, *M_coll_filtered;

float *vx, *vy, *vz, *deltax_hires, *xH, *deltax, *Fcoll, *delta_T, *v, *smoothed_box, *smoothed_box_stored;
double *p_box, *k_ave;
unsigned long long *in_bin_ct;

float k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor;
int NUM_BINS;