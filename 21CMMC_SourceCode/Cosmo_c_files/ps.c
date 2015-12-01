#ifndef _PS_
#define _PS_

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "../Parameter_files/COSMOLOGY.H"
#include "../Parameter_files/INIT_PARAMS.H"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "cosmo_progs.c"
#include "misc.c"

/* New in v1.1 */
#define ERFC_NPTS (int) 75
#define ERFC_PARAM_DELTA (float) 0.1
static double log_erfc_table[ERFC_NPTS], erfc_params[ERFC_NPTS];
static gsl_interp_accel *erfc_acc;
static gsl_spline *erfc_spline;

double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;


/*****     FUNCTION PROTOTYPES     *****/
double init_ps(); /* initialize global variables, MUST CALL THIS FIRST!!! returns R_CUTOFF */
void free_ps(); /* deallocates the gsl structures from init_ps */
double splined_erfc(double); /* returns erfc for x>=0, using cubic spline in logy-x space */
double deltolindel(float del, float z); /* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double lindeltodel(float lindel, float z); /* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double power_in_k(double k); /* Returns the value of the linear power spectrum density (i.e. <|delta_k|^2>/V) at a given k mode at z=0 */
double RtoM(double); /* R in Mpc, M in Msun */
double MtoR(double); /* R in Mpc, M in Msun */
double M_J_WDM(); /* returns the "effective Jeans mass" corresponding to the gas analog of WDM ; eq. 10 in BHO 2001 */
double sheth_delc(double del, double sig);
double dNdM_st(double z, double M);
double dNdM(double z, double M);
double dnbiasdM(double M, float z, double M_o, float del_o); /* dnbiasdM */
double FgtrM(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z
double FgtrM_st(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, with Sheth-Tormen correction
double FgtrM_bias(double z, double M, double del_bias, double sig_bias);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias);/* Uses sigma parameters instead of Mass for scale */

double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias); // as above, but this version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate integral version of eq. 2 in Barkana & Loeb 2008)
double dicke(double z); //calculates the dicke growth function at redshift z
double ddickedz(double z); /* Redshift derivative of the growth function at z */
double ddickedt(double z); /* Time derivative of the growth function at z */
double sigma_z0(double M); //calculates sigma at z=0 (no dicke)
double dsigmasqdm_z0(double M); //calculates d(sigma^2)/dm at z=0 (i.e. does not include dicke growth)
double TFmdm(double k); //Eisenstien & Hu power spectrum transfer function
void TFset_parameters();
float get_R_c();  // returns R_CUTOFF
double get_M_min_ion(float z);
/***************************************/

/* Returns the minimum source mass for ionizing sources, according to user specifications */
double get_M_min_ion(float z){
  double MMIN;
  if (ION_M_MIN < 0){ // use the virial temperature for Mmin
    if (ION_Tvir_MIN < 9.99999e3) // neutral IGM
      MMIN = TtoM(z, ION_Tvir_MIN, 1.22);
    else // ionized IGM
      MMIN = TtoM(z, ION_Tvir_MIN, 0.6);
  }
  else if (ION_Tvir_MIN < 0){ // use the mass
    MMIN = ION_M_MIN;
  }
  else{
    fprintf(stderr, "You have to \"turn-off\" either the ION_M_MIN or \
                     the ION_Tvir_MIN option in ANAL_PARAMS.H\nAborting...\n");
    return -1;
  }
  // check for WDM
  if (P_CUTOFF && ( MMIN < M_J_WDM()))
    MMIN = M_J_WDM();
  //  printf("Mmin is %e\n", MMIN);
  return MMIN;
}

/* Returns the minimum source mass for x-ray sources, according to user specifications */
double get_M_min_xray(float z){
  double MMIN;
  if (X_RAY_Tvir_MIN < 9.99999e3) //neutral IGM
    MMIN = TtoM(z, X_RAY_Tvir_MIN, 1.22);
  else // ionized IGM
    MMIN = TtoM(z, X_RAY_Tvir_MIN, 0.6);

  // check for WDM
  if (P_CUTOFF && ( MMIN < M_J_WDM()))
    MMIN = M_J_WDM();
  //  printf("Mmin is %e\n", MMIN);
  return MMIN;
}


/* returns the "effective Jeans mass" in Msun
   corresponding to the gas analog of WDM ; eq. 10 in Barkana+ 2001 */
double M_J_WDM(){
  double z_eq, fudge=60;
  if (!P_CUTOFF) 
    return 0;
  z_eq = 3600*(OMm-OMb)*hlittle*hlittle/0.15;
  return fudge*3.06e8 * (1.5/g_x) * sqrt((OMm-OMb)*hlittle*hlittle/0.15) * pow(M_WDM, -4) * pow(z_eq/3000.0, 1.5);
}


/* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double deltolindel(float del, float z){
  float onepdel = 1.0+del;
  return ( 1.68647 - 1.35*pow(onepdel,-2/3.0) + 0.78785*pow(onepdel,-0.58661) - 1.12431*pow(onepdel,-0.5) )/dicke(z);
}

/* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double lindeltodel(float lindel, float z){
  float prev_lindelguess, delcrit, delguess;
  float lindelguess, delmin, delmax, epsilon = 1.0e-7;

  // set the critical density corresponding to virialization
  // this will be maximum allowed del
  delcrit = Deltac_nonlinear(z)*rho_critz(z)/(OMm*RHOcrit*pow(1+z, 3)) - 1;

  delmin = -1;
  delmax = 500;
  prev_lindelguess = -1e10;
  while (1){
    delguess = 0.5*(delmax+delmin);
    lindelguess = deltolindel(delguess, z);
    //fprintf(stderr, "%e\t%e\n", delmin, delmax);
	      //    fprintf(stderr, "%e\t%e\t%e\n\n", delguess, lindelguess, lindel);
    if ((fabs((lindelguess-lindel)/lindel) < epsilon ) ||
	(fabs(lindelguess-lindel) < epsilon ) || 
	(fabs(prev_lindelguess - lindelguess) < TINY ))// close enough, or resolution loop
      return delguess;

    if (lindelguess > lindel)
      delmax = delguess;
    else
      delmin = delguess;

    // check if we are above delcrit (see above)
    if (delmin > delcrit){
      //      printf("exced max at lindel=%e\n", lindel);
      return delcrit;
    }

    prev_lindelguess = lindelguess;
  }
}


/* R in Mpc, M in Msun */
double RtoM(double R){
  // set M according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return (4.0/3.0)*PI*pow(R,3)*(OMm*RHOcrit);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow(2*PI, 1.5) * OMm*RHOcrit * pow(R, 3);
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}

/* R in Mpc, M in Msun */
double MtoR(double M){
  // set R according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return pow(3*M/(4*PI*OMm*RHOcrit), 1.0/3.0);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow( M/(pow(2*PI, 1.5) * OMm * RHOcrit), 1.0/3.0 );
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}


/* equation (5) from jenkis et al. (2001) */
double f_jenkins(float del, double sigsq){
  if (del < 0){  fprintf(stderr, "ERROR:  In function f_jenkins del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }

  //  fprintf(stderr, "%f\t%f\n", del, sqrt(sigsq));
  return sqrt(2/PI) * del/sqrt(sigsq) * pow(E, -0.5*del*del/sigsq);
}


float get_R_c(){
  return R_CUTOFF;
}

/* sheth correction to delta crit */
double sheth_delc(double del, double sig){
  return sqrt(SHETH_a)*del*(1 + SHETH_b*pow(sig*sig/(SHETH_a*del*del), SHETH_c));
}


/* dnbiasdM */
double dnbiasdM(double M, float z, double M_o, float del_o){
  double sigsq, del, sig_one, sig_o;

  if ((M_o-M) < TINY){
    fprintf(stderr, "WARNING:  In function dnbiasdM: M must be less than M_o!\nAborting...\n");
    return -1;
  }
  del = Deltac/dicke(z) - del_o;
  if (del < 0){  fprintf(stderr, "ERROR:  In function dnbiasdM: del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }
  sig_o = sigma_z0(M_o);
  sig_one = sigma_z0(M);
  sigsq = sig_one*sig_one - sig_o*sig_o;
  return -(RHOcrit*OMm)/M /sqrt(2*PI) *del*pow(sigsq,-1.5)*pow(E, -0.5*del*del/sigsq)*dsigmasqdm_z0(M);
}


/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function with Sheth-Torman correction for ellipsoidal collapse at 
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Sheth, Mo, Torman 2001
*/
double dNdM_st(double z, double M){
  double sigma, dsigmadm, nuhat, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * dicke_growth*dicke_growth/(2.0*sigma);
  nuhat = sqrt(SHETH_a) * Deltac / sigma;

  return (-OMm*RHOcrit/M) * (dsigmadm/sigma) * sqrt(2/PI)*SHETH_A * (1+ pow(nuhat, -2*SHETH_p)) * nuhat * pow(E, -nuhat*nuhat/2.0);
}



/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function at 
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Padmanabhan, pg. 214
*/
double dNdM(double z, double M){
  double sigma, dsigmadm, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * (dicke_growth*dicke_growth/(2*sigma));

  return (-OMm*RHOcrit/M) * sqrt(2/PI) * (Deltac/(sigma*sigma)) * dsigmadm * pow(E, -(Deltac*Deltac)/(2*sigma*sigma));
}


/*
  FUNCTION FgtrM_st(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
  Uses Sheth-Torman correction
*/
double dFdlnM_st (double lnM, void *params){
  double z = *(double *)params;
  double M = exp(lnM);
  return dNdM_st(z, M) * M * M;
}
double FgtrM_st(double z, double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = 0.001; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);

  F.function = &dFdlnM_st;
  F.params = &z;
  lower_limit = log(M);
  upper_limit = log(FMAX(1e16, M*100));
			   
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);
}


/*
  FUNCTION FgtrM(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
*/
double FgtrM(double z, double M){
  double del, sig;

  del = Deltac/dicke(z); //regular spherical collapse delta
  sig = sigma_z0(M);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
 */
double FgtrM_bias(double z, double M, double del_bias, double sig_bias){
  double del, sig, sigsmallR;

  sigsmallR = sigma_z0(M);

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

   return splined_erfc(del / (sqrt(2)*sig));
}


/* Uses sigma parameters instead of Mass for scale */
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias){
  double del, sig;

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "local_FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  Calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias.
  This version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate
  integral version of eq. 2 in Barkana & Loeb 2008)
 */
double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias){
  return FgtrM_st(z, M) / FgtrM(z, M) * FgtrM_bias(z, M, del_bias, sig_bias);
}


/*
  FUNCTION dicke(z)
  Computes the dicke growth function at redshift z, i.e. the z dependance part of sigma

  References: Peebles, "Large-Scale...", pg.53 (eq. 11.16). Includes omega<=1
  Nonzero Lambda case from Liddle et al, astro-ph/9512102, eqs. 6-8.
  and quintessence case from Wang et al, astro-ph/9804015

  Normalized to dicke(z=0)=1
*/
double dicke(double z){
  double omegaM_z, dick_z, dick_0, x, x_0;
  double tiny = 1e-4;

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return 1.0/(1.0+z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    dick_z = 2.5*omegaM_z / ( 1.0/70.0 + omegaM_z*(209-omegaM_z)/140.0 + pow(omegaM_z, 4.0/7.0) );
    dick_0 = 2.5*OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );
    return dick_z / (dick_0 * (1.0+z));
  }
  else if ( (OMtot < (1+tiny)) && (fabs(OMl) < tiny) ){ //open, zero lambda case (peebles, pg. 53)
    x_0 = 1.0/(OMm+0.0) - 1.0;
    dick_0 = 1 + 3.0/x_0 + 3*log(sqrt(1+x_0)-sqrt(x_0))*sqrt(1+x_0)/pow(x_0,1.5);
    x = fabs(1.0/(OMm+0.0) - 1.0) / (1+z);
    dick_z = 1 + 3.0/x + 3*log(sqrt(1+x)-sqrt(x))*sqrt(1+x)/pow(x,1.5);
    return dick_z/dick_0;
  }
  else if ( (OMl > (-tiny)) && (fabs(OMtot-1.0) < tiny) && (fabs(wl+1) > tiny) ){
    fprintf(stderr, "IN WANG\n");
    return -1;
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}


/* redshift derivative of the growth function at z */
double ddicke_dz(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;

   return (dicke(z+dz)-dicke(z))/dz;
}

/* Time derivative of the growth function at z */
double ddickedt(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;
  double tiny = 1e-4;

   return (dicke(z+dz)-dicke(z))/dz/dtdz(z); // lazy non-analytic form getting

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return -pow(1+z,-2)/dtdz(z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    domegaMdz = omegaM_z*3/(1+z) - OMm*pow(1+z,3)*pow(OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4), -2) * (3*OMm*(1+z)*(1+z) + 4*OMr*pow(1+z,3));
    dick_0 = OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );

    ddickdz = (domegaMdz/(1+z)) * (1.0/70.0*pow(omegaM_z,-2) + 1.0/140.0 + 3.0/7.0*pow(omegaM_z, -10.0/3.0)) * pow(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0) , -2);
    ddickdz -= pow(1+z,-2)/(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0));

    return ddickdz / dick_0 / dtdz(z);
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}

/*
  FUNCTION sigma_z0(M)
  Returns the standard deviation of the normalized, density excess (delta(x)) field,
  smoothed on the comoving scale of M (see filter definitions for M<->R conversion).
  The sigma is evaluated at z=0, with the time evolution contained in the dicke(z) factor,
  i.e. sigma(M,z) = sigma_z0(m) * dicke(z)

  normalized so that sigma_z0(M->8/h Mpc) = SIGMA8 in ../Parameter_files/COSMOLOGY.H

  NOTE: volume is normalized to = 1, so this is equvalent to the mass standard deviation

  M is in solar masses

  References: Padmanabhan, pg. 210, eq. 5.107
*/
double dsigma_dk(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, kR;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  // NOTE: only use top hat for SIGMA8 normalization
  kR = k*R;
  if ( (FILTER == 0) || (sigma_norm < 0) ){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0;} // w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  return k*k*p*w*w;
}
double sigma_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;

  R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return sigma_norm * sqrt(result);
}


/*
  Returns the value of the linear power spectrum DENSITY (i.e. <|delta_k|^2>/V)
  at a given k mode linearly extrapolated to z=0
*/
double power_in_k(double k){
  double p, T, gamma, q, aa, bb, cc;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  return p*TWOPI*PI*sigma_norm*sigma_norm;
}


/*
  FUNCTION dsigmasqdm_z0(M)
  returns  d/dm (sigma^2) (see function sigma), in units of Msun^-1
*/
double dsigmasq_dm(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, dwdr, drdm, kR;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu ApJ, 1999, 511, 5
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / (hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + pow(bb*k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  kR = k * R;
  if (FILTER == 0){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0; }// w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}

    // now do d(w^2)/dm = 2 w dw/dr dr/dm
    if ( (kR) < 1.0e-10 ){  dwdr = 0;}
    else{ dwdr = 9*cos(kR)*k/pow(kR,3) + 3*sin(kR)*(1 - 3/(kR*kR))/(kR*R);}
	    //3*k*( 3*cos(kR)/pow(kR,3) + sin(kR)*(-3*pow(kR, -4) + 1/(kR*kR)) );}
      //     dwdr = -1e8 * k / (R*1e3);
    drdm = 1.0 / (4.0*PI * OMm*RHOcrit * R*R);
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
    dwdr = - k*kR * w;
    drdm = 1.0 / (pow(2*PI, 1.5) * OMm*RHOcrit * 3*R*R);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  //  printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", k, R, p, w, dwdr, drdm, dsigmadk[1]);
  return k*k*p*2*w*dwdr*drdm * d2fact;
}
double dsigmasqdm_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;

  R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);
  d2fact = M*10000/sigma_z0(M);

  F.function = &dsigmasq_dm;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return sigma_norm * sigma_norm * result /d2fact;
}



/*
  FUNCTION TFmdm is the power spectrum transfer function from Eisenstein & Hu ApJ, 1999, 511, 5
*/
double TFmdm(double k){
  double q, gamma_eff, q_eff, TF_m, q_nu;

  q = k*pow(theta_cmb,2)/omhh;
  gamma_eff=sqrt(alpha_nu) + (1.0-sqrt(alpha_nu))/(1.0+pow(0.43*k*sound_horizon, 4));     
  q_eff = q/gamma_eff;
  TF_m= log(E+1.84*beta_c*sqrt(alpha_nu)*q_eff);
  TF_m /= TF_m + pow(q_eff,2) * (14.4 + 325.0/(1.0+60.5*pow(q_eff,1.11)));
  q_nu = 3.92*q/sqrt(f_nu/N_nu);
  TF_m *= 1.0 + (1.2*pow(f_nu,0.64)*pow(N_nu,0.3+0.6*f_nu)) / 
    (pow(q_nu,-1.6)+pow(q_nu,0.8));

  //   printf("%f  %e  %f  %f  %f  %f\n",omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu);
  // printf("%f  %f  %f  %f\n", beta_c,sound_horizon,theta_cmb,z_equality);
  //printf("%f  %e  %f  %f  %f\n\n",q, k, gamma_eff, q_nu, TF_m);


  return TF_m;
}


void TFset_parameters(){
  double z_drag, R_drag, R_equality, p_c, p_cb, f_c, f_cb, f_nub, k_equality;

  z_equality = 25000*omhh*pow(theta_cmb, -4) - 1.0;
  k_equality = 0.0746*omhh/(theta_cmb*theta_cmb);

  z_drag = 0.313*pow(omhh,-0.419) * (1 + 0.607*pow(omhh, 0.674));
  z_drag = 1 + z_drag*pow(OMb*hlittle*hlittle, 0.238*pow(omhh, 0.223));
  z_drag *= 1291 * pow(omhh, 0.251) / (1 + 0.659*pow(omhh, 0.828));

  y_d = (1 + z_equality) / (1.0 + z_drag);

  R_drag = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_drag);
  R_equality = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_equality);

  sound_horizon = 2.0/3.0/k_equality * sqrt(6.0/R_equality) * 
    log( (sqrt(1+R_drag) + sqrt(R_drag+R_equality)) / (1.0 + sqrt(R_equality)) );

  p_c = -(5 - sqrt(1 + 24*(1 - f_nu-f_baryon)))/4.0;
  p_cb = -(5 - sqrt(1 + 24*(1 - f_nu)))/4.0;
  f_c = 1 - f_nu - f_baryon;
  f_cb = 1 - f_nu;
  f_nub = f_nu+f_baryon;

  alpha_nu = (f_c/f_cb) * (2*(p_c+p_cb)+5)/(4*p_cb+5.0);
  alpha_nu *= 1 - 0.553*f_nub+0.126*pow(f_nub,3);
  alpha_nu /= 1-0.193*sqrt(f_nu)+0.169*f_nu;
  alpha_nu *= pow(1+y_d, p_c-p_cb);
  alpha_nu *= 1+ (p_cb-p_c)/2.0 * (1.0+1.0/(4.0*p_c+3.0)/(4.0*p_cb+7.0))/(1.0+y_d);
  beta_c = 1.0/(1.0-0.949*f_nub);
}




double init_ps(){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;
  int i;
  double x;

  // Set cuttoff scale for WDM (eq. 4 in Barkana et al. 2001) in comoving Mpc
  R_CUTOFF = 0.201*pow((OMm-OMb)*hlittle*hlittle/0.15, 0.15)*pow(g_x/1.5, -0.29)*pow(M_WDM, -1.15);

//  fprintf(stderr, "For M_DM = %.2e keV, R_CUTOFF is: %.2e comoving Mpc\n", M_WDM, R_CUTOFF);
  if (!P_CUTOFF)
//    fprintf(stderr, "But you have selected CDM, so this is ignored\n");

  omhh = OMm*hlittle*hlittle;
  theta_cmb = T_cmb / 2.7;

  // Translate Parameters into forms GLOBALVARIABLES form
  f_nu = OMn/OMm;
  f_baryon = OMb/OMm;
  if (f_nu < TINY) f_nu = 1e-10;
  if (f_baryon < TINY) f_baryon = 1e-10;


  TFset_parameters();

  sigma_norm = -1;

  R = 8.0/hlittle;
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);

  sigma_norm = SIGMA8/sqrt(result); //takes care of volume factor


  /* initialize the lookup table for erfc */
  /*
  for (i=0; i<=ERFC_NPTS; i++){
    erfc_params[i] = i*ERFC_PARAM_DELTA;
    log_erfc_table[i] = log(erfcc(erfc_params[i]));
  }
  // Set up spline table
  erfc_acc   = gsl_interp_accel_alloc ();
  erfc_spline  = gsl_spline_alloc (gsl_interp_cspline, ERFC_NPTS);
  gsl_spline_init(erfc_spline, erfc_params, log_erfc_table, ERFC_NPTS);
  */

  return R_CUTOFF;
}

void free_ps(){
  /*    gsl_spline_free (erfc_spline);
    gsl_interp_accel_free(erfc_acc);
  */
  return;
}

double splined_erfc(double x){
  if (x < 0){
    //    fprintf(stderr, "WARNING: Negative value %e passed to splined_erfc. Returning 1\n", x);
    return 1;
  }
  return erfcc(x); // the interpolation below doesn't seem to be stable in Ts.c
  if (x > ERFC_PARAM_DELTA*(ERFC_NPTS-1))
    return erfcc(x);
  else
    return exp(gsl_spline_eval(erfc_spline, x, erfc_acc));
}

#endif
