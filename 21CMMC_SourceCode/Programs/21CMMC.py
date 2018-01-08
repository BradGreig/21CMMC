import os
import numpy
import math
from scipy import interpolate
from decimal import *
import string
import pickle
import time

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')
SIXPLACES = Decimal(10) ** -6       # same as Decimal('0.000001')

from CosmoHammer_21CMMC.sampler.CosmoHammerSampler import CosmoHammerSampler
from CosmoHammer_21CMMC.likelihood.chain.LikelihoodComputationChain import LikelihoodComputationChain
from CosmoHammer_21CMMC.likelihood.module.Likelihood21cmFast import Likelihood21cmFast_multiz

from CosmoHammer_21CMMC.likelihood.module.Likelihood21cmFast import McGreer_Redshift
from CosmoHammer_21CMMC.likelihood.module.Likelihood21cmFast import QSO_Redshift

if __name__ == '__main__':


	# New version of 21CMMC allows for the full computation of the spin temperature fluctuations during the X-ray heating epoch. 21CMMC now includes a modified and 
	# streamlined version of Ts.c from 21cmFAST. 

	#### NOTE ####
	# At present the full spin temperature calculation does not support a mass dependent ionising efficiency scaling. Therefore, only works for IncludeAlpha = False.
	# Additionally, it assumes the same minimum ionising mass for the ionising and X-ray sources (ION_TVIR_MIN = X_RAY_TVIR_MIN).
	# Eventually these assumptions will be removed...



	################### Setting up the flag options for the main computation and requisite interpolation tables ####################################



	# Whether or not one wants the inital conditions to be varied within the MCMC (needs to be True if one wants to jointly sample the cosmological parameters)
	#### NOTE ####
	# If this is set to true, none of the interpolation tables can be used as they are valid only for a single set of initial conditions.
	GenerateNewICs = False

	# Create the list of cosmological values to vary. Parameter names entered here must match the naming convention below, otherwise they will not be detected and thus varied
	# Fiducial values and parameter bounds for the cosmological parameters can be changed below
#	CosmologyToVary = ['SIGMA_8','littleh','OMEGA_M','OMEGA_b','NS']
	CosmologyToVary = []

	# Performs the full evolution (Ts.c) of the IGM during reionisation and heating epoch. Setting to false reverts to saturated spin temperature limit (Ts >> Tcmb).
	Include_Ts_fluc = True

	# If the full spin temperature computation is to be performed, a redshift must be provided to which to perform the evolution down to.
	TsCalc_z = 6.0

	# Decide whether to use light-cone boxes or co-eval boxes
	# Note that the light-cone can only be generated along the z-direction (21cmFAST could do any arbitrary direction, this only does the z-direction). Should be 
	# trivial if one wants to add support for light-cones along any direction.
	IncludeLightCone = False

	# Use an interpolation table for the full box collapsed fraction for the computation of the IGM spin temperature. 
	UseFcollTable = False
	# The full Ts.c calculation can be ~20% faster if an Interpolation table for the collapsed fraction is used. This collapsed fraction is
	# averaged across the full box (necessary for Ts.c) as a function of the minimum halo mass for the full range set for this parameter below.
	# This table can be created by calling CreateFcollTable from the command line, and providing the necessary arguments (See CreateFcollTable.c for further details).
	### NOTE ### 
	# The advantage of this is that it only has to be done once per cosmology/box size and mininum halo mass range. If any of these are changed this table
	# must be regenerated
	### NOTE ### 
	# This is only really worth doing for boxes of size HII_DIM >= 128.			
	### NOTE ### 
	# Be careful that you are using a valid interpolation table (i.e cosmology etc. hasn't changed)	

	### NOTE ###
	# Always double check the code against the version without the interpolation table to ensure the resolution (binning) of the interpolation table
	# is sufficiently high to obtain sufficiently high accuracy results

	# Use an interpolation table for the full box as a function of smoothing radius for the computation of the IGM neutral fraction (find_HII_bubbles part of the code)
	USE_IONISATION_FCOLL_TABLE = False	
	# This is generated using 'Create_ionisaton_fcoll_tables.py'. What this does is compute the averaged collapsed fraction for each light-cone redshift, as a function
	# of smoothing radius and Tvir_min. Setting this flag can boost code efficiency by ~15-20%. This must be generated for the same upper limit for the mean free path
	# (i.e. UpperBound_MFP) and for the same lower and upper limits for the virial temperature of the haloes (X_RAY_TVIR_LB and	X_RAY_TVIR_UB)
	### NOTE ###
	# This must be called outside of this code. See 'Create_ionisaton_fcoll_tables.py' for more details.	
	### NOTE ### 
	# This is different to 'CreateFcollTable' above

	# Whether to include inhomogeneous recombinations in the computation
	USE_INHOMO_RECO = True
	# This will enable inhomogeneous recombinations to be used, as included in the latest version of 21cmFAST
	# This uses the Sobacchi & Mesinger (2014) approach for computing the sub-grid recombinations
	### NOTE ###
	# In setting this as true, the mean free path parameter (R_MFP) is no longer a free parameter. R_MFP becomes fixed at 50 Mpc, as is default in 21cmFAST
	# To actually change this default choice of 50 Mpc, the user needs to change the value for INHOMO_RECO_R_BUBBLE_MAX in drive_21cmMC_streamlined.c as this is 
	# where the value is hard coded (needs to be there for when drive_21cmMC_streamlined is run independent of the MCMC sampler)
	### NOTE ###
	# When USE_INHOMO_RECO is True, then the code reverts to performing find_HII_bubbles for all the redshifts sampled by the spin temperature calculation (Ts.c)
	# This is because to perform the inhomogeneous recombinations, the previous timesteps need to be stored to accurately track the sub-grid recombinations


	# Whether to include line of sight (z-direction only) redshift space distortions (RSDs)
	INCLUDE_SUBCELL_RSDS = False
	# This is applied along the line of sight (z-direction)
	# The preamble of drive_21cmMC_streamlined.c contains a few other parameters that are associated with RSDs. Change within drive_21cmMC_streamlined.c

	# The number of sub-cells to be used for determining the line-of-sight redshift space distortions
	N_RSD_SUBCELLS = 20
	
	# Line-of-sight direction. 0: x 1: y 2: z
	### NOTE ###
	# I do not know if it works for the x or y directions, I have not checked (I think it should though). I always work in the z-direction.
	LOS_direction = 2

	USE_GLOBAL_SIGNAL = False
	# Set to true if one wants to use the global signal rather than the 21cm PS to parameter sampling. Note, the mock observations will need to be changed accordingly
	### NOTE ###
	# If one wants to use the global signal, set IncludeLightCone = True. Setting this to true, will output the global signal at the redshift sampling of the spin temperature
	# calculation, i.e. a large number of data points with fine delta-z (frequency) resolution
	# If the computation is too long, one can set IncludeLightCone = False, and provide the desired redshift sampling through the "Redshift" list below. Note however, that a 
	# fine enough sampling should be provided to ensure reasonable interpolation of the global signal curve. I have arbitrarily set a requirement of a minimum of 20 redshift
	# (frequency) bins to measure and interpolate the global signal (in the absence of setting IncludeLightCone = True).
	### NOTE ###
	# The global signal is interpolated, therefore it is important that if IncludeLightCone = False, that the redshifts selected in the lists Redshift or Redshifts_For_Prior
	# include the endpoints of the redshift (frequency) range of the global signal measurement (i.e. cannot guarantee the robustness of the extrapolation). Additionally, be wary 
	# that the chosen TsCalc_z (redshift of the spin temperature calculation) corresponds to a redshift (frequency) consistent with the global signal measurement. For example,
	# the default TsCalc_z = 6.0 corresponds basically to 200 MHz (in-built (default) Z_HEAT_MAX is z = 35 which corresponds to 40 MHz)	

	# Define whether a fixed error on the 21cm global signal is to be used, or whether to read from file
	GLOBAL_SIGNAL_FIXED_ERROR = False


	# Setting this to true will keep all relevant statistical data (i.e. tau, xH vs etc., PS vs k at all redshift etc.)
	# Separating the accepted/rejected points from the MCMC output can be done in post-processing (a separate script is provided to do so "ReadAllData.py"). 
	# It was a bit too unwieldly to do internally, so I opted for externally dealing with separating the data.
	KEEP_ALL_DATA = False



	################### Setting up variables for performing the full spin temperature calculation (i.e. Ts.c) ####################################



	#### NOTE: Here, in 21CMMC some of the variables are different to the usual Ts.c variables. However this will change with an updated Ts.c coming soon for 21cmFAST ####	
	# In this version, I have replaced the old ZETA_X with a soft-band luminosity L_X. The soft band luminosity now determines the normalisation, and is available as a free
	# parameter. This is defined as the 2keV luminosity, hence NU_X_BAND_MAX = 2000.00.
	# The full X-ray emissivity calculation is evaluated over the full 10keV band, and hence NU_X_MAX = 10000.0.
	#### NOTE: these are not varied, but are fixed. No need to change either of these ####
	NU_X_BAND_MAX = 2000.0
	NU_X_MAX = 10000.0

	# For the creation of the total collapsed fraction interpolation table, it is important to know the lower and upper bounds for the full range of X_RAY_TVIR_MIN. Can 
	# set these quantities here. X_RAY_TVIR_LB (lower bound) and X_RAY_TVIR_UB (upper bound). Note, these are log10 quantities of the lower and upper limits.
	X_RAY_TVIR_LB = 4.0
	X_RAY_TVIR_UB = 6.0

	# Setting IncludeAlpha = True, results in the sampling of a mass dependent ionising efficiency.
	# This mass dependent ionising efficiency is defined as (Mass / M_vir)**(Alpha)
	# - M_vir is the mass corresponding to the Virial Temperature (T_Vir) evaluated at the respective redshift (T_vir is defined to be redshift independent)
	# - Alpha is a the power law. Alpha = 0 corresponds to a mass independent ionising efficiency, corresponding to the default 3 parameter reionisation model
	# IncludeAlpha = True allows for a 4 parameter model to be sampled.
	# *** NOTE *** As yet there is no support for the 4 parameter model. This needs to be added at some point down the line. I don't know what you will get if
	# you set this to 'True' at the present time! (I think it works if Include_Ts_fluc is set to false, i.e. only co-eval boxes in the saturated spin temperature limit)
	IncludeAlpha = False 

	# Reionisation redshifts for the multi-z 21cm Fast "observations"
	# Need to make sure that these boxes exist in the "Boxes" folder. If not, please generate new boxes
	# Shifted the range to lower-z, motivated by the lower Planck (2016) values of tau
	# Expanded the default range to allow for an approximate computation of tau, to be applied as a prior

	####### Note: If spin temperature fluctuations are to be performed, select redshifts that **will** be sampled by the spin temperature algorithm. This will require
	# running drive_21cmMC_streamlined once, and outputting the variable 'zp' to see the redshift sampling ########
	
	Redshift = []

	# This list allows the user to add additional redshifts to the list (co-eval boxes only) to improve the sampling for any of the priors (is not used for the likelihood calucation)
	# Note: Adding any additional redshifts adds to the computation time, so this will make the code slower
	Redshifts_For_Prior = []

	# If the light-cone is being directly sampled, it outputs across the full redshift range, so don't need to pass it a redshift. Will populate the list 'Redshift'
	# with the filenames of the mock observations constructing the light-cone to ensure the correct number of 21cm PS are used for the likelihood
	### NOTE ### Redshifts added here must always be added in increasing order.
	if IncludeLightCone is False:

		Redshift = ['6.429094', '7.041489', '7.533692', '8.610322']

#		Redshift = ['6.429094', '7.202319', '8.237142', '9.402521', '10.714930', '12.456770', '14.457590', '17.111031','20.219959','24.359810']
#		Redshift = ['6.000594', '7.041489', '8.056021', '8.998578', '10.949220', '13.000420', '15.082080', '17.111023']	
#		Redshift = ['6.000594', '7.041489', '8.056021', '8.998578']	
#		Redshift = ['6.000594']

		Redshifts_For_Prior = []

	
	### NOTE ###
	# if Include_Ts_fluc = True the redshifts listed above must match with the redshifts sampled by the spin temperature algorithm. If not, then the code will fail.
	# Generally speaking, the sampling should be fine enough that one redshift will be close enough to the corresponding redshift of interest.
	# If not, and it is important for accuracy purposes, one can lower ZPRIME_STEP_FACTOR in HEAT_PARAMS.H
	# Can determine the redshift sampling by running the test instance of the driver (./drive_21cmMC_streamlined 1.000000 1.000000 1 1 0 6.0), using the provided Walker file.


	################### Enabling the addition of some observational priors ####################################


	# Added some basic support for observational priors on the IGM neutral fraction, or reionisation history from Planck tau
	# (can use the form of these priors to include your own priors on individual parameters, or for global quantites)
	# Provide three possible observational priors:
	# 1) Constraints on tau from Planck (2016) (https://arxiv.org/abs/1605.03507)
	# 2) Limit on the IGM neutral fraction at z = 5.9, from dark pixels by I. McGreer et al. (2015) (http://adsabs.harvard.edu/abs/2015MNRAS.447..499M)
	# 3) Constraints on the IGM neutral fraction at z = 7.1 from the IGM damping wing of ULASJ1120+0641 Greig et al (2016) (http://arxiv.org/abs/1606.00441)
	
	# NOTE: Constraints are set in CosmoHammer_21CMMC/likelihood/module/Likelihood21cmFast.py
	# Go here to change any of the prior behaviour, or to add additional priors not provided below.

	# NOTE: If any priors are to be added, it is recommended that at least three redshifts are considered. 
	# PROVIDING LESS THAN THREE REDSHIFTS WILL CAUSE THE USE OF PRIORS TO BE TURNED OFF (i.e PRIORS WILL NOT BE INCLUDED FOR LESS THAN THREE REDSHIFTS)
	# These priors require interpolation/extrapolation of the reionisation history, therefore, to aid accuracy at least three should be selected.

	# NOTE: If the redshifts selected by the user ** is ** at the redshift of the prior, then interpolation is not required, and the prior will be used!

	# 1) The Planck prior is modelled as a Gaussian: tau = 0.058 \pm 0.012
	IncludePlanck = False
	# 2) The McGreer et al. prior is a upper limit on the IGM neutral fraction at 5.9
	# Modelled as a flat, unity prior at x_HI <= 0.06, and a one sided Gaussian at x_HI > 0.06 ( Gaussian of mean 0.06 and one sigma of 0.05 )
	IncludeMcGreer = False
	# 3) The Greig et al. prior is computed directly from the PDF of the IGM neutral fraction for the Small HII reionisation simulation
	IncludeGreig = False
	
	# A check to see if additional redshifts have been added to the computation that are not necessary
	if IncludePlanck is False and IncludeMcGreer is False and IncludeGreig is False and len(Redshifts_For_Prior) > 0:
		print 'Additional redshifts have been added (in Redshifts_For_Prior) that are not necessary'
		print 'Setting Redshifts_For_Prior to an empty list'

		Redshifts_For_Prior = []

	

	PriorLegend = dict()
	PriorLegend['PlanckPrior'] = IncludePlanck
	PriorLegend['McGreerPrior'] = IncludeMcGreer
	PriorLegend['GreigPrior'] = IncludeGreig
	
	# If the QSO damping wing constraint is set, need to read in the PDF data from the provided text file.
	NFVals_QSODamping = []
	PDFVals_QSODamping = []

	if IncludeGreig is True:

		with open('PriorData/NeutralFractionsForPDF.out','rb') as handle:
			NFVals_QSODamping =	pickle.loads(handle.read())

		with open('PriorData/NeutralFractionPDF_SmallHII.out','rb') as handle:
			PDFVals_QSODamping = pickle.loads(handle.read())

		# Normalising the PDF to have a peak probability of unity (consistent with how other priors are treated)
		# Ultimately, this step does not matter
		normalisation = numpy.amax(PDFVals_QSODamping)

		PDFVals_QSODamping = PDFVals_QSODamping/normalisation



	###### Read in the data for the mock observation ######



	# Read in the mock 21cm PS observation. Read in both k and the dimensionless PS. These are needed for performing the chi^2 statistic for the likelihood
	# Note: To calculate the likelihood statistic a spline is performed for each of the mock PS, simulated PS and the Error PS	
	multi_z_mockobs_k = []
	multi_z_mockobs_PS = []

	if USE_GLOBAL_SIGNAL is True:
	
		MockObsFileName = 'FaintGalaxies_GlobalSignal'
#		MockObsFileName = 'BrightGalaxies_GlobalSignal'
		ModelName = 'FaintGalaxies'
#		ModelName = 'BrightGalaxies'

		### NOTE ###
		# The noise and mock observation of the global signal must be provided in order of decreasing redshift (it'll convert to frequency, which will then be
		# in increasing order for the spline interpolations)
		mockobs_k_values = numpy.loadtxt('MockObs/%s/GlobalSignal/%s.txt'%(ModelName,MockObsFileName), usecols=(0,))
		mockobs_PS_values = numpy.loadtxt('MockObs/%s/GlobalSignal/%s.txt'%(ModelName,MockObsFileName), usecols=(2,))
		
		multi_z_mockobs_k.append(mockobs_k_values)
		multi_z_mockobs_PS.append(mockobs_PS_values)

	else:				

		if IncludeLightCone is True:

			# For the light-cone version of the code, use a text file containing the 21cm PS from the light-cones.	
			# *** NOTE *** This text file should contain the number of light-cone 21cm in ** increasing ** redshift order (i.e. lowest-z to highest-z)
			# *** NOTE *** Only need to include up until the redshift of interest. The likelihood will only be computed for this number of 21cm PS (the c-code creates
			# its own text file in increasing order, and should be sampled along the L.o.S the same as these mock observations. i.e. the Delta-z or frequency bandwidth
			# between the mock observations and sampled 21cm PS should be the same)

			# For the light-cone version of the code, use a text file containing the 21cm PS from the light-cones.		
			MockObsFileName = 'LightCone21cmPS_FaintGalaxies_600Mpc_400'
		
			# Note here, we are populating the list 'Redshift' with the filenames. The length of this is needed for ensuring the correct number of 21cm PS are used for the likelihood
			# Re-using the same list filename means less conditions further down this script. The likelihood correctly accounts for it with the 'IncludeLightCone' flag.
			filename = open('MockObs/%s.txt'%(MockObsFileName), 'r') 
			Redshift = [line.rstrip('\n') for line in filename]

			for i in range(len(Redshift)):
				mockobs_k_values = numpy.loadtxt('MockObs/%s'%(Redshift[i]), usecols=(0,))
				mockobs_PS_values = numpy.loadtxt('MockObs/%s'%(Redshift[i]), usecols=(1,))

				multi_z_mockobs_k.append(mockobs_k_values)
				multi_z_mockobs_PS.append(mockobs_PS_values)		

		else:

			### NOTE ###
			# If Include_Ts_fluc is set, the user must ensure that the co-eval redshift to be sampled (set by the Redshift list above) is to be sampled by the code.			

#			MockObsFileName = 'MockObs_FaintGalaxies_PS_600Mpc'
			MockObsFileName = 'MockObs_PS_200Mpc_EOS_FaintGalaxies'
#			MockObsFileName = 'MockObs_BrightGalaxies_PS_600Mpc'
#			ModelName = 'FaintGalaxies'
			ModelName = 'EOS_FaintGalaxies'
#			ModelName = 'BrightGalaxies'
#			BoxType = 'Co-Eval'
			BoxType = 'Co-Eval'			

			for i in range(len(Redshift)):		
				mockobs_k_values = numpy.loadtxt('MockObs/%s/%s/%s_z%s.txt'%(ModelName,BoxType,MockObsFileName,Redshift[i]), usecols=(0,))
				mockobs_PS_values = numpy.loadtxt('MockObs/%s/%s/%s_z%s.txt'%(ModelName,BoxType,MockObsFileName,Redshift[i]), usecols=(1,))

				multi_z_mockobs_k.append(mockobs_k_values)
				multi_z_mockobs_PS.append(mockobs_PS_values)

	multi_z_mockobs_k = numpy.array(multi_z_mockobs_k)
	multi_z_mockobs_PS = numpy.array(multi_z_mockobs_PS)



	###### Read in the data for the telescope sensitivites ######



	# Set for the desired telescope ['SKA_halveddipoles_compact', 'HERA331'] corresponding to the file structure in "NoiseData"
#	Telescope_Name = 'SKA'
	Telescope_Name = 'HERA331'
#	Telescope_Name = 'GlobalSignal_ConstantError'

	ObsDuration = '1000hr'

	# Total noise sensitivity as computed from 21cmSense. This total noise combines both the Poisson component of the observationally measured PS 
	# (the high res 21cm FAST simulation) and the instrument thermal noise profile.
	# Note: These will be interpolated in the computation of the likelihood statistic
	multi_z_Error_k = []
	multi_z_Error_PS = []

	if USE_GLOBAL_SIGNAL is True:

		ModelName = 'FaintGalaxies'
#		ModelName = 'BrightGalaxies'
		
		ErrorFileName = 'TotalError_FaintGalaxies_%s_%s'%(Telescope_Name,ObsDuration)
#		ErrorFileName = 'TotalError_BrightGalaxies_%s_%s'%(Telescope_Name,ObsDuration)

		# For the global signal, one must define the frequency region over which the global signal will be fit. The user must ensure that the frequency coverage
		# matches the sampled range of redshifts in the MCMC
		# The minimum and maximum frequencies (in MHz) for which the global signal will be fit
		FrequencyMin = 40.
		FrequencyMax = 200.

		# The user can either choose a fixed error on the global signal, or read in from a text file
		if GLOBAL_SIGNAL_FIXED_ERROR is True:			
			# Fixed error (in mK) on the 21cm global signal
			ErrorOnGlobal = 10.

			# The binning (frequency bandwidth) of the global signal. If the user chooses to read from file, this value is not used
			Bandwidth = 4.

			# Fill the data with empty lists. It will not be used anyway
			Error_k_values = []
			Error_PS_values = []

		else:				
			# Read in the global signal noise from file.
			### NOTE ###
			# By reading in from file, this adopts the redshift sampling used in the file. Therefore, make sure that the text file containing the error
			# includes the requisite redshift sampling that you want (sampling of the full spectrum using spline interpolation, not redshift sampling of the 21cmFAST output)
			### NOTE ###
			# It will use the values set for FrequencyMin and FrequencyMax, so set these values accordingly

			Error_k_values = numpy.loadtxt('NoiseData/%s/GlobalSignal/%s.txt'%(ModelName,ErrorFileName), usecols=(0,))
			Error_PS_values = numpy.loadtxt('NoiseData/%s/GlobalSignal/%s.txt'%(ModelName,ErrorFileName), usecols=(1,))

		multi_z_Error_k.append(Error_k_values)
		multi_z_Error_PS.append(Error_PS_values)

	else:
		
		if IncludeLightCone is True:

			# *** NOTE *** Again, the names of these files should be placed in a text-file in ** increasing ** redshift order
			MockObsFileName = 'LightCone21cmPS_Error_FaintGalaxies_%s_%s_600Mpc_400'%(Telescope_Name,ObsDuration)
		
			filename = open('NoiseData/%s.txt'%(MockObsFileName), 'r') 
			LightConeErrors = [line.rstrip('\n') for line in filename]

			# Use LightConeSnapShots here to ensure it crashes if the number of error files is less than the number or observations
			for i in range(len(Redshift)):
				Error_k_values = numpy.loadtxt('NoiseData/%s'%(LightConeErrors[i]), usecols=(0,))
				Error_PS_values = numpy.loadtxt('NoiseData/%s'%(LightConeErrors[i]), usecols=(1,))

				multi_z_Error_k.append(Error_k_values)
				multi_z_Error_PS.append(Error_PS_values)

		else:

#			NoiseFileName = 'TotalError_%s_PS_600Mpc'%(Telescope_Name)
			NoiseFileName = 'TotalError_%s_PS_200Mpc'%(Telescope_Name)

#			ModelName = 'FaintGalaxies'
			ModelName = 'EOS_FaintGalaxies'
			BoxType = 'Co-Eval'

			for i in range(len(Redshift)):
#				Error_k_values = numpy.loadtxt('NoiseData/FaintGalaxies/%s/%s_z%s_%s.txt'%(BoxType,NoiseFileName,Redshift[i],ObsDuration), usecols=(0,))
#				Error_PS_values = numpy.loadtxt('NoiseData/FaintGalaxies/%s/%s_z%s_%s.txt'%(BoxType,NoiseFileName,Redshift[i],ObsDuration), usecols=(1,))
				Error_k_values = numpy.loadtxt('NoiseData/%s/%s/%s_z%s_%s_%s.txt'%(ModelName,BoxType,NoiseFileName,Redshift[i],ModelName,ObsDuration), usecols=(0,))
				Error_PS_values = numpy.loadtxt('NoiseData/%s/%s/%s_z%s_%s_%s.txt'%(ModelName,BoxType,NoiseFileName,Redshift[i],ModelName,ObsDuration), usecols=(1,))

				multi_z_Error_k.append(Error_k_values)
				multi_z_Error_PS.append(Error_PS_values)
		
	multi_z_Error_k = numpy.array(multi_z_Error_k)
	multi_z_Error_PS = numpy.array(multi_z_Error_PS)	

	# k-space cut to be made to the likelihood fitting corresponding to the removal of foregrounds
	# NOTE: Be careful that the choices of "foreground_cut" and "shot_noise_cut" are contained within the ranges of 1) The mock observation PS 2) Error PS and 
	# 3) the boxes to be sampled by the MCMC. For all, a spline interpolation is performed, therefore it is imperative that you do not define these ranges outside
	# of the bounds of your data.
	foreground_cut = 0.15
	shot_noise_cut = 1.0
	# Number of spline points to be used for the likelihood computation
	NSplinePoints = 8

	# Whether to sample the 21cm PS (k-space) in linear or log space (more important as one lowers the "foreground_cut" limit)
	Log_k_sampling = False


	# Including a modelling uncertainty. This accounts for differences between semi-numerical and N-body/Hydro/RT prescriptions for EoR simulations. Setting to 0.0
	# will remove the uncertainty. 0.10 corresponds to 10 per cent modelling uncertainty. If set > 0.0, this error adds in quadrature with the telescope sensitivity 
	# for each k-mode.
	ModUncert = 0.2

	# To be added to the string name for the output MCMC data
	if IncludeLightCone is True:
		multiz_flag = 'lightcone'
	else:
		if len(Redshift) == 1:
			multiz_flag = 'Co-eval_singlez'
		else:
			multiz_flag = 'Co-eval_multiz'



	############### Creating the MCMC parameter ranges for each of the EoR and X-ray heating parameters to be sampled ##################



	# A dictionary to define whether or not to vary a parameter for the MCMC
	param_legend = dict()

	# Create list of the parameter names and the lower and upper limits for each variable.
	# This is a little clumsy, but makes things a little easier. Probably could be made more sophisticated if I thought about it
	param_string_names = []
	param_lower_limits = []
	param_upper_limits = []

	# Setting up a dictionary of the available parameters to be sampled. Set "True" to allow this parameter to be varied, set "False" if it is to be held fixed

	# Alpha (power-law). Note: Current version does not support a mass-dependent power law ionising efficiency to be used with the spin temperature fluctuations
	param_legend['ALPHA'] = False	

	# Set a fiducial value for the power law index alpha, and its lower and upper bounds. Not all will be used, depends on what options are set.
	Fiducial_Alpha = 0.0
	LowerBound_Alpha = -2.0
	UpperBound_Alpha = 2.0

	param_string_names.append('ALPHA')
	param_lower_limits.append(LowerBound_Alpha)
	param_upper_limits.append(UpperBound_Alpha)

	# Ionising efficiency, Zeta
	param_legend['ZETA'] = True

	# Set a fiducial value for Zeta, and its lower and upper bounds. Not all will be used, depends on what options are set.
	Fiducial_Zeta = 30.0
	LowerBound_Zeta = 10.0
	UpperBound_Zeta = 250.0

	param_string_names.append('ZETA')
	param_lower_limits.append(LowerBound_Zeta)
	param_upper_limits.append(UpperBound_Zeta)

	# Mean ionising photon horizon, R_mfp
	param_legend['MFP'] = True

	if USE_INHOMO_RECO is True:
		# If inhomogeneous recombinations are set, then the mean free path is not varied, but fixed to 50 Mpc
		param_legend['MFP'] = False

	# Set a fiducial value for R_mfp, and its lower and upper bounds. Not all will be used, depends on what options are set.
	if param_legend['MFP'] is True:
		Fiducial_MFP = 15.0
		LowerBound_MFP = 5.0
		UpperBound_MFP = 25.0

		param_string_names.append('MFP')
		param_lower_limits.append(LowerBound_MFP)
		param_upper_limits.append(UpperBound_MFP)	

	else:
		Fiducial_MFP = 50.0

	# Minimum halo mass hosting sources contributing to X-ray heating or ionisation. (Note: Currently version assumes ION_TVIR_MIN = X-RAY_TVIR_MIN a future release
	# will remove this assumption)
	param_legend['TVIR_MIN'] = True

	# Set a fiducial value for ION_TVIR_MIN (and X-RAY_TVIR_MIN), and its lower and upper bounds. Not all will be used, depends on what options are set.
	# Defined as log10(ION_TVIR_MIN). E.g. 4.69897 = log10(5x10^4)
	Fiducial_TVIR = 4.69897
	# Under the assumption ION_TVIR_MIN == X-RAY_TVIR_MIN it is important that this range is set to the lower and upper limits for X-ray heating
	# Allowed to be written like this as if you are not performing the spin temperature fluctuations calculation, easier to set here.
	LowerBound_TVIR = X_RAY_TVIR_LB
	UpperBound_TVIR = X_RAY_TVIR_UB

	param_string_names.append('TVIR_MIN')
	param_lower_limits.append(LowerBound_TVIR)
	param_upper_limits.append(UpperBound_TVIR)	

	# Soft band X-ray luminosity, L_X. Used for determining the number of X-ray photons produced per stellar baryon
	param_legend['L_X'] = True

	# Set a fiducial value for L_X, and its lower and upper bounds. Not all will be used, depends on what options are set.
	# Defined as log10(L_X). E.g. 40 = log10(10^40)
	Fiducial_LX = 40.0
	LowerBound_LX = 38.0
	UpperBound_LX = 42.0

	param_string_names.append('L_X')
	param_lower_limits.append(LowerBound_LX)
	param_upper_limits.append(UpperBound_LX)	

	# Minimum frequency X-ray photon contributing to IGM heating and ionization
	param_legend['NU_X_THRESH'] = True

	# Set a fiducial value for NU_X_THRESH, and its lower and upper bounds. Not all will be used, depends on what options are set.
	# Defined in eV. E.g 500 = 0.5keV
	Fiducial_NU_X_THRESH = 500.0
	LowerBound_NU_X_THRESH = 100.0
	UpperBound_NU_X_THRESH = 1500.0

	param_string_names.append('NU_X_THRESH')
	param_lower_limits.append(LowerBound_NU_X_THRESH)
	param_upper_limits.append(UpperBound_NU_X_THRESH)	

	# X-Ray spectral index at frequencies higher than NU_X_THRESH
	param_legend['X_RAY_SPEC_INDEX'] = True

	# Set a fiducial value for X_RAY_SPEC_INDEX, and its lower and upper bounds. Not all will be used, depends on what options are set.
	Fiducial_X_RAY_SPEC_INDEX = 1.0
	LowerBound_X_RAY_SPEC_INDEX = -1.0
	UpperBound_X_RAY_SPEC_INDEX = 3.0

	param_string_names.append('X_RAY_SPEC_INDEX')
	param_lower_limits.append(LowerBound_X_RAY_SPEC_INDEX)
	param_upper_limits.append(UpperBound_X_RAY_SPEC_INDEX)	

	############### Condition: If Include_Ts_fluc = False, then X-ray heating parameters cannot be varied.  ###############
	############### Overwrite any param_legend values for the X-ray parameters to False 					###############

	if Include_Ts_fluc is False:
			param_legend['L_X'] = False
			param_legend['NU_X_THRESH'] = False
			param_legend['X_RAY_SPEC_INDEX'] = False		



	############### Now set up the cosmological values #######################



	# Set up ranges and values for any of the cosmological values to be varied

	# Sigma 8 (normalisation of the 21cm PS)
	Fiducial_Sigma8 = 0.820000
	
	LowerBound_SIGMA_8 = 0.7
	UpperBound_SIGMA_8 = 0.95

	if 'SIGMA_8' in CosmologyToVary:
		param_legend['SIGMA_8'] = True

		param_string_names.append('SIGMA_8')
		param_lower_limits.append(LowerBound_SIGMA_8)
		param_upper_limits.append(UpperBound_SIGMA_8)
	else:
		param_legend['SIGMA_8'] = False		

	# little h (value of H_0)
	Fiducial_littleh = 0.680000

	LowerBound_littleh = 0.65
	UpperBound_littleh = 0.71

	if 'littleh' in CosmologyToVary:
		param_legend['littleh'] = True
	
		param_string_names.append('littleh')
		param_lower_limits.append(LowerBound_littleh)
		param_upper_limits.append(UpperBound_littleh)
	else:
		param_legend['littleh'] = False

	# Omega matter, fraction of mass, and correspondingly, the dark energy fraction
	# Note here that for Omega M, we implicitly assume Omega L = 1 - Omega M
	Fiducial_Omega_M = 0.310000

	LowerBound_OMEGA_M = 0.25
	UpperBound_OMEGA_M = 0.35

	if 'OMEGA_M' in CosmologyToVary:
		param_legend['OMEGA_M'] = True
	
		param_string_names.append('OMEGA_M')
		param_lower_limits.append(LowerBound_OMEGA_M)
		param_upper_limits.append(UpperBound_OMEGA_M)
	else:
		param_legend['OMEGA_M'] = False
	# Omega baryon, the baryon component
	Fiducial_Omega_b = 0.048

	LowerBound_OMEGA_b = 0.0470
	UpperBound_OMEGA_b = 0.0490

	if 'OMEGA_b' in CosmologyToVary:
		param_legend['OMEGA_b'] = True

		param_string_names.append('OMEGA_b')
		param_lower_limits.append(LowerBound_OMEGA_b)
		param_upper_limits.append(UpperBound_OMEGA_b)
	else:
		param_legend['OMEGA_b'] = False

	Fiducial_ns = 0.97000

	LowerBound_NS = 0.9
	UpperBound_NS = 1.0

	if 'NS' in CosmologyToVary:
		param_legend['NS'] = True

		param_string_names.append('NS')
		param_lower_limits.append(LowerBound_NS)
		param_upper_limits.append(UpperBound_NS)
	else:
		param_legend['NS'] = False

	# Some other parameters (which in future can be varied)

	# Star-formation time-scale as a fraction of the Hubble time
	Fiducial_t_STAR = 0.5
	# The fraction of baryons converted to stars
	Fiducial_F_STAR = 0.05



	############### Collate the data and options set above for the MCMC #######################



	# This dictionary stored all the default (fiducial) parameters set by the user. Will use these values for any X-ray or EoR parameter that is held fixed
	Fiducial_Params = dict()

	# All these parameters are passed to a text-file for which the 21CMMC driver can read in and use these parameters to perform any computation.
	Fiducial_Params['ALPHA'] = Fiducial_Alpha
	Fiducial_Params['ZETA'] = Fiducial_Zeta
	Fiducial_Params['MFP'] = Fiducial_MFP
	Fiducial_Params['TVIR_MIN'] = Fiducial_TVIR
	Fiducial_Params['L_X'] = Fiducial_LX
	Fiducial_Params['NU_X_THRESH'] = Fiducial_NU_X_THRESH
	Fiducial_Params['X_RAY_SPEC_INDEX'] = Fiducial_X_RAY_SPEC_INDEX

	Fiducial_Params['NU_X_BAND_MAX'] = NU_X_BAND_MAX
	Fiducial_Params['NU_X_MAX'] = NU_X_MAX
	Fiducial_Params['X_RAY_TVIR_MIN'] = Fiducial_Params['TVIR_MIN']
	Fiducial_Params['X_RAY_TVIR_LB'] = X_RAY_TVIR_LB
	Fiducial_Params['X_RAY_TVIR_UB'] = X_RAY_TVIR_UB

	Fiducial_Params['t_STAR'] = Fiducial_t_STAR
	Fiducial_Params['F_STAR'] = Fiducial_F_STAR

	Fiducial_Params['N_RSD_SUBCELLS'] = N_RSD_SUBCELLS
	Fiducial_Params['LOS_direction'] = LOS_direction	

	Fiducial_Params['SIGMA_8'] = Fiducial_Sigma8	
	Fiducial_Params['littleh'] = Fiducial_littleh
	Fiducial_Params['OMEGA_M'] = Fiducial_Omega_M
	Fiducial_Params['OMEGA_b'] = Fiducial_Omega_b
	Fiducial_Params['NS'] = Fiducial_ns

	### Some Global signal parameters ###
	if USE_GLOBAL_SIGNAL is True:
		Fiducial_Params['MIN_FREQ'] = FrequencyMin
		Fiducial_Params['MAX_FREQ'] = FrequencyMax

		if GLOBAL_SIGNAL_FIXED_ERROR is True:			
			Fiducial_Params['CONST_ERROR'] = ErrorOnGlobal
			Fiducial_Params['BANDWIDTH'] = Bandwidth

	numpy.array(param_lower_limits)
	numpy.array(param_upper_limits)

	# Determine which variables are to be varied, and populate the array to be passed to the MCMC algorithm to be sampled
	params = []
	for i in range(len(param_string_names)):
		if param_legend[param_string_names[i]] is True:
			params_new = [param_lower_limits[i],param_upper_limits[i]]

			params.append(params_new)


	params = numpy.array(params)

	FlagOptions = dict()
	FlagOptions['GENERATE_NEW_ICS'] = GenerateNewICs
	FlagOptions['INCLUDE_RSDS'] = INCLUDE_SUBCELL_RSDS
	FlagOptions['USE_INHOMO_RECO'] = USE_INHOMO_RECO
	FlagOptions['USE_FCOLL_TABLE'] = UseFcollTable
	FlagOptions['CALC_TS_FLUC'] = Include_Ts_fluc
	FlagOptions['KEEP_GLOBAL_DATA'] = USE_GLOBAL_SIGNAL
	FlagOptions['USE_IONISATION_FCOLL_TABLE'] = USE_IONISATION_FCOLL_TABLE
	FlagOptions['INCLUDE_POWERLAW'] = IncludeAlpha
	FlagOptions['USE_GS_FIXED_ERROR'] = GLOBAL_SIGNAL_FIXED_ERROR
	
	FlagOptions['KEEP_ALL_DATA'] = KEEP_ALL_DATA

	FlagOptions['LOG_LINEAR_K_SAMPLING'] = Log_k_sampling

	Create_Output_Directory = 'KEEP_MCMC_DATA_%s'%(time.strftime("%a_%d_%b_%Y_%Hh_%Mm_%Ss"))

	if KEEP_ALL_DATA is True:
		command = "mkdir %s"%(Create_Output_Directory)
		os.system(command)
		command = "mkdir %s/AveData"%(Create_Output_Directory)
		os.system(command)
		command = "mkdir %s/StatisticalData"%(Create_Output_Directory)
		os.system(command)
		command = "mkdir %s/TauData"%(Create_Output_Directory)
		os.system(command)
		command = "mkdir %s/WalkerData"%(Create_Output_Directory)
		os.system(command)

	FlagOptions['KEEP_ALL_DATA_FILENAME'] = Create_Output_Directory


	############ Rudimentary error checking ###############



	# Almost certainly a more sophisticated method can be used in Python, but I can't be bothered looking into it

	ErrorMessage = False
	ErrorString = []

	if GenerateNewICs is False and len(CosmologyToVary) > 0.0:
		ErrorString.append("ERROR: Cosmological parameters have been set to vary but the initial conditions (density field) is not set to be generated on the fly")
		ErrorString.append("(If the user wants to vary cosmology, the GenerateNewICs flag to true).")
		ErrorMessage = True

	if USE_GLOBAL_SIGNAL is True and USE_IONISATION_FCOLL_TABLE is True:
		ErrorString.append("ERROR: Cannot use the interpolation table for the find_HII_bubbles part when generating the global signal")
		ErrorString.append("(In principle you should be able to, but, instead one can decrease the redshift sampling of the global signal instead)")
		ErrorString.append("(i.e. set co-eval, and put in every second/third redshift that Ts should sample.")
		ErrorMessage = True

	if USE_GLOBAL_SIGNAL is True and IncludeLightCone is False and (len(Redshift) + len(Redshifts_For_Prior)) < 20:
		ErrorString.append("ERROR: Too few redshift bins to accurately determine the global signal (i.e. spline interpolation)")
		ErrorString.append("Either set IncludeLightCone = True or add more redshifts to either of Redshift or Redshifts_For_Prior lists")
		ErrorString.append("(Or if you are bold, reduce this arbitrary number from 20 to something smaller!)")
		ErrorString.append("(Ensure the chosen redshifts include the endpoints of the global signal redshift (frequency) coverage)")
		ErrorMessage = True		

	if GenerateNewICs is True and (UseFcollTable or USE_IONISATION_FCOLL_TABLE is True):
		ErrorString.append("ERROR: Cannot use interpolation tables when generating new initial conditions on the fly.")
		ErrorMessage = True
	
	if USE_INHOMO_RECO is True and USE_IONISATION_FCOLL_TABLE is True:
		ErrorString.append("ERROR: Cannot use the f_coll ionisation table (from find_HII_bubbles) in conjuction with inhomogeneous recombinations")
		ErrorString.append("Require the full boxes for keeping track of the recombinations")
		ErrorMessage = True

	if USE_INHOMO_RECO is True and Include_Ts_fluc is False:
		ErrorString.append("ERROR: Inhomogeneous recombinations can only be used in combination with the full computation of the IGM spin temperature")
		ErrorString.append("This differs from 21cmFAST, but arises as no intermediate information is stored")
		ErrorMessage = True

	if Include_Ts_fluc is True and param_legend['ALPHA'] is True:
		ErrorString.append("ERROR: Current version does not support a mass-dependent ionising efficiency to be used with the full spin temperature fluctuations.")
		ErrorMessage = True

	if (len(Redshift) + len(Redshifts_For_Prior)) < 3:

		if PriorLegend['PlanckPrior'] is True and IncludeLightCone is False:
			ErrorString.append("ERROR: planck prior cannot be used as insufficient redshift boxes have been chosen.")
			ErrorString.append("A minimum of three independent redshifts are required to estimate tau")
			ErrorString.append("Consider adding additional redshifts either to:")
			ErrorString.append("1) the likelihood calculation (i.e. the Redshift list)")
			ErrorString.append("2) the redshifts for prior only list (i.e. Redshifts_For_Prior)")
			ErrorMessage = True

		if PriorLegend['McGreerPrior'] is True and IncludeLightCone is False:
			if '%s'%(McGreer_Redshift) not in Redshift:
				ErrorString.append("ERROR: The McGreer et al. prior cannot be used as insufficient redshift boxes have been chosen for the interpolation/extrapolation.")
				ErrorString.append("A minimum of three independent redshifts are required (unless the sampler samples the exact redshift of the prior)")
				ErrorString.append("Consider adding additional redshifts either to:")
				ErrorString.append("1) the likelihood calculation (i.e. the Redshift list)")
				ErrorString.append("2) the redshifts for prior only list (i.e. Redshifts_For_Prior)")
				ErrorString.append("3) the exact redshift of the McGreer et al. prior")
				ErrorMessage = True

		if PriorLegend['GreigPrior'] is True and IncludeLightCone is False:		
			if '%s'%(QSO_Redshift) not in Redshift:
				ErrorString.append("ERROR: The Greig et al. QSO Damping Wing prior cannot be used as insufficient redshift boxes have been chosen for the interpolation/extrapolation.")
				ErrorString.append("A minimum of three independent redshifts are required (unless the sampler samples the exact redshift of the prior)")
				ErrorString.append("Consider adding additional redshifts either to:")
				ErrorString.append("1) the likelihood calculation (i.e. the Redshift list)")
				ErrorString.append("2) the redshifts for prior only list (i.e. Redshifts_For_Prior)")
				ErrorString.append("3) the exact redshift of the Greig et al. prior")
				ErrorMessage = True
	


	############## Setup MCMC #################



	# Start the MCMC sampling. The only things that need to be modified are "walkersRatio", "burninIterations", "sampleIterations" and "threadCount". "filethin" can be used to thin the size of the stored
	# MCMC walker positions (I included this for a different application). Should be left at 1.
	# - "walkersRatio" must be set so that "walkersRatio"*"len(params)" is divisible by two. This is by construction of the EMCEE sampling.
	# - "threadCount" defines the number of threads to be used. The maximum number of threads that can be used is walkersRatio*len(params)/2. Therefore, to maximise resources
	# of your machine, increase walkersRatio.
	# - The stop criteria is set to max iterations, therefore the sampler stops only at the end of iterations. Make sure to choose a sufficiently high number for convergence. 
	# Using reuseBurnin=True will read in from file the last locations of the walkers to "continue" the sampling if the sampler hasn't sufficently converged.
	# From experience I find it preferable to increase "walkersRatio" over "sampleIterations" as "walkersRatio" allows the parameter space to be sampled at more locations (more walkers)

	# Redshift takes the lowest (latest) redshift and uses that in the initial condition generator (can remove this from the IC generator)

	# Typically for 3 parameters, I would perform the following:
	# walkersRatio = 16, burninIterations = 250, sampleIterations = 3000, threadCount = 24 (for a 24 core machine. Note 3*16/2 = 24)

	chain = LikelihoodComputationChain()
	
	Likelihoodmodel21cmFast = Likelihood21cmFast_multiz(multi_z_mockobs_k,multi_z_mockobs_PS,multi_z_Error_k,multi_z_Error_PS,
			Redshift,Redshifts_For_Prior,param_legend,Fiducial_Params,FlagOptions,param_string_names,NSplinePoints,TsCalc_z,foreground_cut,shot_noise_cut,IncludeLightCone,
			ModUncert,PriorLegend,NFVals_QSODamping,PDFVals_QSODamping)	

	chain.addLikelihoodModule(Likelihoodmodel21cmFast)

	chain.setup()

	File_String = 'ReionModel_21cmFast_%s_%s'%(Telescope_Name,multiz_flag)
	
	sampler = CosmoHammerSampler(
                    params = params,
                    likelihoodComputationChain=chain,
                    filePrefix="%s"%(File_String),
                    walkersRatio=2,
                    FiducialParams=Fiducial_Params,
                    param_legend=param_legend,
                    LowerBound_XRAY=X_RAY_TVIR_LB,
                    UpperBound_XRAY=X_RAY_TVIR_UB,
                    SpinTz=TsCalc_z,
                    burninIterations=0,
                    sampleIterations=1,
                    filethin = 1,
                    threadCount=4,
	                reuseBurnin=False
	           	)

	"""
	sampler = CosmoHammerSampler(
                    params = params,
                    likelihoodComputationChain=chain,
                    filePrefix="%s"%(File_String),
                    walkersRatio=72,
                    CreateFFTData=StoreFFTData,
                    CreateFCollData=CreateFcollTable,
                    LowerBound_XRAY=X_RAY_TVIR_LB,
                    UpperBound_XRAY=X_RAY_TVIR_UB,
                    SpinTz=TsCalc_z,
                    ThreadsForTable=FcollTableThreads,
                    burninIterations=50,
                    sampleIterations=350,
                    filethin = 1,
                    threadCount=216,
	                reuseBurnin=False
	           	)
	
	"""



	if ErrorMessage is True:

		for i in range(len(ErrorString)):
			print ErrorString[i]
	
	else:
		print 'Start sampling'
		sampler.startSampling()            
			