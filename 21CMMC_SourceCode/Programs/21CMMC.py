import numpy
import math
from scipy import interpolate
from decimal import *
import string
import pickle

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')
SIXPLACES = Decimal(10) ** -6       # same as Decimal('0.000001')

from CosmoHammer_21CMMC.sampler.CosmoHammerSampler import CosmoHammerSampler
from CosmoHammer_21CMMC.likelihood.chain.LikelihoodComputationChain import LikelihoodComputationChain
from CosmoHammer_21CMMC.likelihood.module.Likelihood21cmFast import Likelihood21cmFast_multiz

from CosmoHammer_21CMMC.likelihood.module.Likelihood21cmFast import McGreer_Redshift
from CosmoHammer_21CMMC.likelihood.module.Likelihood21cmFast import QSO_Redshift

if __name__ == '__main__':

	# Setting IncludeAlpha = True, results in the sampling of a mass dependent ionising efficiency.
	# This mass dependent ionising efficiency is defined as (Mass / M_vir)**(Alpha)
	# - M_vir is the mass corresponding to the Virial Temperature (T_Vir) evaluated at the respective redshift (T_vir is defined to be redshift independent)
	# - Alpha is a the power law. Alpha = 0 corresponds to a mass independent ionising efficiency, corresponding to the default 3 parameter reionisation model
	# IncludeAlpha = True allows for a 4 parameter model to be sampled.
	IncludeAlpha = False 

	# Reionisation redshifts for the multi-z 21cm Fast "observations"
	# Need to make sure that these boxes exist in the "Boxes" folder. If not, please generate new boxes
	# Shifted the range to lower-z, motivated by the lower Planck (2016) values of tau
	# Expanded the default range to allow for an approximate computation of tau, to be applied as a prior
	Redshift = [6., 7., 8., 9., 10.]
#	Redshift = [8., 9., 10.]
	
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
	IncludePlanck = True
	# 2) The McGreer et al. prior is a upper limit on the IGM neutral fraction at 5.9
	# Modelled as a flat, unity prior at x_HI <= 0.06, and a one sided Gaussian at x_HI > 0.06 ( Gaussian of mean 0.06 and one sigma of 0.05 )
	IncludeMcGreer = True
	# 3) The Greig et al. prior is computed directly from the PDF of the IGM neutral fraction for the Small HII reionisation simulation
	IncludeGreig = True
	
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


	# Set for the desired telescope ['SKA_halveddipoles_compact', 'HERA331'] corresponding to the file structure in "NoiseData"
	Telescope_Name = 'HERA331'
	ErrorFileTelescope_Name = '%s'%(Telescope_Name)

	# Fiducial model parameters. Only used if these are saved in the file names for the error files
	Zeta_val = 20.0
	MFP_val = 15.0
	TVir_val = 30000.0
	Alpha_val = 0.0
	
#	Zeta_val = 15.0
#	MFP_val = 15.0
#	TVir_val = 50000.0
#	Alpha_val = 0.4	

	multi_z_mockobs_k = []
	multi_z_mockobs_PS = []

	separator = "_"
	seq = []
	seq.append("%s"%(Zeta_val))
	seq.append("%s"%(MFP_val))
	seq.append("%s"%(TVir_val))
	seq.append("%s"%(Alpha_val))

	StringArgument = string.join(seq,separator)

	# Read in the mock 21cm PS observation. Read in both k and the dimensionless PS. These are needed for performing the chi^2 statistic for the likelihood
	# Note: To calculate the likelihood statistic a spline is performed for each of the mock PS, simulated PS and the Error PS
	for i in range(len(Redshift)):		
		mockobs_k_values = numpy.loadtxt('NoiseData/MockObs_PS_250Mpc_z%s_%s.txt'%(Redshift[i],StringArgument), usecols=(0,))
		mockobs_PS_values = numpy.loadtxt('NoiseData/MockObs_PS_250Mpc_z%s_%s.txt'%(Redshift[i],StringArgument), usecols=(1,))

		multi_z_mockobs_k.append(mockobs_k_values)
		multi_z_mockobs_PS.append(mockobs_PS_values)

	multi_z_mockobs_k = numpy.array(multi_z_mockobs_k)
	multi_z_mockobs_PS = numpy.array(multi_z_mockobs_PS)

	# k-space cut to be made to the likelihood fitting corresponding to the removal of foregrounds
	foreground_cut = 0.15
	shot_noise_cut = 1.0

	# Including a modelling uncertainty. This accounts for differences between semi-numerical and N-body/Hydro/RT prescriptions for EoR simulations. Setting to 0.0
	# will remove the uncertainty. 0.10 corresponds to 10 per cent modelling uncertainty
	ModUncert = 0.10

	multi_z_Error_k = []
	multi_z_Error_PS = []

	# Total noise sensitivity as computed from 21cmSense. This total noise combines both the Poisson component of the observationally measured PS 
	# (the high res 21cm FAST simulation) and the instrument thermal noise profile.
	# Note: These will be interpolated in the computation of the likelihood statistic
	for i in range(len(Redshift)):
		Error_k_values = numpy.loadtxt('NoiseData/TotalError_%s_PS_250Mpc_z%s_%s_1000hr.txt'%(ErrorFileTelescope_Name,Redshift[i],StringArgument), usecols=(0,))
		Error_PS_values = numpy.loadtxt('NoiseData/TotalError_%s_PS_250Mpc_z%s_%s_1000hr.txt'%(ErrorFileTelescope_Name,Redshift[i],StringArgument), usecols=(1,))

		multi_z_Error_k.append(Error_k_values)
		multi_z_Error_PS.append(Error_PS_values)

	multi_z_Error_k = numpy.array(multi_z_Error_k)
	multi_z_Error_PS = numpy.array(multi_z_Error_PS)

	# Bad value flag returns a value far away from the highest likelihood value to ensure the MCMC never enters into this region. This flag will be typically 
	# raised for a measured PS from a completely ionised simulation box. Hence, is just the chi^{2} for the PS itself (i.e. PS zero) multiplied by 100 to really
	# ensure this parameter combination can never be chosen
	bad_value_flag = -0.5*(100.*sum(numpy.square((multi_z_mockobs_PS[0])/multi_z_Error_PS[0])))
	
	if len(Redshift) == 1:
		multiz_flag = 'singlez'
	else:
		multiz_flag = 'multiz'

	# Creating the list of parameters to be handed to the CosmoHammer (EMCEE) sampler
	params = [[]]

	""" This needs to be modified to use the dictionary class to make it more intelligent as to which parameters are free or fixed. 
	At this point it is hard coded into the likelihood calculation itself """

#   Set the parameter range for Zeta (ionising efficiency). CosmoHammer takes 4 values, 2nd and 3rd are min and max allowed for the parameter, first and last
#	are not important. The first is set to the middle of the range and the last is like a one-sigma val. I just set this to be about 1-2 sigma of my parameter range.
#	It is an MCMC, so will sample the entire parameter space. It is most important for the initial locations of the "walkers", but I switched the initial position generator
	params = [[50.,5.,100.,20.]]

#	Parameter range for mean free path
	params_new = [10.,2.,20.,3.]
	params.append(params_new)

#   Parameter range for Virial Temperature (note to sample the parameter space better, we deal with the log10 of the virial temperature)
	params_new = [numpy.log10(80000.),numpy.log10(10000.),numpy.log10(200000.),numpy.log10(50000.)]
	params.append(params_new)        

#	If Alpha_val is zero, then only considering the three parameter model (constant ionising efficiency)	
	if IncludeAlpha is True:
		
#  		Parameter range for the power law alpha, for a mass dependendant ionising efficiency. Alpha == 0, corresponds to mass independent ionising efficiency
		params_new = [0.,-3.,3.,2.]
		params.append(params_new)        		

	params = numpy.array(params)

	chain = LikelihoodComputationChain()

	num_redshifts = len(Redshift)

	if num_redshifts < 3:

		if PriorLegend['PlanckPrior'] is True:
			print 'The planck prior will not be used as insufficient redshift boxes have been chosen'
			print 'A minimum of three independent redshifts are required to estimate tau'

		if PriorLegend['McGreerPrior'] is True:
			if McGreer_Redshift not in Redshift:
				print 'The McGreer et al. prior will not be used as insufficient redshift boxes have been chosen'
				print 'A minimum of three independent redshifts are required to interpolate/extrapolate the reionisation history'				

		if PriorLegend['GreigPrior'] is True:		
			if QSO_Redshift not in Redshift:
				print 'The Greig et al. QSO Damping Wing prior will not be used as insufficient redshift boxes have been chosen'
				print 'A minimum of three independent redshifts are required to interpolate/extrapolate the reionisation history'				

	Likelihoodmodel21cmFast = Likelihood21cmFast_multiz(multi_z_mockobs_k,multi_z_mockobs_PS,multi_z_Error_k,multi_z_Error_PS,
		Redshift,foreground_cut,shot_noise_cut,ModUncert,PriorLegend,NFVals_QSODamping,PDFVals_QSODamping)
	
	chain.addLikelihoodModule(Likelihoodmodel21cmFast)

	chain.setup()

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

	File_String = 'ReionModel_21cmFast_%s_%s'%(Telescope_Name,multiz_flag)
	
	sampler = CosmoHammerSampler(
                    params = params,
                    likelihoodComputationChain=chain,
                    filePrefix="%s"%(File_String),
                    walkersRatio=2,
                    Redshift = Redshift[0],
                    burninIterations=1,
                    sampleIterations=1,
                    filethin = 1,
                    threadCount=1,
	                reuseBurnin=False
	           	)

	sampler.startSampling()            
	
