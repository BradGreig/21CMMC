import numpy
import math
from scipy import interpolate
from decimal import *

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')

from CosmoHammer_21CMMC.sampler.CosmoHammerSampler import CosmoHammerSampler
from CosmoHammer_21CMMC.likelihood.chain.LikelihoodComputationChain import LikelihoodComputationChain
from CosmoHammer_21CMMC.likelihood.module.Likelihood21cmFast import Likelihood21cmFast_multiz

if __name__ == '__main__':

	# Reionisation redshifts for the multi-z 21cm Fast "observations"
	# Need to make sure that these boxes exist in the "Boxes" folder. If not, please generate new boxes
#	Redshift = [8., 9., 10.]
	Redshift = [9.]

	# Set for the desired telescope ['SKA_halveddipoles_compact', 'HERA331'] corresponding to the file structure in "NoiseData"
	Telescope_Name = 'HERA331'
	ErrorFileTelescope_Name = '%s'%(Telescope_Name)

	# Fiducial model parameters. Only used if these are saved in the file names for the error files
	Zeta_val = 30.0
	MFP_val = 15.0
	TVir_val = 30000.0

	# Read in the length of the PS files for numpy arrays. The length of the files should be consistent
	# Note: This should be removed from final version.
	k_values = numpy.loadtxt('NoiseData/MockObs_%s_PS_250Mpc_z%s_%s_%s_%s.txt'%(Telescope_Name,Redshift[0],Zeta_val,MFP_val,TVir_val), usecols=(0,))
	
	multi_z_mockobs_k = numpy.zeros((len(Redshift),len(k_values)))
	multi_z_mockobs_PS = numpy.zeros((len(Redshift),len(k_values)))
	multi_z_PS_Error = numpy.zeros((len(Redshift),len(k_values)))

	# Read in the mock 21cm PS observation. Read in both k and the dimensionless PS
	# Note: The mock observations need to be binned "exactly" the same as the resolution of the sampled boxes by 21cmFAST. If not, the chi^2 statistic is still going to 
	# be computed, but it will be wrong. I need to modify 21CMMC to interpolate the 21cm PS of the "observation" in the computation step.
	for i in range(len(Redshift)):
		multi_z_mockobs_k[i] = numpy.loadtxt('NoiseData/MockObs_%s_PS_250Mpc_z%s_%s_%s_%s.txt'%(Telescope_Name,Redshift[i],Zeta_val,MFP_val,TVir_val), usecols=(0,))
		multi_z_mockobs_PS[i] = numpy.loadtxt('NoiseData/MockObs_%s_PS_250Mpc_z%s_%s_%s_%s.txt'%(Telescope_Name,Redshift[i],Zeta_val,MFP_val,TVir_val), usecols=(1,))

	# k-space cut to be made to the likelihood fitting corresponding to the removal of foregrounds
	foreground_cut = 0.15
	shot_noise_cut = 1.0

	# Including a modelling uncertainty. This accounts for differences between semi-numerical and N-body/Hydro/RT prescriptions for EoR simulations. Setting to 0.0
	# will remove the uncertainty. 0.10 corresponds to 10 per cent modelling uncertainty
	ModUncert = 0.10

	# Total noise sensitivity as computed from 21cmSense. This total noise combines both the Poisson component of the observationally measured PS 
	# (the high res 21cm FAST simulation) and the instrument thermal noise profile.
	for i in range(len(Redshift)):
		multi_z_PS_Error[i] = numpy.loadtxt('NoiseData/TotalError_%s_PS_250Mpc_z%s_%s_%s_%s.txt'%(ErrorFileTelescope_Name,Redshift[i],Zeta_val,MFP_val,TVir_val), usecols=(1,))

	# Bad value flag returns a value far away from the highest likelihood value to ensure the MCMC never enters into this region. This flag will be typically 
	# raised for a measured PS from a completely ionised simulation box. Hence, is just the chi^{2} for the PS itself (i.e. PS zero) multiplied by 100 to really
	# ensure this parameter combination can never be chosen
	bad_value_flag = -0.5*(100.*sum(numpy.square((multi_z_mockobs_PS[0])/multi_z_PS_Error[0])))
	
	if len(Redshift) == 1:
		multiz_flag = 'singlez'
	else:
		multiz_flag = 'multiz'

	# Creating the list of parameters to be handed to the CosmoHammer (EMCEE) sampler
	params = [[]]

	""" This needs to be modified to use the dictionary class to make it more intelligent as to which parameters are free or fixed. 
	At this point it is hard coded into the likelihood calculation itself """

#   Set the parameter range for Zeta (ionising efficiency). CosmoHammer takes 4 values, 2nd and 3rd are min and max allowed for the parameter, first and last
#	are not really important. The first is set to the middle of the range and the last is like a one-sigma val. I just set this to be about 1-2 sigma of my parameter range.
#	It is an MCMC, so will sample the entire parameter space. It is most important for the initial locations of the "walkers", but I switched the initial position generator
	params = [[50.,5.,100.,20.]]

#	Parameter range for mean free path
	params_new = [10.,2.,20.,3.]
	params.append(params_new)

#   Parameter range for Virial Temperature (note to sample the parameter space better, we deal with the log10 of the virial temperature)
	params_new = [numpy.log10(80000.),numpy.log10(10000.),numpy.log10(200000.),numpy.log10(50000.)]

	params.append(params_new)        

	params = numpy.array(params)

	chain = LikelihoodComputationChain()

	Likelihoodmodel21cmFast = Likelihood21cmFast_multiz(multi_z_mockobs_k,multi_z_mockobs_PS,multi_z_PS_Error,Redshift,foreground_cut,shot_noise_cut,ModUncert)
	
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
                    walkersRatio=4,
                    Redshift = Redshift[0],
                    burninIterations=2,
                    sampleIterations=2,
                    filethin = 1,
                    threadCount=6,
	                reuseBurnin=False
	           	)

	sampler.startSampling()            
	    
