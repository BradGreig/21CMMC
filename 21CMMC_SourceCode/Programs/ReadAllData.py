import os
import numpy as np
from scipy import interpolate
import scipy.stats as stats 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.patches as mpatches
import scipy.ndimage as ndimage
import pickle
import math
import scipy
import scipy.signal
import string
import sys



np.seterr(invalid='ignore', divide='ignore')
from decimal import *

FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')
SIXPLACES = Decimal(10) ** -6       # same as Decimal('0.0001')

import cosmolopy.constants as cc

#def main():
if __name__ == '__main__':

	# Script to read in all the data and separate out the MCMC walkers into initial conditions, burn-in and actual data
	# This script only reads all the data and all associated quantities together into memory.
	# Can modify this script to do whatever you choose

	#### NOTE ###
	# This only works on the actual MCMC data (i.e. not the burn-in data).

	# 1) Must manually provide the total number of redshifts used by 21CMMC.py (i.e. the sum total of Redshift and Redshifts_For_Prior if the co-eval option is set). 
	#	 If the light-cone option is set, must provide the total number of redshifts (number of redshifts in the AveData_* textfile)

	n_redshifts = 6

	# 2) Must manually provide the number of sample iterations and the walkers ratio chosen in 21CMMC.py

	walkers_ratio = 2
	sample_iterations = 4

	# 3) Provide the directory name for the stored data (21CMMC.py will create the directory with a timestamp)
	KeepDataString = 'KEEP_MCMC_DATA_Fri_06_Oct_2017_09h_22m_38s'

	Stored_AveData = 'AveData'
	Stored_TauData = 'TauData'
	Stored_PSData = 'StatisticalData'

	# Provide the remaining string to fill the file name (it'll read in all data, burnin and the actual data)
	# Switch the file name to include the directory if that is preferred
	FileNames = 'HERA331_Co-eval_multiz'

	# First the burn-in. It'll need the burn-in data to determine the accepted/rejected points of the actual MCMC data (i.e not the burn-in data)

	f = open('ReionModel_21cmFast_%sburnin.out'%(FileNames),'r')
	parameters_burnin = [[float(v) for v in line.rstrip('\n').split('\t')] for line in f.readlines()]
	f.close()

	f2 = open('ReionModel_21cmFast_%sburninprob.out'%(FileNames),'r')
	probs_burnin = [float(v) for v in f2.readlines()]
	f2.close()

	f3 = open('ReionModel_21cmFast_%sNFvalsburnin.out'%(FileNames),'r')
	NFvals_Burnin = [[float(v) for v in line.rstrip('\n').split('\t')] for line in f3.readlines()]
	f3.close()


	# Now the MCMC data

	f = open('ReionModel_21cmFast_%s.out'%(FileNames),'r')
	parameters = [[float(v) for v in line.rstrip('\n').split('\t')] for line in f.readlines()]
	f.close()

	f2 = open('ReionModel_21cmFast_%sprob.out'%(FileNames),'r')
	probs = [float(v) for v in f2.readlines()]
	f2.close()

	f3 = open('ReionModel_21cmFast_%sNFvals.out'%(FileNames),'r')
	NFvals = [[float(v) for v in line.rstrip('\n').split('\t')] for line in f3.readlines()]
	f3.close()


	# Number of total iterations for the burn-in and main MCMC. Useful for lengths of lists/arrays
	nsamples_burnin = len(probs_burnin)
	nsamples = len(probs)

	# Total number of sampled astrophysical/cosmological parameters
	n_parameters = len(parameters[0])

	# Total number of iterations
	n_iterations = len(parameters)

	# Total number of walkers
	n_walkers = walkers_ratio*n_parameters

	# Rudimentary check to determine if the inputs match the expected amount of data
	if walkers_ratio*n_parameters*sample_iterations != n_iterations:
		print 'ERROR: The inputs (walkers_ratio, sample_iterations) have not been added correctly. Make sure these numbers are correct'


	# Read in the file names of all the walker files. These aren't needed, but could come in use depending on the application
	path = '%s/WalkerData'%(KeepDataString)

	WalkerFileNames = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and 'Cosmology' not in f]
	WalkerFileNames_Cosmology = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and 'Cosmology' in f]



	# Set up a list to store True/False data which denotes whether the corresponding walker (sampled parameter set) was accepted/rejected
	LikelihoodCheck = []

	

	burnin_iterations = len(parameters_burnin)/n_walkers

	# A counter, which will be used to store the total number of accepted positions
	total_true = 0


	# This first loop checks whether the first positions in the main MCMC are new, relative to the final position of the burn-in data
	for ii in range(n_walkers):
		if parameters_burnin[ii+(burnin_iterations-1)*n_walkers] == parameters[ii]:
			# New proposed position was rejected
			LikelihoodCheck.append(False)
		else:
			# New proposed position was accepted
			LikelihoodCheck.append(True)

			total_true += 1

	
	# Now go through the entire main MCMC data and determine the positions accepted/rejected
	for i in range(sample_iterations-1):
		for ii in range(n_walkers):

			if parameters[ii+i*n_walkers] == parameters[ii+(i+1)*n_walkers]:
				LikelihoodCheck.append(False)
				# New proposed position was rejected
			else:
				LikelihoodCheck.append(True)
				# New proposed position was accepted

				total_true += 1



	separator = "_"
	seq = []
	# Add the random thread ID
	seq.append("%s"%(Decimal(repr(NFvals[0][n_redshifts])).quantize(SIXPLACES)))
	# Add the second ID
	seq.append("%s"%(Decimal(repr(NFvals[0][n_redshifts+1])).quantize(SIXPLACES)))

	StringArgument = string.join(seq,separator)





	# Read the k-values from the sampled PS. This remains constant for all redshifts and parameter samplings.
	AllData_kvals = np.loadtxt('%s/%s/TotalPSData_%s.txt'%(KeepDataString,Stored_PSData,StringArgument), usecols=(0,))



	# Create the arrays to contain all the data from the accepted points
	AllData_PSvals = np.zeros((total_true,n_redshifts,len(AllData_kvals)))

	AllData_zvals = np.zeros(n_redshifts)
	AllData_xHvals = np.zeros((total_true,n_redshifts))
	AllData_Tbvals = np.zeros((total_true,n_redshifts))

	AllData_TauVals = np.zeros(total_true)


	AllData_zvals = np.loadtxt('%s/%s/AveData_%s.txt'%(KeepDataString,Stored_AveData,StringArgument), usecols=(0,))


	index_true = 0

	for i in range(nsamples):

		if LikelihoodCheck[i] is True:
			# This astrophysical/cosmological parameter set has been accepted by the MCMC sampler. Read in and store the data in memory

			# Create the file_ending, which is specific to each individual MCMC walker
			seq = []
			# Add the random thread ID
			seq.append("%s"%(Decimal(repr(NFvals[0][n_redshifts])).quantize(SIXPLACES)))
			# Add the second ID
			seq.append("%s"%(Decimal(repr(NFvals[0][n_redshifts+1])).quantize(SIXPLACES)))

			StringArgument = string.join(seq,separator)


			# Read in all the power spectrum data (all k and PS)
			for j in range(n_redshifts):

				AllData_PSvals[index_true][j] = np.loadtxt('%s/%s/TotalPSData_%s.txt'%(KeepDataString,Stored_PSData,StringArgument), usecols=(j+1,))

			# Read in all the global data

			# Read in the IGM neutral fraction
			AllData_xHvals[index_true] = np.loadtxt('%s/%s/AveData_%s.txt'%(KeepDataString,Stored_AveData,StringArgument), usecols=(1,))

			# Read in the average brightness temperature contrast
			AllData_Tbvals[index_true] = np.loadtxt('%s/%s/AveData_%s.txt'%(KeepDataString,Stored_AveData,StringArgument), usecols=(2,))
	
			# Read in the electron scattering optical depth, \tau
			AllData_TauVals[index_true] = np.loadtxt('%s/%s/Tau_e_%s.txt'%(KeepDataString,Stored_TauData,StringArgument), usecols=(0,))

			index_true += 1
