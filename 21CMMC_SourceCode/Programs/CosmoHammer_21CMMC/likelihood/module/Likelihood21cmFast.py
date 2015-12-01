#!/usr/bin/env python

import os
import numpy as np
np.seterr(invalid='ignore', divide='ignore')
from decimal import *

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')

class Likelihood21cmFast_multiz(object):
    
    def __init__(self, k_values, PS_values, PS_Error, Redshift, Foreground_cut, Shot_Noise_cut, ModUncert):
        self.k_values = k_values
        self.PS_values = PS_values
        self.PS_Error = PS_Error
        self.Redshift = Redshift
        self.Foreground_cut = Foreground_cut
        self.Shot_Noise_cut = Shot_Noise_cut
        self.ModUncert = ModUncert

    def Likelihood(self,ctx):
        params = ctx.getParams()

        nf_vals = np.zeros(len(self.Redshift))

        # Generate a unique ID for each thread by sampling a randomly seeded distribution.
        # Given than file I/O needs to be unique to each thread, it is beneficial to provide a unique ID in the off chance that two different threads 
        # end up with the same walker position (same parameter set)
        np.random.seed()
        
        random_number = np.random.normal(size=1.0)

        # rounded values of the new parameters
        ZetaVal = Decimal(repr(params[0])).quantize(FOURPLACES)
        RmfpVal = Decimal(repr(params[1])).quantize(FOURPLACES)
        TvirVal = Decimal(repr(params[2])).quantize(FOURPLACES)

        Individual_ID = Decimal(repr(random_number[0])).quantize(FOURPLACES)

        total_sum = 0;
        for i in range(len(self.Redshift)):

            # Run 21CMMC (the 21cmFAST code) to perform the reionisation simulation and recover the neutral fraction and 21cm PS
            command = "./drive_21cmMC_streamlined %g %s %s %s %s 1"%(self.Redshift[i],Individual_ID,ZetaVal,RmfpVal,TvirVal) 
            os.system(command)

            # Read in the neutral fraction and 21cm PS for this parameter set and redshift
            k_values_estimate = np.loadtxt('delTps_estimate_%s_%s_%s_%s.txt'%(Individual_ID,ZetaVal,RmfpVal,TvirVal), usecols=(0,))
            PS_values_estimate = np.loadtxt('delTps_estimate_%s_%s_%s_%s.txt'%(Individual_ID,ZetaVal,RmfpVal,TvirVal), usecols=(1,))
            nf_value = np.loadtxt('NeutralFraction_%s_%s_%s_%s.txt'%(Individual_ID,ZetaVal,RmfpVal,TvirVal), usecols=(0,))

            nf_vals[i] = nf_value

            if nf_value == 0.:
                # neutral fraction of zero, means no 21cm PS. recover a large number to steer clear of this region of parameter space
                total_sum += 100.*sum(np.square((self.PS_values[i])/self.PS_Error[i]))
            else:
                # compute the chi-squared value compated to the mock observation for this redshift
                for ii in range(len(k_values_estimate)):
                    if k_values_estimate[ii] >= self.Foreground_cut and k_values_estimate[ii] < self.Shot_Noise_cut: 
                        total_sum += np.square((self.PS_values[i][ii] - PS_values_estimate[ii])/(np.sqrt(self.PS_Error[i][ii]**2. + (self.ModUncert*PS_values_estimate[ii])**2.))) 
       
            # remove the temporary files
            command = "rm delTps_estimate_%s_%s_%s_%s.txt"%(Individual_ID,ZetaVal,RmfpVal,TvirVal)
            os.system(command)
            command = "rm NeutralFraction_%s_%s_%s_%s.txt"%(Individual_ID,ZetaVal,RmfpVal,TvirVal)
            os.system(command)

        return -0.5*total_sum,nf_vals

    def computeLikelihood(self, ctx):

        return self.Likelihood(ctx)

    def setup(self):
        print "Likelihood Fitting for 21cm Fast (3 parameters, Zeta, Rmfp and Tvir)" 
