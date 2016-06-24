#!/usr/bin/env python

import os
import numpy as np
np.seterr(invalid='ignore', divide='ignore')
from decimal import *
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
import string

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')
SIXPLACES = Decimal(10) ** -6       # same as Decimal('0.000001')

McGreer_Redshift = 5.9

# The redshift of the QSO
QSO_Redshift = 7.0842

class Likelihood21cmFast_multiz(object):
    
    def __init__(self, k_values, PS_values, Error_k_values, PS_Error, Redshift, Foreground_cut, Shot_Noise_cut, ModUncert, PriorLegend, NFValsQSO, PDFValsQSO):
        self.k_values = k_values
        self.PS_values = PS_values
        self.Error_k_values = Error_k_values
        self.PS_Error = PS_Error
        self.Redshift = Redshift
        self.Foreground_cut = Foreground_cut
        self.Shot_Noise_cut = Shot_Noise_cut
        self.ModUncert = ModUncert
        self.PriorLegend = PriorLegend
        self.NFValsQSO = NFValsQSO
        self.PDFValsQSO = PDFValsQSO

    def Likelihood(self,ctx):
        params = ctx.getParams()

        nf_vals = np.zeros(len(self.Redshift))

        # Generate a unique ID for each thread by sampling a randomly seeded distribution.
        # Given than file I/O needs to be unique to each thread, it is beneficial to provide a unique ID in the off chance that two different threads 
        # end up with the same walker position (same parameter set)
        np.random.seed()
        
        random_number = np.random.normal(size=1.0)
    
        Individual_ID = Decimal(repr(random_number[0])).quantize(SIXPLACES)

        separator = " "
        separator_other = "_"
        seq = []
        # Add the random thread ID
        seq.append("%s"%(Individual_ID))
        # Add all rounded variables
        for i in range(len(params)):
            # rounded values of the new parameters
            seq.append("%s"%(Decimal(repr(params[i])).quantize(SIXPLACES)))

        StringArgument = string.join(seq,separator)
        StringArgument_other = string.join(seq,separator_other)

        total_sum = 0;
        for i in range(len(self.Redshift)):

            # Run 21CMMC (the 21cmFAST code) to perform the reionisation simulation and recover the neutral fraction and 21cm PS
            command = "./drive_21cmMC_streamlined %g %s 1"%(self.Redshift[i],StringArgument) 
            os.system(command)

            # Read in the neutral fraction and 21cm PS for this parameter set and redshift
            k_values_estimate = np.loadtxt('delTps_estimate_%s.txt'%(StringArgument_other), usecols=(0,))
            PS_values_estimate = np.loadtxt('delTps_estimate_%s.txt'%(StringArgument_other), usecols=(1,))
            nf_value = np.loadtxt('NeutralFraction_%s.txt'%(StringArgument_other), usecols=(0,))

            nf_vals[i] = nf_value

            if nf_value == 0.:
                # neutral fraction of zero, means no 21cm PS. recover a large number to steer clear of this region of parameter space
                total_sum += 100.*sum(np.square((self.PS_values[i])/self.PS_Error[i]))
            else:
                # compute the chi-squared value compated to the mock observation for this redshift

                # Interpolating the mock and error PS in log space
                PSError_Spline = interpolate.splrep(self.Error_k_values[i],np.log10(self.PS_Error[i]),s=0)
                MockPS_Spline = interpolate.splrep(self.k_values[i],np.log10(self.PS_values[i]),s=0)

                for ii in range(len(k_values_estimate)):
                    if k_values_estimate[ii] >= self.Foreground_cut and k_values_estimate[ii] < self.Shot_Noise_cut: 

                        # As the interpolation is performed in log space
                        MockPS_val = 10**(interpolate.splev(k_values_estimate[ii],MockPS_Spline,der=0))
                        ErrorPS_val = 10**(interpolate.splev(k_values_estimate[ii],PSError_Spline,der=0))

                        total_sum += np.square((MockPS_val - PS_values_estimate[ii])/(np.sqrt(ErrorPS_val**2. + (self.ModUncert*PS_values_estimate[ii])**2.))) 
       
            # remove the temporary files
            command = "rm delTps_estimate_%s.txt"%(StringArgument_other)
            os.system(command)
            command = "rm NeutralFraction_%s.txt"%(StringArgument_other)
            os.system(command)

        if self.PriorLegend['PlanckPrior'] is True and len(self.Redshift) > 2:

            # Mean and one sigma errors for the Planck constraints
            # The Planck prior is modelled as a Gaussian: tau = 0.058 \pm 0.012 (https://arxiv.org/abs/1605.03507)
            PlanckTau_Mean = 0.058
            PlanckTau_OneSigma = 0.012

            # Simple linear extrapolation of the redshift range provided by the user, to be able to estimate the optical depth
            nZinterp = 15

            # The minimum of the extrapolation is chosen to 5.9, to correspond to the McGreer et al. prior on the IGM neutral fraction.
            # The maximum is chosed to be z = 18., which is arbitrary.
            ZExtrap_min = 5.9
            ZExtrap_max = 18.0

            ZExtrapVals = np.zeros(nZinterp)
            XHI_ExtrapVals = np.zeros(nZinterp)

            # Perform only a linear interpolation/extrapolation
            order = 1

            # The linear interpolation/extrapolation function, taking as input the redshifts supplied by the user and the corresponding neutral fractions
            # recovered for the specific EoR parameter set
            LinearInterpolationFunction = InterpolatedUnivariateSpline(self.Redshift, nf_vals, k=order)

            for i in range(nZinterp):
                ZExtrapVals[i] = ZExtrap_min + (ZExtrap_max - ZExtrap_min)*float(i)/(nZinterp - 1)
    
                XHI_ExtrapVals[i] = LinearInterpolationFunction(ZExtrapVals[i])
            
                # Ensure that the neutral fraction does not exceed unity, or go negative
                if XHI_ExtrapVals[i] > 1.0:
                    XHI_ExtrapVals[i] = 1.0
                if XHI_ExtrapVals[i] < 0.0:
                    XHI_ExtrapVals[i] = 0.0

            # Set up the arguments for calculating the estimate of the optical depth. Once again, performed using command line code.
            separator_Planck = " "
            seq_Planck = []
            for i in range(nZinterp):
                seq_Planck.append("%s"%(ZExtrapVals[i])) 
                seq_Planck.append("%s"%(XHI_ExtrapVals[i]))    

            StringArgument_Planck = string.join(seq_Planck,separator_Planck)

            # Perform the computation of tau
            command = './ComputingTau_e %s %s %s'%(Individual_ID,Decimal(repr(params[0])).quantize(SIXPLACES),StringArgument_Planck)
            os.system(command)

            # Read tau from file
            tau_value = np.loadtxt('Tau_e_%s_%s.txt'%(Individual_ID,Decimal(repr(params[0])).quantize(SIXPLACES)), usecols=(0,))

            # remove the temporary files
            command = "rm Tau_e_%s_%s.txt"%(Individual_ID,Decimal(repr(params[0])).quantize(SIXPLACES))
            os.system(command)

            # As the likelihood is computed in log space, the addition of the prior is added linearly to the existing chi^2 likelihood
            total_sum = total_sum + np.square( ( PlanckTau_Mean - tau_value )/(PlanckTau_OneSigma) )

        if self.PriorLegend['McGreerPrior'] is True:

            # Mean and one sigma errors for the McGreer et al. constraints
            # Modelled as a flat, unity prior at x_HI <= 0.06, and a one sided Gaussian at x_HI > 0.06 ( Gaussian of mean 0.06 and one sigma of 0.05 )
            McGreer_Mean = 0.06
            McGreer_OneSigma = 0.05            

            if McGreer_Redshift in self.Redshift:

                for i in range(len(self.Redshift)):
                    if self.Redshift[i] == McGreer_Redshift:                        
                        McGreer_NF = nf_vals[i]

                if McGreer_NF > 1.:
                    McGreer_NF = 1.
                if McGreer_NF < 0.:
                    McGreer_NF = 0.

                # As the likelihood is computed in log space, the addition of the prior is added linearly to the existing chi^2 likelihood
                if McGreer_NF <= 0.06:
                    total_sum = total_sum + 0.0 # Add zero, as we assume flat (unity) probability at x_HI <= 0.06 (as it is a lower limit)
                else:
                    total_sum = total_sum + np.square( ( McGreer_Mean - McGreer_NF )/(McGreer_OneSigma) )


            elif len(self.Redshift) > 2:

                # Perform only a linear interpolation/extrapolation
                order = 1

                # The linear interpolation/extrapolation function, taking as input the redshifts supplied by the user and the corresponding neutral fractions
                # recovered for the specific EoR parameter set
                LinearInterpolationFunction = InterpolatedUnivariateSpline(self.Redshift, nf_vals, k=order)

                McGreer_NF = LinearInterpolationFunction(McGreer_Redshift)

                if McGreer_NF > 1.:
                    McGreer_NF = 1.
                if McGreer_NF < 0.:
                    McGreer_NF = 0.

                # As the likelihood is computed in log space, the addition of the prior is added linearly to the existing chi^2 likelihood
                if McGreer_NF <= 0.06:
                    total_sum = total_sum + 0.0 # Add zero, as we assume flat (unity) probability at x_HI <= 0.06 (as it is a lower limit)
                else:
                    total_sum = total_sum + np.square( ( McGreer_Mean - McGreer_NF )/(McGreer_OneSigma) )

        if self.PriorLegend['GreigPrior'] is True:

            # Interpolate the QSO damping wing PDF
            spline_QSODampingPDF = interpolate.splrep(self.NFValsQSO,self.PDFValsQSO,s=0)

            if QSO_Redshift in self.Redshift:

                for i in range(len(self.Redshift)):
                    if self.Redshift[i] == QSO_Redshift:                        
                        NF_QSO = nf_vals[i]

                # Ensure that the neutral fraction does not exceed unity, or go negative
                if NF_QSO > 1.0:
                    NF_QSO = 1.0
                if NF_QSO < 0.0:
                    NF_QSO = 0.0

                QSO_Prob = interpolate.splev(NF_QSO,spline_QSODampingPDF,der=0)

                # Interpolating the PDF from the QSO damping wing might cause small negative values at the edges (i.e. x_HI ~ 0 or ~1)
                # In case it is zero, or negative, set it to a very small non zero number (we take the log of this value, it cannot be zero)
                if QSO_Prob <= 0.0:
                    QSO_Prob = 0.000006

                # We work with the log-likelihood, therefore convert the IGM Damping wing PDF to log space
                QSO_Prob = -2.*np.log(QSO_Prob)

                total_sum = total_sum + QSO_Prob

            elif len(self.Redshift) > 2:
            
                order = 1

                # Check the redshift range input by the user to determine whether to interpolate or extrapolate the IGM neutral fraction to the QSO redshift
                if QSO_Redshift < np.amin(self.Redshift):
                    # The QSO redshift is outside the range set by the user. Need to extrapolate the reionisation history to obtain the neutral fraction at the QSO redshift

                    # The linear interpolation/extrapolation function, taking as input the redshifts supplied by the user and the corresponding neutral fractions
                    # recovered for the specific EoR parameter set
                    LinearInterpolationFunction = InterpolatedUnivariateSpline(self.Redshift, nf_vals, k=order)

                    NF_QSO = LinearInterpolationFunction(QSO_Redshift)
                            
                else:
                    # The QSO redshift is within the range set by the user. Can interpolate the reionisation history to obtain the neutral fraction at the QSO redshift

                    spline_reionisationhistory = interpolate.splrep(self.Redshift,nf_vals,s=0)

                    NF_QSO = interpolate.splev(QSO_Redshift,spline_reionisationhistory,der=0)

                # Ensure that the neutral fraction does not exceed unity, or go negative
                if NF_QSO > 1.0:
                    NF_QSO = 1.0
                if NF_QSO < 0.0:
                    NF_QSO = 0.0

                QSO_Prob = interpolate.splev(NF_QSO,spline_QSODampingPDF,der=0)

                # Interpolating the PDF from the QSO damping wing might cause small negative values at the edges (i.e. x_HI ~ 0 or ~1)
                # In case it is zero, or negative, set it to a very small non zero number (we take the log of this value, it cannot be zero)
                if QSO_Prob <= 0.0:
                    QSO_Prob = 0.000006

                # We work with the log-likelihood, therefore convert the IGM Damping wing PDF to log space
                QSO_Prob = -2.*np.log(QSO_Prob)

                total_sum = total_sum + QSO_Prob

        return -0.5*total_sum,nf_vals

    def computeLikelihood(self, ctx):

        return self.Likelihood(ctx)

    def setup(self):
        print "Likelihood Fitting for 21cm Fast" 
