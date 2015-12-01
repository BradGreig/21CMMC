
import numpy as np
import os
from decimal import *

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')

class SampleBallPositionGenerator(object):
    """
        Generates samples in a very thight n-dimensional ball 
    """
    
    def __init__(self):
        """
            default constructor
        """
        pass

    def setup(self, sampler):
        """
            setup the generator
        """
        self.sampler = sampler
    
    def generate(self):
        """
            generates the positions
        """
        print('Generate Start Positions')
        start_positions = [self.sampler.paramValues+3.*np.random.normal(size=self.sampler.paramCount)*self.sampler.paramWidths for i in xrange(self.sampler.nwalkers)]
#        start_positions = [np.random.uniform(low=self.sampler.lowerbounds,high=self.sampler.upperbounds,size=self.sampler.paramCount) for i in xrange(self.sampler.nwalkers)]

        for i in xrange(self.sampler.nwalkers):            
            for j in range(len(self.sampler.lowerbounds)):
                if (start_positions[i][j] < self.sampler.lowerbounds[j] or start_positions[i][j] > self.sampler.upperbounds[j]):
                    new_start_parameter_logic = False
                    while new_start_parameter_logic == False:
                        new_start_parameter = self.sampler.paramValues[j]+3.*np.random.normal(size=1.0)*self.sampler.paramWidths[j]
#                        new_start_parameter = np.random.uniform(low=self.sampler.lowerbounds[j],high=self.sampler.upperbounds[j],size=1.0)
                        if (new_start_parameter > self.sampler.lowerbounds[j] and new_start_parameter < self.sampler.upperbounds[j]):
                            new_start_parameter_logic = True

                    start_positions[i][j] = new_start_parameter            

            new_start_parameter_logic = False
            while new_start_parameter_logic == False:

                if self.sampler.Redshift_prior is not None:

#                    command = "../../../ ./drive_21cmMC %g %g %g %g 0"%(self.sampler.Redshift_prior,start_positions[i][0],start_positions[i][1],start_positions[i][2])
#                    command = "./drive_21cmMC %g %g %g %g 0"%(self.sampler.Redshift_prior,start_positions[i][0],start_positions[i][1],start_positions[i][2])
                    command = "./drive_21cmMC %g %s %s %s 1.0 1 0 0 0 0"%(self.sampler.Redshift_prior,Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)) 
#                    print(command)        
                    os.system(command)

#                    nf_value = np.loadtxt('NeutralFraction.txt', usecols=(0,))
                    nf_value = np.loadtxt('NeutralFraction_%s_%s_%s.txt'%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.0).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)), usecols=(0,))

                    if nf_value < 0.05:                
                        
                        new_start_parameter_logic = False

                    else:
                            
#                        command = "./drive_21cmMC %g %g %g %g 1"%(self.sampler.Redshift,start_positions[i][0],start_positions[i][1],start_positions[i][2])
                        command = "./drive_21cmMC %g %s %s %s 1.0 1 0 0 0 0"%(self.sampler.Redshift,Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES))
#                        print(command)
                        os.system(command)                            

                        k_values_estimate = np.loadtxt('delTps_estimate_%s_%s_%s.txt'%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.0).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)), usecols=(0,))
                        PS_values_estimate = np.loadtxt('delTps_estimate_%s_%s_%s.txt'%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.0).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)), usecols=(1,))
                        nf_value = np.loadtxt('NeutralFraction_%s_%s_%s.txt'%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.0).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)), usecols=(0,))

                        if nf_value == 0.:
                            new_start_parameter_logic = False
                        else:
                            new_start_parameter_logic = True

                else:

#                    command = "./drive_21cmMC %g %g %g %g 1"%(self.sampler.Redshift,start_positions[i][0],start_positions[i][1],start_positions[i][2])
                    command = "./drive_21cmMC %g %s %s %s 1.0 1 0 0 0 0"%(self.sampler.Redshift,Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES))
#                    print(command)
                    os.system(command)

                    k_values_estimate = np.loadtxt('delTps_estimate_%s_%s_%s.txt'%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.0).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)), usecols=(0,))
                    PS_values_estimate = np.loadtxt('delTps_estimate_%s_%s_%s.txt'%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.0).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)), usecols=(1,))
                    nf_value = np.loadtxt('NeutralFraction_%s_%s_%s.txt'%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.0).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES)), usecols=(0,))    

                    if nf_value == 0.:
                        new_start_parameter_logic = False
                    else:
                        new_start_parameter_logic = True

                command = "rm delTps_estimate_%s_%s_%s.txt"%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES))
                os.system(command)
                command = "rm NeutralFraction_%s_%s_%s.txt"%(Decimal(start_positions[i][0]).quantize(FOURPLACES),Decimal(15.).quantize(FOURPLACES),Decimal(10**(start_positions[i][1])).quantize(FOURPLACES))
                os.system(command)

                if new_start_parameter_logic == False:
                    for j in range(len(self.sampler.lowerbounds)):
                        start_parameter_logic_brandnew = False
                        while start_parameter_logic_brandnew == False:
                            new_start_parameter = self.sampler.paramValues[j]+3.*np.random.normal(size=1.0)*self.sampler.paramWidths[j]
                            if (new_start_parameter > self.sampler.lowerbounds[j] and new_start_parameter < self.sampler.upperbounds[j]):
                                start_positions[i][j] = new_start_parameter
                                start_parameter_logic_brandnew = True

        print('Start Positions Generated')    
        return start_positions
    
    def __str__(self, *args, **kwargs):
        return "SampleBallPositionGenerator"
