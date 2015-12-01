
import numpy as np
import os
from decimal import *
import multiprocessing
import itertools

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')

def ICposition(start_positions_individual,lowerbounds,upperbounds,paramValues,paramWidths,Redshift):

    np.random.seed()

    # Generate a random number per thread as a unique identifier
    random_number = np.random.normal(size=1.0)

    Individual_ID = Decimal(repr(random_number[0])).quantize(FOURPLACES)

    # Test that the chosen initial conditions are within the initial parameter bounds, if not generate ones that are
    for j in range(len(lowerbounds)):
        if (start_positions_individual[j] < lowerbounds[j] or start_positions_individual[j] > upperbounds[j]):
            new_start_parameter_logic = False
            while new_start_parameter_logic == False:
                new_start_parameter = paramValues[j]+3.*np.random.normal(size=1.0)*paramWidths[j]
                if (new_start_parameter > lowerbounds[j] and new_start_parameter < upperbounds[j]):
                    new_start_parameter_logic = True

            start_positions_individual[j] = new_start_parameter            

    # This step is not really necessary and can be removed
    # Here, I check that the initial locations are not in a region of parameter space where the box is completely ionised (under the assumption that we 
    # are at a redshift with a measured 21 cm PS). The sampler should leave this position anyway, but in some cases (double power law model) the gradient of
    # the parameter space was too sharp to escape (without a sufficiently large (and rare) jump, therefore this improved the sampling in those models)
    new_start_parameter_logic = False
    while new_start_parameter_logic == False:

        neutral_fraction_filename = "NeutralFraction_%s_%s_%s_%s.txt"%(Individual_ID, Decimal(repr(start_positions_individual[0])).quantize(FOURPLACES),
                        Decimal(repr(start_positions_individual[1])).quantize(FOURPLACES), Decimal(repr(start_positions_individual[2])).quantize(FOURPLACES))

        PS_filename = "delTps_estimate_%s_%s_%s_%s.txt"%(Individual_ID,Decimal(repr(start_positions_individual[0])).quantize(FOURPLACES),
                Decimal(repr(start_positions_individual[1])).quantize(FOURPLACES),Decimal(repr(start_positions_individual[2])).quantize(FOURPLACES))

        # Run the 21cmMC driver
        command = "./drive_21cmMC_streamlined %g %s %s %s %s 1"%(Redshift,Individual_ID,Decimal(repr(start_positions_individual[0])).quantize(FOURPLACES),
                Decimal(repr(start_positions_individual[1])).quantize(FOURPLACES),Decimal(repr(start_positions_individual[2])).quantize(FOURPLACES))
        os.system(command)

        k_values_estimate = np.loadtxt('%s'%(PS_filename), usecols=(0,))
        PS_values_estimate = np.loadtxt('%s'%(PS_filename), usecols=(1,))
        nf_value = np.loadtxt('%s'%(neutral_fraction_filename), usecols=(0,))    

        # Test the neutral fraction value
        if nf_value == 0.:
            new_start_parameter_logic = False
        else:
            new_start_parameter_logic = True

        command = "rm %s"%(PS_filename)
        os.system(command)
        command = "rm %s"%(neutral_fraction_filename)
        os.system(command)

        # If we are in a region completely ionised, find a new set of walker positions
        if new_start_parameter_logic == False:
            for j in range(len(lowerbounds)):
                start_parameter_logic_brandnew = False
                while start_parameter_logic_brandnew == False:
                    new_start_parameter = paramValues[j]+3.*np.random.normal(size=1.0)*paramWidths[j]
                    if (new_start_parameter > lowerbounds[j] and new_start_parameter < upperbounds[j]):
                        start_positions_individual[j] = new_start_parameter
                        start_parameter_logic_brandnew = True

    return start_positions_individual

def ICposition_star(all_arguments):

    return ICposition(*all_arguments)

class UniformPosition(object):
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
        # Generate the initial walker positions, first checking they are within the parameter bounds.

        print('Generate Start Positions')
        start_positions = [self.sampler.paramValues+3.*np.random.normal(size=self.sampler.paramCount)*self.sampler.paramWidths for i in xrange(self.sampler.nwalkers)]

        pool = multiprocessing.Pool(self.sampler.threadCount)

        M = pool.map

        returned_list = list(M(ICposition_star,itertools.izip(start_positions, itertools.repeat(self.sampler.lowerbounds), itertools.repeat(self.sampler.upperbounds), itertools.repeat(self.sampler.paramValues), 
        itertools.repeat(self.sampler.paramWidths), itertools.repeat(self.sampler.Redshift))))

        print('Start Positions Generated')    
        return returned_list
    
    def __str__(self, *args, **kwargs):
        return "SampleBallPositionGenerator"
