
import numpy as np
import os
from decimal import *
import multiprocessing
import itertools

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')

def ICposition(start_positions_individual,lowerbounds,upperbounds,FiducialValues,ParamWidths):

    # Test that the chosen initial conditions are within the initial parameter bounds, if not generate ones that are
    for j in range(len(lowerbounds)):
        if (start_positions_individual[j] < lowerbounds[j] or start_positions_individual[j] > upperbounds[j]):
            new_start_parameter_logic = False
            while new_start_parameter_logic == False:
                new_start_parameter = FiducialValues[j]+np.random.normal(size=1.0)*ParamWidths[j]
                if (new_start_parameter > lowerbounds[j] and new_start_parameter < upperbounds[j]):
                    new_start_parameter_logic = True

            start_positions_individual[j] = new_start_parameter            

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

        InputValues = []

        if self.sampler.param_legend['ALPHA'] is True:   
            InputValues.append(self.sampler.FiducialParams['ALPHA'])

        if self.sampler.param_legend['ZETA'] is True:   
            InputValues.append(self.sampler.FiducialParams['ZETA'])    

        if self.sampler.param_legend['MFP'] is True:   
            InputValues.append(self.sampler.FiducialParams['MFP'])    

        if self.sampler.param_legend['TVIR_MIN'] is True:   
            InputValues.append(self.sampler.FiducialParams['TVIR_MIN'])    

        if self.sampler.param_legend['L_X'] is True:   
            InputValues.append(self.sampler.FiducialParams['L_X'])    

        if self.sampler.param_legend['NU_X_THRESH'] is True:   
            InputValues.append(self.sampler.FiducialParams['NU_X_THRESH'])   

        if self.sampler.param_legend['X_RAY_SPEC_INDEX'] is True:   
            InputValues.append(self.sampler.FiducialParams['X_RAY_SPEC_INDEX'])

        if self.sampler.param_legend['SIGMA_8'] is True:
            InputValues.append(self.sampler.FiducialParams['SIGMA_8'])

        if self.sampler.param_legend['littleh'] is True:
            InputValues.append(self.sampler.FiducialParams['littleh'])

        if self.sampler.param_legend['OMEGA_M'] is True:
            InputValues.append(self.sampler.FiducialParams['OMEGA_M'])

        if self.sampler.param_legend['OMEGA_b'] is True:
            InputValues.append(self.sampler.FiducialParams['OMEGA_b'])

        if self.sampler.param_legend['NS'] is True:
            InputValues.append(self.sampler.FiducialParams['NS'])

        ParamWidths = []
        for i in range(len(InputValues)):
            ParamWidths.append( (self.sampler.upperbounds[i] - self.sampler.lowerbounds[i])/3. )

        print('Generate Start Positions')
        start_positions = [InputValues+np.random.normal(size=self.sampler.paramCount)*ParamWidths for i in xrange(self.sampler.nwalkers)]

        pool = multiprocessing.Pool(self.sampler.threadCount)

        M = pool.map

        returned_list = list(M(ICposition_star,itertools.izip(start_positions, itertools.repeat(self.sampler.lowerbounds), itertools.repeat(self.sampler.upperbounds), 
                            itertools.repeat(InputValues), itertools.repeat(ParamWidths))))

        print('Start Positions Generated')    
        return returned_list
    
    def __str__(self, *args, **kwargs):
        return "SampleBallPositionGenerator"
