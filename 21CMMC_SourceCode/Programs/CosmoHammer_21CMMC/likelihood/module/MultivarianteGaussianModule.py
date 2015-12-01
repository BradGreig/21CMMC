#!/usr/bin/env python

import numpy as np

class MultivarianteGaussianModule(object):
    """
    Chain for computing the likelihood of a multivariante gaussian distribution
    """
    def __init__(self, icov, mu):
        self.icov = icov
        self.mu = mu
    
    
    def computeLikelihood(self, ctx):
        x = ctx.getParams()
        diff = x-self.mu
        return -np.dot(diff,np.dot(self.icov,diff))/2.0

    def setup(self):
        print "Multivariante Gaussian setup"
            
        
    