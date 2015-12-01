#!/usr/bin/env python

import numpy as np

class PseudoCmbModule(object):
    """
    Chain for computing the likelihood of a multivariante gaussian distribution
    """
    def __init__(self, icov, mu, params):
        self.icov = icov
        self.mu = mu
        self.min = params[:,1]
        self.max = params[:,2]
        self.a = self.min[-1]
        self.b = self.max[-1]
    
    
    def computeLikelihood(self, ctx):
        x = ctx.getParams()
        if not self.isValid(x):
            return -np.inf

        diff = x[:6]-self.mu
        
        lnprob = -np.dot(diff,np.dot(self.icov,diff))/2.0

        lnprob -= np.log(self.b-self.a)
        return lnprob
    
    
    def isValid(self, p):
        """checks if the given parameters are valid """
        for i in xrange(len(p)):
            if (p[i]<self.min[i] or p[i]>self.max[i]):
                return False
            
        return True    

    def setup(self):
        print "Pseudo cmb setup"
            
        
    