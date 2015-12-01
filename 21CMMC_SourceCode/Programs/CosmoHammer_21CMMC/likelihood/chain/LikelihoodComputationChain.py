#!/usr/bin/env python
from CosmoHammer_21CMMC.likelihood.ChainContext import ChainContext
import numpy as np
from collections import deque
import logging

class LikelihoodComputationChain(object):
    """
    Implementation of a likelihood computation chain.
    """

    def __init__(self, min=None, max=None):
        """
        Constructor for the likelihood chain

        :param min: array 
            lower bound for the parameters
        :param max: array
            upper bound for the parameters

        """
        self.min = min
        self.max = max
        self._likelihoodModules = deque();
        self._coreModules = deque();
            

    def getCoreModules(self):
        """pointer to the likelihood module list """
        return self._coreModules

    def getLikelihoodModules(self):
        """pointer to the core module list """
        return self._likelihoodModules

    def addLikelihoodModule(self, module):
        """
        adds a module to the likelihood module list
        
        :param module: callable
            the callable module to add for the likelihood computation
        """
        self.getLikelihoodModules().append(module)
        
    def addCoreModule(self, module):
        """
        adds a module to the likelihood module list
        
        :param module: callable
            the callable module to add for the computation of the data
        """
        self.getCoreModules().append(module)
        
        
    def isValid(self, p):
        """
        checks if the given parameters are valid 
        """
        if(self.min is not None):
            for i in xrange(len(p)):
                if (p[i]<self.min[i]):
                    logging.debug("Params out of bounds i="+str(i)+" params "+str(p))
                    return False
        
        if(self.max is not None):
            for i in xrange(len(p)):
                if (p[i]>self.max[i]):
                    logging.debug("Params out of bounds i="+str(i)+" params "+str(p))
                    return False
        
        return True

    
    def setup(self):
        """sets up the chain and its modules """
        for cModule in self.getCoreModules():
            cModule.setup()
            
        for cModule in self.getLikelihoodModules():
            cModule.setup()
            
    
    def __call__(self, p):
        """
        Computes the log likelihood by calling all the core and likelihood modules.
        
        :return: the current likelihood and a dict with additional data
        """
        if not self.isValid(p):
            return -np.inf, []
        
        ctx = self.createChainContext(p)

        self.invokeCoreModules(ctx)

        likelihood = self.computeLikelihoods(ctx)
        return likelihood, ctx.getData()
    
    def createChainContext(self, p):
        """
        Returns a new instance of a chain context 
        """
        return ChainContext(self, p)
    
    def invokeCoreModules(self, ctx):
        """
        Iterates thru the core modules and invokes them
        """
        for cModule in self.getCoreModules():
            self.invokeCoreModule(cModule, ctx)
            
    
    def invokeCoreModule(self, coreModule, ctx):
        """
        Invokes the given module with the given ChainContext
        """
        coreModule(ctx)
        
    
    def computeLikelihoods(self, ctx):
        """
        Computes the likelihoods by iterating thru all the modules.
        Sums up the log likelihoods.
        """
        likelihood = 0      
        stuff = 0  
        
        for lModule in self.getLikelihoodModules():            
            likelihood_i, stuff_i = self.invokeLikelihoodModule(lModule, ctx)
            likelihood += likelihood_i
            stuff += stuff_i

        return likelihood, stuff
    
    def invokeLikelihoodModule(self, likelihoodModule, ctx):
        """
        Invokes the given module with the given ChainContext
        """
        val1, val2 = likelihoodModule.computeLikelihood(ctx)
        return val1, val2
    
    def __str__(self, *args, **kwargs):
        s = "Core Modules: \n  "
        s = s + "\n  ".join(map(lambda o: type(o).__name__, self.getCoreModules()))

        s = s + "\nLikelihood Modules: \n  "
        s = s + "\n  ".join(map(lambda o: type(o).__name__, self.getLikelihoodModules()))
        return s
