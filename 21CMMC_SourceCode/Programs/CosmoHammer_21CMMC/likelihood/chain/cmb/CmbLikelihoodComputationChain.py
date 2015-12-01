#!/usr/bin/env python

from cosmoHammer.likelihood.chain.LikelihoodComputationChain import LikelihoodComputationChain
from cosmoHammer.core.camb.CambCoreModule import CambCoreModule
from cosmoHammer.likelihood.module.wmap.WmapLikelihoodModule import WmapLikelihoodModule
from cosmoHammer.likelihood.module.wmap.WmapExtLikelihoodModule import WmapExtLikelihoodModule

class CmbLikelihoodComputationChain(LikelihoodComputationChain):
    """
        Chain to compute the likelihood using cmb power spectrum. 
        Configures the likelihood chain by default with the 
        CambCoreModule for the theory prediction and either CmbWmapLikelihoodComputationModule
        CmbWmapExtLikelihoodComputationModule for the likelihood computation
    """
    
    def __init__(self, min, max, version, path=None, aszIndex=None):
        """
            Arguments:
            min: lower bound for the parameters
            max: upper bound for the parameters
            path: path to the wmap data
            aszIndex: Index of the asz parameter in the walker position sequence.
        """
        super(CmbLikelihoodComputationChain, self).__init__(min, max)
        
        self.addCoreModule(CambCoreModule())
        
        if(aszIndex is None):
            wmapLikelihood = WmapLikelihoodModule(version=version, path=path);
        else:
            wmapLikelihood = WmapExtLikelihoodModule(version=version, 
                        path=path, aszIndex=aszIndex) 
        
        self.addLikelihoodModule(wmapLikelihood)
        
    
