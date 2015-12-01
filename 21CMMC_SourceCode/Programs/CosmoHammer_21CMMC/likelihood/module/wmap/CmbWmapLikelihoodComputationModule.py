#!/usr/bin/env python

from wmap7Wrapper import wmapWrapperManager
from cosmoHammer.likelihood.ChainConstants import *

class CmbWmapLikelihoodComputationModule(object):
    
    """
    Module for the wmap likelihood computation
    
    :param path: the path to the WMAP data
    
    """
    def __init__(self, path=None):
        """
        Constructor
        """
        self.path = path
    
    

    def computeLikelihood(self, ctx):
        """
        call the native code to compute log likelihood
        """
        cl_tt =ctx.get(CL_TT_KEY)
        cl_te =ctx.get(CL_TE_KEY)
        cl_ee =ctx.get(CL_EE_KEY)
        cl_bb =ctx.get(CL_BB_KEY)
        
        loglike = -wmapWrapperManager.computewmaplikelihood(cl_tt,cl_te,cl_ee,cl_bb)
        return loglike

    def setup(self):
        """
        Sets up the cmb likelihood wrapper
        """

        wmapWrapperManager.setup(self.path)
