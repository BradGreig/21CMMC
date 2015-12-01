#!/usr/bin/env python
"""

"""
from .WmapLikelihoodModule import WmapLikelihoodModule
from cosmoHammer.likelihood.ChainConstants import *
from pkg_resources import resource_filename
import cosmoHammer
import numpy as np

FILE_PATH  = "data/WMAP_SZ_VBand.dat"

class WmapExtLikelihoodModule(WmapLikelihoodModule):
    """
    Extension for the WMAP likelihood computation using the sz template to post process the power spectrum
    
    :param year: the desired WMAP year
    :param path: the path to the WMAP data
    :param aszIndex: Index of the asz parameter in the walker position sequence
    """
    def __init__(self, year="9", aszIndex=6, path=None):
        """
        Constructor
        """
        super(WmapExtLikelihoodModule, self).__init__(year, path)
        self.aszIndex = aszIndex

    def setup(self):
        super(WmapExtLikelihoodModule, self).setup()
        
        path = resource_filename(cosmoHammer.__name__, FILE_PATH)
        
        print "Loading sz template from: " + path
        
        self.sz = np.loadtxt(path)[:,1]
        
        
    def computeLikelihood(self, ctx):
        self.postProcessPowerSpectrum(ctx)
        
        return super(WmapExtLikelihoodModule, self).computeLikelihood(ctx)
        
    def postProcessPowerSpectrum(self, ctx):
        """Add data from sz template times Asz"""
        cl_tt = ctx.get(CL_TT_KEY)
        p = ctx.getParams()
        
        l = len(cl_tt)
        cl_tt[0:l] = cl_tt[0:l] + p[self.aszIndex]*self.sz[0:l]
        