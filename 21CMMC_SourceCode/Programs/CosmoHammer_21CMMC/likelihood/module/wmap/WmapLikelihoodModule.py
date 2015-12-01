#!/usr/bin/env python
from cosmoHammer.likelihood.ChainConstants import *
import importlib

class WmapLikelihoodModule(object):
    
    """
    Module for the wmap likelihood computation
    
    :param year: the desired WMAP year
    :param path: the path to the WMAP data
    
    """
    def __init__(self, year="9", path=None):
        """
        Constructor
        """
        self.path = path
        factory = WmapWrapperFactory()
        self.wmapWrapper = factory.createWrapper(year)
    
    

    def computeLikelihood(self, ctx):
        """
        call the native code to compute log likelihood
        """
        cl_tt =ctx.get(CL_TT_KEY)
        cl_te =ctx.get(CL_TE_KEY)
        cl_ee =ctx.get(CL_EE_KEY)
        cl_bb =ctx.get(CL_BB_KEY)
        
        loglike = -self.wmapWrapper.computewmaplikelihood(cl_tt,cl_te,cl_ee,cl_bb)
        return loglike

    def setup(self):
        """
        Sets up the wmap likelihood wrapper
        """
        self.wmapWrapper.setup(self.path)



class WmapWrapperFactory(object):
    
    _supportedVersions = ["3", "5", "7", "9"]
    _importStmt = "wmap{0}Wrapper.wmapWrapperManager"
    
    def createWrapper(self, version):
        if(version in self._supportedVersions):
            stmt = self._importStmt.format(version)
            return importlib.import_module(stmt)
        else:
            raise WmapVersionUnsupportedError(version)
        
        
class WmapVersionUnsupportedError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return "The given WmapWrapper version {0} is not supported!".format(self.value)