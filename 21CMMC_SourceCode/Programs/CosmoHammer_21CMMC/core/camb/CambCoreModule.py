from cambWrapper import cambWrapperManager
from cosmoHammer.likelihood.ChainConstants import *

class CambCoreModule(object):
    """
    Core Module for the delegation of the computation of the cmb power spectrum to the cambWrapperManager.
    
    :param filePath: (optional)
        path to the "param.ini" to be used
    """
    def __init__(self, filePath=None):
        """
        Constructor
        """
        self.filePath = filePath
        
    def __call__(self, ctx):
        self.computeCmbPowerSpectrum(ctx)

    def computeCmbPowerSpectrum(self, ctx):
        """
        Calls the native code to compute the cmb power spectrum
        """
        p1 = ctx.getParams()[0:6]
        cl_tt,cl_te,cl_ee,cl_bb = cambWrapperManager.computecmbpowerspectrum(p1)

        ctx.add(CL_TT_KEY, cl_tt)
        ctx.add(CL_TE_KEY, cl_te)
        ctx.add(CL_EE_KEY, cl_ee)
        ctx.add(CL_BB_KEY, cl_bb)

    def setup(self):
        """
        Sets up the cmb likelihood wrapper
        """
        cambWrapperManager.setup(self.filePath)
