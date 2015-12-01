#!/usr/bin/env python

from cosmoHammer.likelihood.ChainConstants import *
import numpy as np

class PseudoCambCoreModule(object):
    """
    Core Module for the delegation of the computation of the cmb power spectrum to the cambWrapperManager.

    """
    
    def __init__(self, filePath=None):
        self.filePath = filePath
        self.cl_tt = np.zeros(1199)
        self.cl_te = np.zeros(1199)
        self.cl_ee = np.zeros(1199)
        self.cl_bb = np.zeros(1199)
        
    def __call__(self, ctx):
        self.computeCmbPowerSpectrum(ctx)

    def computeCmbPowerSpectrum(self, ctx):
        """
        Calls the native code to compute the cmb power spectrum
        """
        ctx.add(CL_TT_KEY, self.cl_tt)
        ctx.add(CL_TE_KEY, self.cl_te)
        ctx.add(CL_EE_KEY, self.cl_ee)
        ctx.add(CL_BB_KEY, self.cl_bb)

    def setup(self):
        """
        Sets up the cmb likelihood wrapper
        """
        data = np.loadtxt(self.filePath)
        self.copyArray(data, self.cl_tt, 1)
        self.copyArray(data, self.cl_te, 4)
        self.copyArray(data, self.cl_ee, 2)
        self.copyArray(data, self.cl_bb, 3)

    def copyArray(self, src, trg, index):
        lenght= min(len(src),len(trg))
        for i in xrange(0,lenght):
            trg[i]=src[i,index]
