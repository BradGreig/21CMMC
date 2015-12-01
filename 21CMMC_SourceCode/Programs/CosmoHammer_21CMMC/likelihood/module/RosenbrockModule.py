import numpy as np

class RosenbrockModule(object):

    def __init__(self):
        self.a1 = 100.0
        self.a2 = 20.0

    def computeLikelihood(self, ctx):
        p = ctx.getParams()
        return -(self.a1 * (p[1]-p[0]**2)**2 + (1-p[0])**2)/self.a2
    
    def setup(self):
        print "Rosenbrock setup"