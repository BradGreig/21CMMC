from .CosmoHammerSampler import CosmoHammerSampler
from collections import namedtuple
from cosmoHammer.util.SampleFileUtil import SampleFileUtil
from mpi4py import MPI
import emcee
import itertools


class MpiCosmoHammerSampler(CosmoHammerSampler):
    """
    A sampler implementation extending the regular sampler in order to allow for distributing 
    the computation with MPI.

    :param kwargs:  
        key word arguments passed to the CosmoHammerSampler
    
    """
    def __init__(self, **kwargs):
        """
        CosmoHammer sampler implementation
        
        """
        self._rank = MPI.COMM_WORLD.Get_rank()
        
        super(MpiCosmoHammerSampler, self).__init__(**kwargs)
        
        
        self.M = self._getMapFunction()
        
    def _getMapFunction(self):
        """
        Returns the build in map function
        """
        return map
    
    def createSampleFileUtil(self):
        """
        Returns a new instance of a File Util
        """
        return SampleFileUtil(self.filePrefix, self.isMaster(), reuseBurnin=self.reuseBurnin)
    
       
    def sampleBurnin(self, p0):
        """
        Starts the sampling process. The master node (mpi rank = 0) persists the result to the disk
        """
        p0 = self.mpiBCast(p0)
        
        self.log("MPI Process rank "+ str(self._rank)+" starts sampling")
        return super(MpiCosmoHammerSampler, self).sampleBurnin(p0);
   
    def sample(self, burninPos, burninProb, burninRstate, datas):
        """
        Starts the sampling process. The master node (mpi rank = 0) persists the result to the disk
        """
        burninPos = self.mpiBCast(burninPos)
        burninProb = self.mpiBCast(burninProb)
        burninRstate = self.mpiBCast(burninRstate)
        
        self.log("MPI Process rank "+ str(self._rank)+" starts sampling")
        super(MpiCosmoHammerSampler, self).sample(burninPos, burninProb, burninRstate, datas);

            
    def loadBurnin(self):
        """
        loads the burn in form the file system
        """
        if(self.isMaster()):
            pos, prob, rstate = super(MpiCosmoHammerSampler, self).loadBurnin()
        else:
            pos, prob, rstate = []
            
        pos = self.mpiBCast(pos)
        prob = self.mpiBCast(prob)
        rstate = self.mpiBCast(rstate)
        
        self.log("loading done")
        return pos, prob, rstate
    
    def createEmceeSampler(self, callable):
        """
        Factory method to create the emcee sampler
        """
        self.log("Using emcee "+str(emcee.__version__))
        #create a tuple to emulate to pool's map function using our self.mpiParallelizedMap
        pool = namedtuple('pool',['map'])(self.mpiParallelizedMap)
        
        return emcee.EnsembleSampler(self.nwalkers, self.paramCount, callable, 
                                     threads=self.threadCount, pool=pool)

    def createInitPos(self):
        """
        Factory method to create initial positions
        """   
        #bcast the positions to ensure that all mpi nodes start at the same position
        return self.mpiBCast(super(MpiCosmoHammerSampler, self).createInitPos())

    #MPI sync routines
    def mpiBCast(self, value):
        """
        Mpi bcasts the value and Returns the value from the master (rank = 0).
        """
        return MPI.COMM_WORLD.bcast(value)


    def mpiParallelizedMap(self, function,list):
        """
        Emulates a pool map function using Mpi.
        Retrieves the number of mpi processes and splits the list of walker position 
        in order to allow each process its block
        """
        (rank,size) = (MPI.COMM_WORLD.Get_rank(),MPI.COMM_WORLD.Get_size())
        #sync
        list = self.mpiBCast(list)
        #split, process and merge the list
        return self.mergeList(MPI.COMM_WORLD.allgather(self.M(function, self.splitList(list,size)[rank])))


    def splitList(self, list, n):
        """
        Splits the list into block of eqals sizes (listlength/n)
        """
        blockLen = len(list) / float(n)
        return [list[int(round(blockLen * i)): int(round(blockLen * (i + 1)))] for i in range(n)]    
    
    def mergeList(self, lists):
        """
        Merges the lists into one single list
        """
        return list(itertools.chain(*lists))
    
    def isMaster(self):
        """
        Returns true if the rank is 0
        """
        return (self._rank==0)

    
    