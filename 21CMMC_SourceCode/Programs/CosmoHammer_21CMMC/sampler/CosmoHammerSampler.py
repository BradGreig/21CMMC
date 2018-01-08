#!/usr/bin/env python
import CosmoHammer_21CMMC
import CosmoHammer_21CMMC.constants.Constants as c

from CosmoHammer_21CMMC.util.SampleFileUtil import SampleFileUtil
from CosmoHammer_21CMMC.sampler.util.IterationStopCriteriaStrategy import IterationStopCriteriaStrategy
from CosmoHammer_21CMMC.sampler.util.VariousInitialConditionGenerators import UniformPosition

from CosmoHammer_21CMMC import emcee
import numpy as np
import logging
import time

class CosmoHammerSampler(object):
    """
    A complete sampler implementation taking care of correct setup, chain burn in and sampling.

    :param params: the parameter of the priors
    :param likelihoodComputationChain: the callable computation chain
    :param filePrefix: the prefix for the log and output files
    :param walkerRatio: the ratio of walkers and the count of sampled parameters
    :param burninIterations: number of iteration for burn in
    :param sampleIterations: number of iteration to sample
    :param stopCriteriaStrategy: the strategy to stop the sampling. 
        Default is None an then IterationStopCriteriaStrategy is used
    :param initPositionGenerator: the generator for the init walker position. 
        Default is None an then SampleBallPositionGenerator is used
    :param fileUtil: util used to store the results
    :param threadCount: The count of threads to be used for the computation. Default is 1
    :param reuseBurnin: Flag if the burn in should be reused. 
        If true the values will be read from the file System. Default is False

    """

    def __init__(self, params,likelihoodComputationChain, filePrefix, walkersRatio, burninIterations, 
                 sampleIterations, FiducialParams, param_legend, LowerBound_XRAY, UpperBound_XRAY, SpinTz, 
                 filethin = 1, stopCriteriaStrategy=None, initPositionGenerator=None, 
                 storageUtil=None, threadCount=1, reuseBurnin=False):
        """
        CosmoHammer sampler implementation

        """
        self.likelihoodComputationChain = likelihoodComputationChain
        self.walkersRatio = walkersRatio
        self.reuseBurnin = reuseBurnin
        self.filePrefix = filePrefix
        self.threadCount = threadCount
        self.paramCount = len(params[:,0])
        self.nwalkers = self.paramCount*walkersRatio
        self.burninIterations = burninIterations
        self.sampleIterations = sampleIterations
        self.filethin = filethin  

        self.FiducialParams = FiducialParams
        self.param_legend = param_legend

        self.LowerBound_XRAY = LowerBound_XRAY
        self.UpperBound_XRAY = UpperBound_XRAY
        self.SpinTz = SpinTz        

        self.lowerbounds = params[:,0]
        self.upperbounds = params[:,1]
        
        assert sampleIterations > 0, "CosmoHammer needs to sample for at least one iterations"
        
        # setting up the logging
        self._configureLogging(filePrefix+c.LOG_FILE_SUFFIX)
        
        self.log("Using CosmoHammer "+str(CosmoHammer_21CMMC.__version__))        

        # The sampler object
        self._sampler = self.createEmceeSampler(likelihoodComputationChain)
        if(storageUtil is None):
            storageUtil = self.createSampleFileUtil()
            
        self.storageUtil = storageUtil

        if(stopCriteriaStrategy is None):
            stopCriteriaStrategy = self.createStopCriteriaStrategy()
            
        stopCriteriaStrategy.setup(self)
        self.stopCriteriaStrategy = stopCriteriaStrategy
    
        if(initPositionGenerator is None):
            initPositionGenerator = self.createInitPositionGenerator()
            
        initPositionGenerator.setup(self)
        self.initPositionGenerator = initPositionGenerator
    
    def _configureLogging(self, filename):
        logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', 
                            filename=filename, filemode='w', level=logging.INFO)
        
        
    def createStopCriteriaStrategy(self):
        """
        Returns a new instance of a stop criteria stategy
        """
        return IterationStopCriteriaStrategy()
    
    def createSampleFileUtil(self):
        """
        Returns a new instance of a File Util
        """
        return SampleFileUtil(self.filePrefix, reuseBurnin=self.reuseBurnin)
    
    def createInitPositionGenerator(self):
        """
        Returns a new instance of a Init Position Generator
        """
        """
        This needs to me made more intelligent, rather than being hard coded!
        """

        return UniformPosition()
    
#    @profile
    def startSampling(self):
        """
        Launches the sampling
        """
        self.log(self.__str__())
        if(self.burninIterations>0):
            
            if(self.reuseBurnin):            
                pos, prob, NF_Values_old, rstate = self.loadBurnin()
                datas = [None]*len(pos)
                if 'singlez' in self.filePrefix:  
                    NF_Values = [[float(q)] for q in NF_Values_old]
                    NF_Values = np.array(NF_Values)
                else:
                    NF_Values = NF_Values_old

            else:
                pos, prob, NF_Values, rstate, datas = self.startSampleBurnin()
        else:
            pos = self.createInitPos()
            prob = None
            NF_Values = None
            rstate = None
            datas = None
        # Starting from the final position in the burn-in chain, sample for 1000
        # steps.
        self.log("start sampling after burn in")
        start = time.time()

        self.sample(pos, prob, NF_Values, rstate, datas)
        end = time.time()
        self.log("sampling done! Took: " + str(round(end-start,4))+"s")

        # Print out the mean acceptance fraction. In general, acceptance_fraction
        # has an entry for each walker
        self.log("Mean acceptance fraction:"+ str(round(np.mean(self._sampler.acceptance_fraction), 4)))

    def loadBurnin(self):
        """
        loads the burn in form the file system
        """
        self.log("reusing previous burn in")

        pos = self.storageUtil.importFromFile(self.filePrefix+c.BURNIN_SUFFIX)[-self.nwalkers:]

        prob = self.storageUtil.importFromFile(self.filePrefix+c.BURNIN_PROB_SUFFIX)[-self.nwalkers:]

        NF_Values = self.storageUtil.importFromFile(self.filePrefix+c.BURNIN_NF_SUFFIX)[-self.nwalkers:]        

        rstate= self.storageUtil.importRandomState(self.filePrefix+c.BURNIN_STATE_SUFFIX)

        self.log("loading done")
        return pos, prob, NF_Values, rstate
    
#    @profile
    def startSampleBurnin(self):
        """
        Runs the sampler for the burn in
        """
        self.log("start burn in")
        start = time.time()
        p0 = self.createInitPos()

        pos, prob, NF_Values, rstate, data = self.sampleBurnin(p0)
        end = time.time()
        self.log("burn in sampling done! Took: " + str(round(end-start,4))+"s")
        self.log("Mean acceptance fraction for burn in:" + str(round(np.mean(self._sampler.acceptance_fraction), 4)))
        
        self.resetSampler()
        
        return pos, prob, NF_Values, rstate, data
    
    
    def resetSampler(self):
        """
        Resets the emcee sampler in the master node
        """
        if self.isMaster():
            self.log("Reseting emcee sampler")
            # Reset the chain to remove the burn-in samples.
            self._sampler.reset()        
    
    
#    @profile
    def sampleBurnin(self, p0):
        """
        Run the emcee sampler for the burnin to create walker which are independent form their starting position
        """

        counter = 1
        counter_thin = 0
        for pos, prob, neutral_fractions, rstate, datas in self._sampler.sample(p0, iterations=self.burninIterations):

            if self.isMaster():
                counter_thin += 1
                if(counter_thin==self.filethin):
                    self.storageUtil.persistBurninValues(pos, prob, neutral_fractions, datas)
                    counter_thin = 0
                if(counter%10==0):
                    self.log("Iteration finished:" + str(counter))
                
            counter = counter + 1

        if self.isMaster():
            self.log("storing random state")
            self.storageUtil.storeRandomState(self.filePrefix+c.BURNIN_STATE_SUFFIX, rstate)
            
        return pos, prob, neutral_fractions, rstate, datas

#    @profile
    def sample(self, burninPos, burninProb=None, burninNF_Values=None, burninRstate=None, datas=None):
        """
        Starts the sampling process
        """
        counter = 1
        counter_thin = 0
        for pos, prob, NF_Values, _, datas in self._sampler.sample(burninPos, lnprob0=burninProb, neutral_fractions0=burninNF_Values, rstate0=burninRstate, 
                                                        blobs0=datas, iterations=self.sampleIterations):
            if self.isMaster():
                counter_thin += 1
                if(counter_thin==self.filethin):
                    self.storageUtil.persistSamplingValues(pos, prob, NF_Values, datas)
                    counter_thin = 0

                if(self.stopCriteriaStrategy.hasFinished()):
                    break
                
                if(counter%10==0):
                    self.log("Iteration finished:" + str(counter))
                
            counter = counter + 1


    def isMaster(self):
        """
        Returns True. Can be overridden for multitasking i.e. with MPI
        """
        return True

    def log(self, message):
        """
        Logs a message to the logfile
        """
        logging.info(message)
    
    
    def createEmceeSampler(self, callable):
        """
        Factory method to create the emcee sampler
        """
        self.log("Using emcee "+str(emcee.__version__))
#        print 'Create the Emcee Sampler'
        return emcee.EnsembleSampler(self.nwalkers, self.paramCount, callable, lower_bounds=self.lowerbounds, upper_bounds=self.upperbounds, threads=self.threadCount)


    def createInitPos(self):
        """
        Factory method to create initial positions
        """
        return self.initPositionGenerator.generate()


    def getChain(self):
        """
            Returns the sample chain
        """
        return self._sampler.chain()
    
    def __str__(self, *args, **kwargs):
        """
            Returns the string representation of the sampler config
        """
        desc = "Sampler: " + str(type(self))+"\n" \
                "configuration: \n" \
                "  Burnin iterations: " +str(self.burninIterations)+"\n" \
                "  Samples iterations: " +str(self.sampleIterations)+"\n" \
                "  Walkers ratio: " +str(self.walkersRatio)+"\n" \
                "  Reusing burn in: " +str(self.reuseBurnin)+"\n" \
                "  init pos generator: " +str(self.initPositionGenerator)+"\n" \
                "  stop criteria: " +str(self.stopCriteriaStrategy)+"\n" \
                "  storage util: " +str(self.storageUtil)+"\n" \
                "likelihoodComputationChain: \n" + str(self.likelihoodComputationChain) \
                +"\n"
        
        return desc
