#!/usr/bin/env python
import CosmoHammer_21CMMC.constants.Constants as c

import numpy as np
import pickle

class SampleFileUtil(object):
	"""
	Util for handling sample files
	
	:param filePrefix: the prefix to use
	:param master: True if the sampler instance is the master
	:param  reuseBurnin: True if the burn in data from a previous run should be used
	
	"""
	
	def __init__(self, filePrefix, master=True, reuseBurnin=False):
		self.filePrefix = filePrefix
		
		if(master):
			if(reuseBurnin):
				mode = "r"
			else:
				mode = "w"
			self.samplesFileBurnin = open(self.filePrefix+c.BURNIN_SUFFIX, mode)
			self.probFileBurnin = open(self.filePrefix+c.BURNIN_PROB_SUFFIX, mode)
			self.neutral_fractionFileBurnin = open(self.filePrefix+c.BURNIN_NF_SUFFIX, mode)
			
			self.samplesFile = open(self.filePrefix+c.FILE_SUFFIX, "w")
			self.probFile = open(self.filePrefix+c.PROB_SUFFIX, "w")
			self.neutral_fractionFile = open(self.filePrefix+c.NF_SUFFIX, "w")
	
	def importFromFile(self, filePath):
		values = np.loadtxt(filePath, dtype=float)
		return values

	def storeRandomState(self, filePath, randomState):
		pickle.dump(randomState, file(filePath,'w'))

	def importRandomState(self, filePath):
		state = pickle.load(file(filePath,'r'))
		return state

	def persistBurninValues(self, pos, prob, NF_values, data):
		self.persistValues(self.samplesFileBurnin, self.probFileBurnin, self.neutral_fractionFileBurnin, pos, prob, NF_values, data)
		
	def persistSamplingValues(self, pos, prob, NF_values, data):
		self.persistValues(self.samplesFile, self.probFile, self.neutral_fractionFile, pos, prob, NF_values, data)
		

	def persistValues(self, posFile, probFile, NFFile, pos, prob, NF_values, data):
		"""
		Writes the walker positions and the likelihood to the disk
		"""
		posFile.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
		posFile.write("\n")
		posFile.flush()

		NFFile.write("\n".join(["\t".join([str(q) for q in p]) for p in NF_values]))
		NFFile.write("\n")
		NFFile.flush()

		probFile.write("\n".join([str(p) for p in prob]))
		probFile.write("\n")
		probFile.flush();
		
	def __str__(self, *args, **kwargs):
		return "SampleFileUtil"
