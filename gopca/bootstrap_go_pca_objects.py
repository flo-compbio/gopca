# Copyright (c) 2015 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import cPickle as pickle
import time
from copy import deepcopy

import numpy as np

from gopca.go_pca_objects import GOPCAArgumentParser,GOPCAConfig,GOPCA,GOPCAResult

class BootstrapGOPCAArgumentParser(GOPCAArgumentParser):

	def __init__(self,*args,**kwargs):

		if 'description' not in kwargs or kwargs['description'] is None:
			kwargs['description'] = 'Bootstrap-GO-PCA'

		GOPCAArgumentParser.__init__(self,*args,**kwargs)

		parser = self
		parser.add_argument('-T','--repeats',type=int,default=15)
		parser.add_argument('-F','--resample-size',type=float,default=100.0,help='in percent of original sample size')

class BootstrapGOPCAConfig(GOPCAConfig):

	valid_params = set(['repeats','resample_size','seed','pc_permutations','pc_zscore_thresh'])
	all_valid_params = GOPCAConfig.valid_params | valid_params

	def __init__(self,logger,**kwargs):

		gopca_kwargs = dict([p,v] for p,v in kwargs.iteritems() if p in GOPCAConfig.valid_params)
		GOPCAConfig.__init__(self,logger,**gopca_kwargs)

		supplied_params = set(kwargs.keys())
		unknown_params = supplied_params - GOPCAConfig.valid_params - BootstrapGOPCAConfig.valid_params
		for param in sorted(unknown_params):
			logger.warning('Bootstrap GO-PCA parameter "%s" is unknown and will be ignored.' %(param))

		# require that all parameters are provided
		for k in list(BootstrapGOPCAConfig.valid_params):
			assert k in supplied_params

		# make sure the parameters are valid
		assert isinstance(kwargs['repeats'],int) and kwargs['repeats'] >= 2
		assert isinstance(kwargs['resample_size'],float) and (0 < kwargs['resample_size'] <= 100.0)

		kwargs = dict([k,kwargs[k]] for k in BootstrapGOPCAConfig.valid_params)
		self.__dict__.update(kwargs)

	@classmethod
	def read_config_file(cls):
		raise NotImplemented

	def write_config_file(self,output_file):
		raise NotImplemented

	def __repr__(self):
		return '<BootstrapGOPCAConfig object (%s)>' %('; '.join(['%s=%s' %(k,getattr(self,str(k))) for k in sorted(BootstrapGOPCAConfig.all_valid_params)]))

	def __str__(self):
		return '<BootstrapGOPCAConfig object with parameters: %s>' %(', '.join(['%s=%s' %(k,getattr(self,str(k))) for k in sorted(BootstrapGOPCAConfig.all_valid_params)]))

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self):
		if type(self) is not type(other):
			return False
		if repr(self) == repr(other):
			return True
		else:
			return False


class BootstrapGOPCA(GOPCA):

	def __init__(self,logger,config):

		GOPCA.__init__(self,logger,config)

		assert isinstance(config,BootstrapGOPCAConfig)

		self.n_components_copy = config.n_components
		self.samples_full = None
		self.E_full = None
		self.results = None

	def read_expression(self,expression_file):
		"""
		This method first calls the GOPCA function of the same name,
		but then also makes a separate copy of the list of samples
		and the expression matrix. This is important for later
		during bootstrapping, because the native copies of these
		variables get modified.
		"""

		super(BootstrapGOPCA,self).read_expression(expression_file)

		self.samples_full = tuple(self.samples)
		self.E_full = self.E.copy()

	def filter_genes_by_variance(self,n_top):

		super(BootstrapGOPCA,self).filter_genes_by_variance(n_top)
		self.E_full = self.E.copy()

	@property
	def n_full(self):
		return len(self.samples_full)

	@property
	def resample_size(self):
		return self.config.resample_size

	@property
	def repeats(self):
		return self.config.repeats

	def run(self):

		seed_before = self.seed # keep a copy of this
		t0 = time.time()
		np.random.seed(self.seed) # initialize random number generator with seed
		gopca_results = []
		n_resample = int(self.n_full * (self.resample_size/100.0))
		self.message('Resample size = %d' %(n_resample))
		logger = self.logger
		I = np.zeros((self.repeats,n_resample),dtype=np.int64)
		for i in range(self.repeats):
			self.message('Performing repeat %d / %d...' %(i+1,self.repeats))

			# resample
			sel = np.random.choice(self.n_full,size=n_resample,replace=True)
			self.config.seed = np.random.randint(int(1e9)) # overwrite seed (used in estimating number of PCs)
			I[i,:] = sel
			self.samples = tuple([self.samples_full[j] for j in sel])
			self.E = self.E_full[:,sel]
			assert self.E.shape[1] == n_resample
			assert len(self.samples) == n_resample

			# re-estimate number of PCs
			if self.n_components_copy == 0:
				self.estimate_n_components(quiet=True)
				self.message('Estimated %d principal components.' %(self.n_components))

			# run GO-PCA silently (only report errors)
			verbosity_before = logger.verbosity
			logger.verbosity = 1

			gopca_results.append(super(BootstrapGOPCA,self).run())

			# restore original verbosity
			logger.verbosity = verbosity_before

			# output information about result
			self.message('Generated %d signatures.' %(gopca_results[-1].q))

		self.seed = seed_before # restore seed
		result = BootstrapGOPCAResult(deepcopy(self.config),I,gopca_results)
		t1 = time.time()
		self.message('Runtime = %.1fs' %(t1-t0))
		return result

class BootstrapGOPCAResult(object):
	def __init__(self,config,I,gopca_results):

		# checks
		assert isinstance(config,GOPCAConfig)
		for r in gopca_results:
			assert isinstance(r,GOPCAResult)

		self.config = config
		self.I = I
		self.gopca_results = tuple(gopca_results)

	@property
	def T(self):
		return len(self.gopca_results)

	@property
	def resample_size(self):
		return self.config.resample_size

	def save(self,output_file):
		#self.message('Saving GO-PCA result to file "%s"...' %(output_file),endline=False)
		with open(output_file,'wb') as ofh:
			pickle.dump(self,ofh,pickle.HIGHEST_PROTOCOL)
		#self.message("done!")

	@staticmethod
	def load(self,file_name):
		result = None
		with open(file_name,'rb') as fh:
			result = pickle.load(fh)
		return result

	def __repr__(self):
		conf_hash = hash(self.config)
		results_hash = hash(self.results)
		return '<BootstrapGOPCAResult object (%d results; config hash: %d; results hash: %d)' \
				%(self.T,conf_hash,results_hash)

	def __str__(self):
		conf = self.config
		return '<BoostrapGOPCAResult object with %d results, using resampling size of %.1f%%>' \
				%(self.T,self.resample_size)

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self,other):
		if type(self) is not type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False
