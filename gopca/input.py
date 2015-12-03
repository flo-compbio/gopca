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

import os
import sys
import hashlib
import copy
import logging
from ConfigParser import SafeConfigParser

from genometools import misc

"""Module containing the GOPCAInput class.

"""

class GOPCAInput(object):
    """Class representing a complete set of GO-PCA input data.

    GO-PCA input data consists of parameter values and input files. To
    facilitate reproducibility of GO-PCA analyses, this class implements the
    "GO-PCA input hash", a unique identifier for each set of GO-PCA input data.

    GO-PCA parameters can be specified either upon class instantiation, or at
    a later time --- using `set_params`, `update_params`, or by direct
    assignment (e.g., :code:`inpt.mHG_X_min = 5`). The paths of the input files
    are treated as parameters (however, they are ignored in the calculation of
    the input hash; see below).

    After all parameters are specified, the `validate` function should be
    called. This function ensures that all parameter values are valid, meaning
    that they have the right data type, and that their values fall within
    acceptable ranges (e.g., the p-value threshold ``P`` has to be larger than
    0 and not larger than 1). The function also verifies that the specified
    input files exist.

    For valid input data, the `get_hash` function can be called to generate a
    hash value (unique identifier) for the data. The hash value is calculated
    by first computing the MD5 checksums for all input files, and then storing
    storing them in a tuple along with all other parameter values. The hash
    value of this tuple is used as the GO-PCA input hash.

    Any time a parameter value (or input file name) is modified, the hash value
    is reset (to ``None``), and the `get_hash` function needs to be called
    again, in order to compute the new hash. This ensures that any hash
    obtained corresponds to the current set of input data.

    The class furhtermore supports reading and writing of INI-style
    configuration files (see `read_config_file` and `write_config_file`).

    Parameters
    ----------
    params: dict, optional
        Dictionary containing GO-PCA parameter values.

    Attributes
    ----------
    validated: bool
        Whether the current set of parameters has been "validated". If this is
        True, then the ``hash`` property returns a hash value for the current
        GO-PCA configuration.
    hash: int (property)
        Hash value uniquely identifying the current configuration. ``None`` if
        current configuration has not been "validated" with `self.validate`.
    """

    ### static members
    __param_defaults = {
        'expression_file': None,
        'go_annotation_file': None,
        'ontology_file': None,
        'output_file': None,
        'sel_var_genes': 0, # do not apply variance filter
        'n_components': 0, # determine # PCs automatically
        'pval_thresh': 1e-6,
        'sig_corr_thresh': 0.5,
        'mHG_X_frac': 0.25,
        'mHG_X_min': 5,
        'mHG_L': None, # will be set to int(0.125*p), where p is the number of genes
        'escore_pval_thresh': 1e-4,
        'escore_thresh': 2.0,
        'disable_local_filter': False,
        'disable_global_filter': False,
        'seed': 0,
        'pc_permutations': 15, 
        'pc_zscore_thresh': 2.0,
        'go_part_of_cc_only': False,
    }
    """Default values for all GO-PCA parameters, including input/output files.
    """

    __input_file_params = set([
        'expression_file',
        'go_annotation_file',
        'ontology_file',
    ])
    """Names of all parameters corresponding to input files."""

    __output_file_params = set([
        'output_file',
    ])
    """Names of all parameters corresponding to output files."""

    __param_names = set(__param_defaults.keys())

    @staticmethod
    def _get_file_md5sum(path,mode='r'):
        """Get MD5 hash of file content.

        Parameters
        ----------
        path: str
            Path of file.
        
        Returns
        -------
        str
            MD5 hash of file content, represented as a 32-digit hex string.
        """
        digest = None
        with misc.smart_open(path,mode=mode,try_gzip=True) as fh:
            digest = hashlib.md5(fh.read()).hexdigest()
        return digest
    ### end static members

    ### magic functions
    def __init__(self,params = {}):

        # create logger (use parent logger)
        self._logger = logging.getLogger(__name__)

        self.__params = {}

        self.__valid = None
        self.__input_hashes = dict([p,None] for p in self.__input_file_params)
        self.__hash = None

        # set parameters (if params is empty, copy the default values)
        self.__params = copy.deepcopy(self.__param_defaults)
        self.set_params(params)

    def __getattr__(self,name):
        """Redirect lookup to `__params` if ``name`` is a GO-PCA parameter."""
        if name in GOPCAInput.__param_defaults:
            return self.__params[name]
        else:
            # __getattr__ only gets called for unknown attributes. Therefore,
            # if we're not given a GO-PCA parameter name, this is an error
            raise AttributeError('There is no GO-PCA parameter called "%s"!' \
                    %(name))

    def __setattr__(self,name,value):
        """Store value in `__params` if ``name`` is a GO-PCA parameter."""
        if name in GOPCAInput.__param_defaults:
            self.set_param(name,value)
        else:
            self.__dict__[name] = value

    def __repr__(self):
        return '<GOPCAInput object (%s)>' \
                %('; '.join(['%s=%s' %(k,getattr(self,k)) for k in
                sorted(self.__param_defaults.keys())]))

    def __str__(self):
        return '<GOPCAInput object with parameters: %s>' \
                %(', '.join(['%s=%s' %(k,getattr(self,k)) for k in
                sorted(self.__param_defaults.keys())]))

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self):
        if type(self) is not type(other):
            return False
        if repr(self) == repr(other):
            return True
        else:
            return False

    def __deepcopy__(self,memo):
        cp = GOPCAInput()
        cp.update_params(self.__params)
        cp.__valid = self.__valid
        cp.__input_hashes = copy.deepcopy(self.__input_hashes,memo)
        cp.__hash = self.__hash
        return cp

    def __getstate__(self):
        """Called to obtain the data to be pickled.

        We need to remove the logger object before pickling, since pickling
        this object would result in an error.
        """
        d = self.__dict__.copy()
        del d['_logger']
        return d

    def __setstate__(self,state):
        """Called to unpickle the object.

        We restore the logger object that was deleted before pickling.
        """
        state['_logger'] = logging.getLogger(__name__)
        self.__dict__.update(state)
    ### end magic functions

    ### private members

    def __set_param(self,name,value):
        if name not in self.__param_names:
            raise ValueError('No GO-PCA parameter named "%s"!' %(param))
        before = self.__params[name]
        self.__params[name] = value

        if before != value:
            self.__valid = False
            self.__hash = None

            if name in self.__input_file_params:
                # if name of an input file was changed, reset hash
                self.__input_hashes[name] = None

    def __validate(self):
        """Test if current input data is valid.

        The function does not inspect the contents of the input files, it only
        confirms that they exist. The function does not test whether the output
        file is writable.

        Parameters:
        -----------
        None

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If any parameter value is invalid or an input file does not exist.
        """
        # check specification and existence of expression and GO annotation
        # files
        def check(b):
            if not b: raise ValueError('Invalid parameter!')

        try:
            check(isinstance(self.expression_file,str) and
                os.path.isfile(self.expression_file))
            check(isinstance(self.go_annotation_file,str) and
                os.path.isfile(self.go_annotation_file))

            # specification of an ontology file is optional
            # but if it *is* specified (not None), verify its existence
            check(self.ontology_file is None or
                (isinstance(self.ontology_file,str) and
                os.path.isfile(self.ontology_file)))

            check(isinstance(self.sel_var_genes,int))
            check(self.sel_var_genes >= 0)

            check(isinstance(self.mHG_X_frac,(int,float)))
            check(0.0 <= float(self.mHG_X_frac) <= 1.0)

            check(isinstance(self.mHG_X_min,int))
            check(self.mHG_X_min >= 0)

            check(self.mHG_L is None or isinstance(self.mHG_L,int))
            if self.mHG_L is not None:
                check(self.mHG_L >= 0)

            check(isinstance(self.pval_thresh,(int,float)))
            check(0.0 < float(self.pval_thresh) <= 1.0)

            check(isinstance(self.escore_pval_thresh,(int,float)))
            check(0.0 < float(self.escore_pval_thresh) <= 1.0)

            check(isinstance(self.escore_thresh,(int,float)))
            check(float(self.escore_thresh) >= 0.0)

            check(isinstance(self.seed,int))
            check(self.seed >= 0)

            check(isinstance(self.pc_permutations,int))
            check(self.pc_permutations > 0)

            check(isinstance(self.pc_zscore_thresh,float))

        except ValueError:
            self.__valid = False
            raise

        else:
            self.__valid = True

    def __calculate_hash(self):
        invalid_error = 'Cannot calculate GO-PCA input hash, because input ' +\
                'data is not valid!'

        # test if input data is valid
        if self.__valid is None:
            try:
                self.__validate() # raises ValueError if input data is invalid
            except ValueError:
                self._error(invalid_error)
                raise
        
        if not self.__valid:
            self.error(invalid_error)
            raise ValueError(invalid_error)

        # input data is valid, calculate hash
        data = []
        for p in sorted(self.__input_file_params):
            # recalculate input file MD5 hashes, if necessary
            if self.__params[p] is None:
                self.__input_hashes[p] = ''
            elif self.__input_hashes[p] is None:
                fn = self.__params[p]
                md5hash = GOPCAInput._get_file_md5sum(fn)
                self.__input_hashes[p] = md5hash
                self._info('MD5 hash for file "%s": %s', p, md5hash)
            data.append(self.__input_hashes[p])

        param_keys = sorted(set(self.__params.keys()) -
            (self.__input_file_params | self.__output_file_params))

        for k in param_keys:
           data.append(str(self.__params[k]))

        data_str = ','.join(data)
        #self._logger.propagate = False
        self._debug('Input data string: %s', data_str)
        self.__hash = hashlib.md5(data_str).hexdigest()
    ### end private members

    ### logging convenience functions
    def _debug(self,*args):
        self._logger.debug(*args)

    def _info(self,*args):
        """Convenience function for reporting messages."""
        self._logger.info(*args)

    def _warning(self,*args):
        """Convenience function for reporting warnings."""
        self._logger.warning(*args)

    def _error(self,*args):
        """Convenience function for reporting errors."""
        self._logger.error(*args)
    ### end logging convenience functions

    ### public members
    @property
    def valid(self):
        return self.__valid

    @property
    def hash(self):
        return self.__hash

    @property
    def param_names(self):
        """Returns a set of all GO-PCA parameter names."""
        return self.__param_names

    def get_param_strings(self):
        d = []
        for k in sorted(self.__params.keys()):
            d.append('%s: %s' %(k,str(self.__params[k])))
        return d

    def get_default_param(self,name):
        """Return the default value of a GO-PCA parameter.

        Parameters
        ----------
        name: str
            The name of the GO-PCA parameter.

        Returns
        -------
        variable
            The parameter value.
        """
        return self.__param_defaults[name]

    def set_param(self,param,value):
        self.__set_param(param,value)

    def set_params(self,params):
        """Set GO-PCA parameters.

        This function will set unspecified parameters (i.e., parameters not
        contained in ``params``) to their default values, overwriting all
        previously stored configuration data. If this is not desired, use
        `update_params` instead.

        Parameters
        ----------
        params: dict
            Dictionary containing parameter name/value pairs.

        Returns
        -------
        None
        """

        # set unspecified parameters to their default value
        all_params = set(self.__param_defaults.keys())
        supplied_params = set(params.keys())
        unspecified_params = all_params - supplied_params
        for k in unspecified_params:
            self.__set_param(k,self.__param_defaults[k])

        # call self.update_params to set new parameters
        self.update_params(params)

    def update_params(self,params):
        for k,v in params.iteritems():
            self.__set_param(k,v)

    def validate(self):
        """Calls `__validate`."""
        self.__validate()

    def calculate_hash(self):
        """Calls `__calculate_hash`."""
        self.__calculate_hash()

    def read_config_file(self,path):
        raise NotImplemented

    def write_config_file(self,output_file):
        raise NotImplemented
    ### end public members
