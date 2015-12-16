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

"""Module containing the GOPCAConfig class.

"""

import os
import hashlib
import logging
#from ConfigParser import SafeConfigParser

import numpy as np

logger = logging.getLogger(__name__)

class GOPCAConfig(object):
    """GO-PCA configuration class.

    GO-PCA configuration data consists of GO-PCA parameter values and the paths
    (names) of the input and output files. (Internally, the paths of the input
    and output files are also treated as parameters.)

    GO-PCA parameters can be specified upon class instantiation, or at a later
    time --- using the `set_param` and `set_params` functions. Parameters are
    exposed as virtual attributes of the class, using the `__getattr__` magic.

    The class also supports reading and writing of INI-style configuration
    files (see `read_config_file` and `write_config_file`).
    [Not yet implemented!]

    Parameters
    ----------
    params: dict, optional
        Dictionary containing GO-PCA parameter values.

    Attributes
    ----------
    hash: int (property)
        MD5 hash value uniquely identifying the current configuration. The
        paths of the input and output files are excluded from the hash value
        calculation.
    """

    ### static members
    param_defaults = {
        'sel_var_genes': 0, # do not apply variance filter
        'n_components': 0, # determine # PCs automatically
        'pval_thresh': 1e-6,
        'sig_corr_thresh': 0.5,
        'mHG_X_frac': 0.25,
        'mHG_X_min': 5,
        'mHG_L': 0, # will be set to int(0.125*p), where p is the number of genes
        'escore_pval_thresh': 1e-4,
        'escore_thresh': 2.0,
        'no_local_filter': False,
        'no_global_filter': False,
        'pc_seed': 0,
        'pc_permutations': 15, 
        'pc_zscore_thresh': 2.0,
        'go_part_of_cc_only': False,
    }
    """GO-PCA parameter default values."""

    input_file_param_names = set([
        'expression_file',
        'go_annotation_file',
        'ontology_file',
    ])
    """Names of all GO-PCA input file parameters."""

    file_param_names = input_file_param_names | set(['output_file'])
    """Names of all GO-PCA file parameters."""

    param_names = set(param_defaults.keys()) | file_param_names
    """Names of all GO-PCA parameters."""
    ### end static members

    def __init__(self, params = {}):
        self.__params = {}
        # first, set all parameters to their default values
        self.reset_params()
        # then, set specified parameters to the given values
        self.set_params(params)

    def __getattr__(self, name):
        """Redirect lookup to `__params` if ``name`` is a GO-PCA parameter.

        Note: This function is only called for non-existing attributes.
        """
        if name in GOPCAConfig.param_names:
            return self.__params[name]
        else:
            raise AttributeError('There is no GO-PCA parameter called "%s"!' \
                    %(name))

    def __deepcopy__(self, memo):
        cp = GOPCAConfig()
        cp.set_params(self.__params)
        return cp

    ### public members  
    def has_param(self, name):
        """Tests if a GO-PCA parameter exists.

        Parameters
        ----------
        name: str
            The name of the parameter.

        Returns
        -------
        bool
            Whether or not the parameter exists.
        """
        return name in GOPCAConfig.param_names

    def get_param(self, name):
        """Returns the value of a GO-PCA parameter.

        Parameters
        ----------
        name: str
            The name of the parameter.
        
        Returns
        -------
        ?
            The parameter value.
        """
        return self.__params[name]

    def get_param_strings(self):
        d = []
        for k in sorted(self.__params.keys()):
            d.append('%s: %s' %(k,str(self.__params[k])))
        return d

    def get_dict(self):
        """Returns the GO-PCA configuration as a dictionary.

        Parameters
        ----------
        None

        Returns
        -------
        dict
            The configuration.
        """
        return self.__params.copy()

    def set_param(self, name, value):
        """Set a GO-PCA parameter.

        Parameters
        ----------
        name: str
            The parameter name.
        value: ?
            The parameter value.
        """
        if name not in GOPCAConfig.param_names:
            raise ValueError('No GO-PCA parameter named "%s"!' %(param))
        self.__params[name] = value

    def set_params(self, params):
        """Sets multiple GO-PCA parameters using a dictionary.

        Parameters
        ----------
        params: dict
            Dictionary containing the parameter values.

        Returns
        -------
        None
        """
        for k,v in params.iteritems():
            self.set_param(k,v)

    def reset_params(self):
        """Reset all parameters to their default values."""
        self.params = dict([p, None] for p in GOPCAConfig.param_names)
        self.set_params(GOPCAConfig.param_defaults)

    def check(self, test_if_input_exists = True,
                test_if_output_writable = True):
        """Check if the current configuration is valid.

        Parameters:
        -----------
        test_if_input_exists: bool
            If True (default), test if the input files exist.
        test_if_output_writable: bool
            If True (default), test if output directory is writable.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If any parameter value is invalid.
        """
        def check(b):
            if not b:
                raise ValueError('Invalid parameter value!')

        # check input parameters
        # specification of an ontology file is optional
        check(isinstance(self.expression_file,str))
        check(isinstance(self.go_annotation_file,str))
        if self.ontology_file is not None:
            check(isinstance(self.ontology_file,str))

        if test_if_input_exists:
            check(os.path.isfile(self.expression_file))
            check(os.path.isfile(self.go_annotation_file))
            if self.ontology_file is not None:
                check(os.path.isfile(self.ontology_file))
            
        check(isinstance(self.output_file,str))
        if test_if_output_writable:
            output_dir = os.path.dirname(self.output_file)
            if output_dir == '':
                output_dir = '.'
            logger.debug('Output directory: %s', output_dir)
            check(os.access(output_dir, os.W_OK))

        check(isinstance(self.sel_var_genes,int))
        check(self.sel_var_genes >= 0)

        check(isinstance(self.mHG_X_frac,(int,float)))
        check(0.0 <= float(self.mHG_X_frac) <= 1.0)

        check(isinstance(self.mHG_X_min,int))
        check(self.mHG_X_min >= 0)

        check(isinstance(self.mHG_L,int))
        check(self.mHG_L >= 0)

        check(isinstance(self.pval_thresh,(int,float)))
        check(0.0 < float(self.pval_thresh) <= 1.0)

        check(isinstance(self.escore_pval_thresh,(int,float)))
        check(0.0 < float(self.escore_pval_thresh) <= 1.0)

        check(isinstance(self.escore_thresh,(int,float)))
        check(float(self.escore_thresh) >= 0.0)

        check(isinstance(self.pc_seed,int))
        check(self.pc_seed >= 0 and self.pc_seed <= np.iinfo(np.uint32).max)

        if self.n_components == 0:
            check(isinstance(self.pc_permutations,int))
            check(self.pc_permutations > 0)

            check(isinstance(self.pc_zscore_thresh,float))

        #check(isinstance(self.go_part_of_cc_only, bool))

    @property
    def hash(self):
        data = []

        # ignore file names and contents in hash calculation
        hash_params = sorted(GOPCAConfig.param_names -
                GOPCAConfig.file_param_names)

        for p in hash_params:
           data.append(str(self.__params[p]))

        data_str = ','.join(data)
        logger.debug('Configuration data string: %s', data_str)
        return hashlib.md5(data_str).hexdigest()

    def read_config_file(self, path):
        raise NotImplemented

    def write_config_file(self, output_file):
        raise NotImplemented
