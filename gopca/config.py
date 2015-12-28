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
import codecs
import logging
from collections import OrderedDict

from configparser import ConfigParser

import numpy as np

from genometools import misc

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
    param_defaults = OrderedDict([
        ('sel_var_genes', 0), # do not apply variance filter
        ('n_components', -1), # determine # PCs using permutation test
        ('pval_thresh', 1e-6),
        ('sig_corr_thresh', 0.5),
        ('mHG_X_frac', 0.25),
        ('mHG_X_min', 5),
        ('mHG_L', -1), # will be set to p / 8, where p is # genes
        ('escore_pval_thresh', 1e-4),
        ('escore_thresh', 2.0),
        ('no_local_filter', False),
        ('no_global_filter', False),
        ('pc_seed', 0),
        ('pc_permutations', 15), 
        ('pc_zscore_thresh', 2.0),
        ('pc_max', 0), # no limit on number of PCs to test
        ('go_part_of_cc_only', False),
    ])
    """GO-PCA parameter default values."""

    input_file_param_names = OrderedDict.fromkeys([
        'expression_file',
        'gene_ontology_file',
        'go_annotation_file',
    ])
    """Names of all GO-PCA input file parameters."""

    hash_param_names = OrderedDict.fromkeys([n + '_hash'
            for n in input_file_param_names])
    """Names of all GO-PCA file hash parameters."""

    output_file_param_name = 'output_file'

    file_param_names = OrderedDict.fromkeys(input_file_param_names.keys() +
            [output_file_param_name])
    """Names of all GO-PCA file parameters."""

    param_names = OrderedDict.fromkeys(misc.flatten(
            zip(input_file_param_names, hash_param_names)) +
            [output_file_param_name] + param_defaults.keys())
    """Names of all GO-PCA parameters."""
    ### end static members

    def __init__(self, params = {}):
        self.__params = OrderedDict()
        # first, set all parameters to their default values
        self.reset_params()
        # then, set specified parameters to the given values
        self.set_params(params)

    def __getattr__(self, name):
        """Redirect lookup to `__params` if ``name`` is a GO-PCA parameter.

        Note: This function is only called for non-existing attributes.
        """
        if name in self.param_names:
            return self.__params[name]
        else:
            raise AttributeError('There is no GO-PCA parameter called "%s"!' \
                    %(name))

    def __repr__(self):
        return '<GOPCAConfig object (hash=%d)>' %(hash(self))

    def __str__(self):
        param_str = '; '.join(self.get_param_strings())
        return '<GOPCAConfig object (%s)>' %(param_str)

    def __hash__(self):
        return hash(frozenset(self.__params.items()))

    def __deepcopy__(self, memo):
        cp = GOPCAConfig()
        cp.set_params(self.__params)
        return cp

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

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
        return name in self.param_names

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
            d.append('%s=%s' %(k,str(self.__params[k])))
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
        if name not in self.param_names:
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
        self.__params = OrderedDict([p, None]
                for p in self.param_names)
        self.set_params(self.param_defaults)

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
        bool
            True iff no problems were found.

        """

        passed = True

        def check_type(attr, types):
            # checks whether the parameter has a certain type
            val = getattr(self, attr)
            if not isinstance(val, types):
                logger.error('Parameter "%s" = %s: invalid type ' +
                        '(should be %s).', attr, val, str(types))
                passed = False

        def check_file_exists(attr):
            # check whether the specified file exists
            path = getattr(self, attr)
            if not os.path.isfile(path):
                logger.error('Parameter "%s" = %s: file does not exist.',
                        attr, path)
                passed = False

        def check_file_writable(attr):
            # check whether the specified file is writable
            path = getattr(self, attr)
            if not misc.test_file_writable(path):
                logger.error('Parameter "%s" = %s: file not writable.',
                        attr, path)
                passed = False

        def check_range(attr, mn = None, mx = None,
                left_open = False, right_open = False):
            # checks if a GO-PCA parameter falls within a certain numeric range

            val = getattr(self, attr)
            in_range = True

            rel_op = {True: '<', False: '<='}

            if mn is not None:
                left_rel = '%s %s ' %(str(mn), rel_op[left_open])
                if left_open:
                    if not mn < val: in_range = False
                else:
                    if not mn <= val: in_range = False

            if mx is not None:
                right_rel = ' %s %s' %(rel_op[right_open], str(mx))
                if right_open:
                    if not val < mx: in_range = False
                else:
                    if not val <= mx: in_range = False

            if not in_range:
                logger.error('Parameter "%s" = %s: out of range ' +
                        '(should be %s %s %s).',
                        attr, val, left_rel, attr, right_rel)
                passed = False

        # check if input files are strings
        # specification of gene ontology file is optional
        check_type('expression_file', (str, unicode))
        check_type('go_annotation_file', (str, unicode))
        if self.gene_ontology_file is not None:
            check_type('gene_ontology_file', (str, unicode))

        if test_if_input_exists:
            # check if input files exist
            check_file_exists('expression_file')
            check_file_exists('go_annotation_file')
            if self.gene_ontology_file is not None:
                check_file_exists('gene_ontology_file')
            
        # check if hash values are strings
        if self.expression_file_hash is not None:
            check_type('expression_file_hash', str)
        if self.go_annotation_file_hash is not None:
            check_type('go_annotation_file_hash', str)
        if self.gene_ontology_file is not None and \
                self.gene_ontology_file_hash is not None:
            check_type('gene_ontology_file_hash', str)

        if self.output_file is not None:
            # check if output file is a string
            check_type('output_file', (str, unicode))
            if test_if_output_writable:
                # check if output files are writable
                check_file_writable('output_file')

        # check types and ranges of GO-PCA parameters
        check_type('n_components', int)
        check_range('n_components', -1)

        check_type('sel_var_genes', int)
        check_range('sel_var_genes', 0)

        check_type('mHG_X_frac', (int,float))
        check_range('mHG_X_frac', 0, 1)

        check_type('mHG_X_min', int)
        check_range('mHG_X_min', 0)

        check_type('mHG_L', int)
        check_range('mHG_L', -1)

        check_type('pval_thresh', (int, float))
        check_range('pval_thresh', 0, 1, left_open = True)

        check_type('escore_pval_thresh', (int, float))
        check_range('escore_pval_thresh', 0, 1, left_open = True)

        check_type('escore_thresh', (int, float))
        check_range('escore_thresh', 0)

        if self.n_components == -1:
            check_type('pc_seed', int)
            check_range('pc_seed', -1, np.iinfo(np.uint32).max)

            check_type('pc_permutations', int)
            check_range('pc_permutations', 0, left_open = True)

            check_type('pc_zscore_thresh', (int,float))

            check_type('pc_max', int)
            check_range('pc_max', 0, np.iinfo(np.uint32).max)

        #check(isinstance(self.go_part_of_cc_only, bool))
        return passed

    @property
    def hash(self):
        """MD5 hash value for the current configuration."""
        data = []
        # ignore file names and contents in hash calculation
        hash_params = sorted(self.param_names - self.file_param_names)
        for p in hash_params:
           data.append(str(self.__params[p]))
        data_str = ','.join(data)
        logger.debug('Configuration data string: %s', data_str)
        return hashlib.md5(data_str).hexdigest()

    @classmethod
    def read_config(cls, path):
        """Reads GO-PCA configuration data form an INI-style text file.

        Parameters
        ----------
        path: str or unicode
            The file path.

        Returns
        -------
        `gopca.GOPCAConfig`
            The GO-PCA configuration data.
        """
        params = {}
        with codecs.open(path, 'rb', encoding = 'utf-8') as fh:
            config = ConfigParser()
            config.optionxform = lambda x: x
            config.read_file(fh)
            if 'GO-PCA' not in config:
                logger.error('Config file has no [GO-PCA] section!')
                return None
            d = config['GO-PCA']
            for p,v in d.iteritems():
                if p in cls.param_defaults:
                    t = type(cls.param_defaults[p])
                    if t == bool:
                        v = d.getboolean(p)
                    elif t == float:
                        v = d.getfloat(p)
                    elif t == int:
                        v = d.getint(p)
                params[p] = v
        return cls(params)

    def write_config(self, path):
        """Write configuration data to an INI-style text file.

        Parameters
        ----------
        path: str or unicode
            The file path

        Returns
        -------
        None
        """
        config = ConfigParser()
        config.optionxform = lambda x: x
        config['GO-PCA'] = OrderedDict()
        g = config['GO-PCA']
        for p,v in self.__params.iteritems():
            if v is not None:
                g[p] = unicode(v)

        with codecs.open(path, 'wb', encoding = 'utf-8') as ofh:
            config.write(ofh)
