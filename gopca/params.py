# Copyright (c) 2015, 2016 Florian Wagner
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

"""Module containing the `GOPCAParams` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

# import os
import io
import hashlib
import logging
from collections import OrderedDict

from configparser import ConfigParser

import numpy as np

# from genometools import misc

logger = logging.getLogger(__name__)


class GOPCAParams(object):
    """A set of GO-PCA parameters.

    These parameters are for use in combination with a specific collection of
    gene sets. They exclude "global" GO-PCA parameters --- those that are
    independent of the gene sets used, like the number of principal components
    to test.

    Parameters
    ----------
    params: dict, optional
        Dictionary containing GO-PCA parameter values.

    Notes
    -----
    Parameter values can be specified upon class instantiation, or at a later
    time --- using the `set_param` and `set_params` functions. Parameters that
    are left unspecified are assigned default values.

    Parameters are exposed as virtual class attributes (using the
    `__getattr__` magic).

    The class also supports reading and writing of INI-style configuration
    files (see :func:`read_ini` and :func:`write_ini`).
    """
    __param_defaults = OrderedDict([
        ('pval_thresh', 1e-6),
        ('mHG_X_frac', 0.25),
        ('mHG_X_min', 5),
        ('mHG_L', -1),  # will be set to p / 8, where p is # genes
        ('escore_pval_thresh', 1e-4),
        ('escore_thresh', 2.0),
        ('no_local_filter', False),
        ('no_global_filter', False),
        ('sig_corr_thresh', 0.5),
        ('sig_min_genes', 5),
        ('go_part_of_cc_only', False),
    ])
    """Configuration-specific GO-PCA parameter default values."""

    @staticmethod
    def get_param_defaults():
        return GOPCAParams.__param_defaults.copy()

    def __init__(self, params=None):

        if params is None:
            params = {}

        assert isinstance(params, dict)

        self.__params = OrderedDict()
        # first, set all parameters to their default values
        self.reset_params()
        # then, set specified parameters to the given values
        self.set_params(params)

    def __getattr__(self, name):
        """
        Note: This function is only called for non-existing attributes.
        """
        if name == '_GOPCAParams__params':
            raise AttributeError()

        try:
            return self.__params[name]
        except KeyError:
            raise AttributeError()

    def __contains__(self, item):
        return item in self.__param_defaults

    def __getitem__(self, key):
        """Redirect lookup to `__params` if ``name`` is a GO-PCA parameter.
        """
        try:
            return self.__params[key]
        except KeyError:
            raise AttributeError('There is no GO-PCA parameter named "%s"!'
                                 % key)

    def __setitem__(self, key, value):
        if key not in self.__param_defaults:
            raise AttributeError('There is no GO-PCA parameter named "%s"!'
                                 % key)
        self.__params[key] = value

    def __repr__(self):
        return '<%s object (hash="%s")>' \
               % (self.__class__.__name__, self.hash)

    def __str__(self):
        param_str = ', '.join(self.param_strings)
        return '<%s object (%s)>' % (self.__class__.__name__, param_str)

    def __deepcopy__(self, memo):
        cp = GOPCAParams()
        cp.set_params(self.__params)
        return cp

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return self.__dict__ == other.__dict__
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    # public members
    @property
    def hash(self):
        data_str = ';'.join([str(repr(v)) for k, v in self.__params.items()])
        data = data_str.encode('UTF-8')
        return str(hashlib.md5(data).hexdigest())

    @property
    def param_names(self):
        return list(self.__params.keys())

    @property
    def params(self):
        """Returns the GO-PCA configuration as a dictionary."""
        return self.__params.copy()

    @property
    def param_defaults(self):
        """Returns the default GO-PCA configuration as a dictionary."""
        return self.__param_defaults.copy()

    def get_default(self, name):
        if name not in self.__param_defaults:
            raise AttributeError('%s has no parameter "%s".'
                                 % (self.__class__.__name__, name))
        return self.__param_defaults[name]

    @property
    def param_strings(self):
        d = []
        for k in sorted(self.__params.keys()):
            d.append('%s=%s' % (str(k), str(self.__params[k])))
        return d

    def set_param(self, name, value):
        """Set a GO-PCA parameter.

        Parameters
        ----------
        name: str
            The parameter name.
        value: ?
            The parameter value.
        """
        self[name] = value

    def set_params(self, params):
        assert isinstance(params, dict)
        for k, v in params.items():
            self[k]  = v

    def reset_params(self):
        """Reset all parameters to their default values."""
        self.set_params(self.__param_defaults)

    def check_params(self):
        """Check if the current configuration is valid.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True iff no problems were found.

        """
        passed = [True]

        def check_type(attr, types):
            # checks whether the parameter has a certain type
            val = getattr(self, attr)
            if not isinstance(val, types):
                logger.error('Parameter "%s" = %s: invalid type '
                             '(should be %s).', attr, val, str(types))
                passed[0] = False

        """
        def check_file_exists(attr):
            # check whether the specified file exists
            path = getattr(self, attr)
            if not os.path.isfile(path):
                logger.error('Parameter "%s" = %s: file does not exist.',
                        attr, path)
                passed[0] = False

        def check_file_writable(attr):
            # check whether the specified file is writable
            path = getattr(self, attr)
            if not misc.test_file_writable(path):
                logger.error('Parameter "%s" = %s: file not writable.',
                        attr, path)
                passed[0] = False
        """

        def check_range(attr, mn=None, mx=None,
                        left_open=False, right_open=False):
            # checks if a GO-PCA parameter falls within a certain numeric range

            val = getattr(self, attr)
            in_range = True

            rel_op = {True: '<', False: '<='}

            left_rel = ''
            if mn is not None:
                left_rel = '%s %s ' % (str(mn), rel_op[left_open])
                if left_open:
                    if not mn < val:
                        in_range = False
                else:
                    if not mn <= val:
                        in_range = False

            right_rel = ''
            if mx is not None:
                right_rel = ' %s %s' % (rel_op[right_open], str(mx))
                if right_open:
                    if not val < mx:
                        in_range = False
                else:
                    if not val <= mx:
                        in_range = False

            if not in_range:
                logger.error('Parameter "%s" = %s: out of range '
                             '(should be %s %s %s).',
                             attr, val, left_rel, attr, right_rel)
                passed[0] = False

        # check types and ranges of GO-PCA parameters
        check_type('mHG_X_frac', (int, float))
        check_range('mHG_X_frac', 0, 1)

        check_type('mHG_X_min', int)
        check_range('mHG_X_min', 0)

        check_type('mHG_L', int)
        check_range('mHG_L', -1)

        check_type('pval_thresh', (int, float))
        check_range('pval_thresh', 0, 1, left_open=True)

        check_type('escore_pval_thresh', (int, float))
        check_range('escore_pval_thresh', 0, 1, left_open=True)

        check_type('escore_thresh', (int, float))
        check_range('escore_thresh', 0)

        check_type('sig_corr_thresh', (int, float))
        check_range('sig_corr_thresh', 0, 1)

        check_type('sig_min_genes', int)
        check_range('sig_min_genes', 1, self.mHG_X_min)

        # check(isinstance(self.go_part_of_cc_only, bool))
        return passed[0]

    @classmethod
    def read_ini(cls, path):
        """Reads GO-PCA configuration data form an INI-style text file.

        Parameters
        ----------
        path: str
            The file path.

        Returns
        -------
        `gopca.GOPCAParams`
            The GO-PCA configuration data.
        """
        params = {}
        with io.open(path, encoding='UTF-8') as fh:
            config = ConfigParser()
            config.optionxform = lambda x: x
            config.read_file(fh)
            if 'GO-PCA' not in config:
                logger.error('Config file has no [GO-PCA] section!')
                return None
            d = config['GO-PCA']
            for p, v in d.items():
                if p in cls.__param_defaults:
                    t = type(cls.__param_defaults[p])
                    if t == bool:
                        v = d.getboolean(p)
                    elif t == float:
                        v = d.getfloat(p)
                    elif t == int:
                        v = d.getint(p)
                    params[p] = v
                else:
                    logger.warning('Ignoring parameter with unknown name '
                                   '"%s".', p)
        return cls(params)

    def write_ini(self, path):
        """Write configuration data to an INI-style text file.

        Parameters
        ----------
        path: str
            The file path

        Returns
        -------
        None
        """
        config = ConfigParser()
        config.optionxform = lambda x: x
        config['GO-PCA'] = OrderedDict()
        g = config['GO-PCA']
        for p, v in self.__params.items():
            if v is not None:
                g[p] = str(v)

        with io.open(path, mode='w', encoding='UTF-8') as ofh:
            config.write(ofh)
