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

"""Module containing the `GOPCAOutput` class.
"""

# TO-DO: implement output hash

import logging
import hashlib
from copy import deepcopy
import cPickle as pickle

import numpy as np

from gopca import GOPCAConfig, GOPCASignature

logger = logging.getLogger(__name__)

class GOPCAOutput(object):
    """Class representing a complete set of GO-PCA output data.

    Parameters
    ----------
    user_config: `gopca.GOPCAConfig`
        The configuration data provided by the user.
    config: `gopca.GOPCAConfig`
        The full GO-PCA configuration data.
    genes: tuple or list
        The list of genes (gene symbols) in the analysis.
    samples: tuple or list
        The list of samples (sample names) in the analysis.
    W: `numpy.ndarray`
        The PCA loading matrix; shape = (len(genes) x # PCs).
    Y: `numpy.ndarray`
        The PC score matrix; shape = (len(samples) x # PCs).
    signatures: list or tuple of `go_pca.GOPCASignature`
        The GO-PCA signatures.
    S: `numpy.ndarray`
        The GO-PCA signature matrix; shape = (len(signatures) x len(samples)).
    """

    def __init__(self, user_config, config,
                genes, samples,
                W, Y,
                signatures, S):

        # W = PCA loading matrix
        # Y = PCA score matrix
        # S = GO-PCA signature matrix

        # checks
        assert isinstance(user_config, GOPCAConfig)
        assert isinstance(config, GOPCAConfig)
        assert isinstance(genes, (list,tuple))
        assert isinstance(samples, (list,tuple))
        assert isinstance(W, np.ndarray)
        assert isinstance(Y, np.ndarray)
        assert isinstance(signatures, (list,tuple))
        for s in signatures:
            assert isinstance(s, GOPCASignature)
        assert isinstance(S, np.ndarray)

        assert W.shape[0] == len(genes)
        assert Y.shape[0] == len(samples)
        assert W.shape[1] == Y.shape[1]
        assert S.shape[0] == len(signatures)

        # initialization
        self.user_config = deepcopy(user_config)
        self.config = deepcopy(config)

        self.signatures = tuple(signatures)
        self.S = S.copy()
        self.S.flags.writeable = False

        self.genes = tuple(genes)
        self.samples = tuple(samples)
        self.W = W.copy()
        self.W.flags.writeable = False
        self.Y = Y.copy()
        self.Y.flags.writeable = False

    ### magic functions
    def __repr__(self):
        hash_str = 'signatures hash=%d; samples hash=%d; ' \
                %(hash(self.signatures), hash(self.samples)) + \
                \
                'S hash=%d' %(hash(self.S.data)) + \
                \
                'user_config hash=%d; config hash=%d; ' \
                %(hash(self.user_config), hash(self.config)) + \
                \
                'genes hash=%d; samples hash=%d; ' \
                %(hash(self.genes), hash(self.samples)) + \
                \
                'W hash=%d; Y hash=%d' \
                %(hash(self.W.data), hash(self.Y.data))

        return '<GOPCAOutput: %d signatures (%s)>' %(self.q, hash_str)

    def __str__(self):
        param_str = 'expression matrix: %d genes / %d samples / %d PCs tested' \
                %(self.p, self.n, self.D)

        return '<GOPCAOutput with %d signatures (%s) - output hash: %s>' \
                %(self.q, param_str, self.hash)

    def __eq__(self,other):
        if type(self) is not type(other):
            return False
        elif repr(self) == repr(other):
            return True
        else:
            return False

    def __hash__(self):
        return hash(repr(self))

    def __setstate__(self, d):
        self.__dict__ = d
        self.W.flags.writeable = False
        self.Y.flags.writeable = False
        self.S.flags.writeable = False
    ### end magic functions

    @property
    def p(self):
        """The number of genes in the analysis."""
        return len(self.genes)

    @property
    def n(self):
        """The number of samples in the analysis."""
        return len(self.samples)

    @property
    def D(self):
        """The number of principal components tested."""
        return self.W.shape[1]

    @property
    def q(self):
        """The number of signatures generated."""
        return len(self.signatures)

    def get_param(self, name):
        return getattr(self.input, name)

    @property
    def hash(self):
        """Calculates a MD5 hash value for the GO-PCA output.

        This explicitly ignores the configuration data, since users might only
        be interested in checking whether the output data is consistent.

        Parameters
        ----------
        None

        Returns
        -------
        str
            MD5 hash as a hex string.
        """
        data = []
        data.append(hash(self.signatures))
        data.append(hash(self.samples))
        data.append(hash(self.S.data))
        data.append(hash(self.genes))
        data.append(hash(self.W.data))
        data.append(hash(self.Y.data))

        hash_str = ','.join(str(d) for d in data)
        h = hashlib.md5(hash_str).hexdigest()

        return h
  
    def save(self,path):
        """Save the current object to a pickle file.

        Parameters
        ----------
        path: str
            The path of the pickle file.

        Returns
        -------
        None
        """
        logger.info('Saving GO-PCA output to file "%s"...', path)
        with open(path,'wb') as ofh:
            pickle.dump(self,ofh,pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(path):
        """Load a GOPCAOutput object from a pickle file.

        Parameters
        ----------
        path: str
            Path of pickle file.

        Returns
        -------
        GOPCAOutput
            The GOPCAOutput object.
        """
        logger.info('Reading GO-PCA output from file "%s"...', path)
        output = None
        with open(path,'rb') as fh:
            output = pickle.load(fh)
        return output
