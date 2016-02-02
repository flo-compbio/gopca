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

"""Module containing the `GOPCAResult` class.
"""

import logging
import hashlib
from copy import deepcopy
import cPickle as pickle

import numpy as np

from gopca import GOPCAConfig, GOPCASignature

logger = logging.getLogger(__name__)

class GOPCAResult(object):
    """Class representing a GO-PCA result.

    Parameters
    ----------
    config: `GOPCAConfig`
        The GO-PCA configuration used.
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

    def __init__(self, config, genes, samples, W, Y, signatures, S):
        # W = PCA loading matrix
        # Y = PCA score matrix
        # S = GO-PCA signature matrix

        assert isinstance(config, GOPCAConfig)
        assert isinstance(genes, (list, tuple))
        assert isinstance(samples, (list, tuple))
        assert isinstance(W, np.ndarray)
        assert isinstance(Y, np.ndarray)
        assert isinstance(signatures, (list, tuple))
        for s in signatures:
            assert isinstance(s, GOPCASignature)
        assert isinstance(S, np.ndarray)

        assert W.shape[0] == len(genes)
        assert Y.shape[0] == len(samples)
        assert W.shape[1] == Y.shape[1]
        assert S.shape[0] == len(signatures)

        self.config = deepcopy(config)

        self.genes = None
        if genes is not None:
            self.genes = tuple(genes)
        self.samples = tuple(samples)
        self.W = None
        if W is not None:
            self.W = W.copy()
            self.W.flags.writeable = False
        self.Y = None
        if Y is not None:
            self.Y = Y.copy()
            self.Y.flags.writeable = False

        self.signatures = tuple(signatures)
        self.S = S.copy()
        self.S.flags.writeable = False

    ### magic functions
    def __repr__(self):
        return '<GOPCAResult: %d signatures (hash=%d)>' %(self.q, hash(self))

    def __str__(self):
        param_str = 'data: %d genes, %d samples, %d PCs tested' \
                %(self.p, self.n, self.D)

        return '<GOPCAReult with %d signatures (%s)>' \
                %(self.q, param_str)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        # internal hash function
        data = []
        data.append(hash(self.signatures))
        data.append(hash(self.samples))
        data.append(hash(self.S.data))
        data.append(hash(self.genes))
        data.append(hash(self.W.data))
        data.append(hash(self.Y.data))
        return hash(tuple(data))

    def __setstate__(self, d):
        # called upon unpickling
        self.__dict__ = d
        # set "writeable" flag to False for all ndarrays, making them hashable
        # (value of writeable flag is not stored in the pickle)
        if self.W is not None:
            self.W.flags.writeable = False
        if self.Y is not None:
            self.Y.flags.writeable = False
        self.S.flags.writeable = False

    def get_hash(self):
        """Calculate MD5 hash value for the result.

        The configuration data is excluded from the hash value calculation.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The MD5 hash value (as hex string).
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

    def write_pickle(self, path):
        """Save the current object to a pickle file.

        Parameters
        ----------
        path: str
            The path of the pickle file.

        Returns
        -------
        None
        """
        logger.info('Writing GO-PCA result to pickle file "%s"...', path)
        with open(path, 'wb') as ofh:
            pickle.dump(self, ofh, pickle.HIGHEST_PROTOCOL)
