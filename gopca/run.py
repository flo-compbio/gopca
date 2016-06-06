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

"""Module containing the `GOPCARun` class.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging
import hashlib
from copy import deepcopy
from collections import Iterable

import six

from . import GOPCAConfig, GOPCASignatureMatrix

if six.PY2:
    import cPickle as pickle
else:
    import pickle

logger = logging.getLogger(__name__)


class GOPCARun(object):
    """A GO-PCA run.

    A GO-PCA "run" consists of metadata (e.g., GO-PCA version,
    timestamp), configuration data, intermediate results, and the signatures
    generated.

    The run does *not* contain the input data itself (i.e., the original
    expression matrix and and the list of gene sets). This is so that the
    file size of the output file can be kept small.

    The run *does* contain all other information necessary to reproduce the
    results from the raw data, including the parameter settings used and hash
    values for the input data. Furthermore, it contains some intermediate
    results (e.g., the PCA loadings matrix and the PC scores) which allow for
    some additional analyses that can help gain insight into the data and the
    signatures generated.

    Parameters
    ----------
    version: str
        The GO-PCA version.
    timestamp: str
        The timestamp.
    user_config: `GOPCAConfig`
        The parameter settings provided by the user.
    final_config: `GOPCAConfig`
        The final parameter settings used.
    expression_hash: str
        Hash value for the expression input data.
    gene_sets_hash: str
        Hash value for the gene set input data.
    ontology_hash: str
        Hash value for the gene ontology input data.
    genes: list or tuple of str
        The genes in the analysis.
    samples: list or tuple of str
        The samples in the analysis.
    W: `numpy.ndarray` of floats
        The PCA loading matrix; shape = (len(genes) x # PCs).
        There must be a 1-to-1 correspondence between `genes` and the rows
        of `W`.
    Y: `numpy.ndarray` of floats
        The PC score matrix; shape = (len(samples) x # PCs).
        There must be a 1-to-1 correspondence between `samples` and the
        rows of `Y`.
    sig_matrix: `GOPCASignatureMatrix`
        The signature matrix generated.
    exec_time: float
        The execution time (in seconds).
    """
    def __init__(self, version, timestamp, user_config, final_config,
                 expression_hash, gene_sets_hash, ontology_hash,
                 genes, samples, W, Y, sig_matrix, exec_time):

        # type checks
        assert isinstance(version, str)
        assert isinstance(timestamp, str)
        assert isinstance(user_config, GOPCAConfig)
        assert isinstance(final_config, GOPCAConfig)
        assert isinstance(expression_hash, str)
        assert isinstance(gene_sets_hash, str)
        if ontology_hash is not None:
            assert isinstance(ontology_hash, str)
        assert isinstance(genes, Iterable)
        assert isinstance(samples, Iterable)
        assert isinstance(sig_matrix, GOPCASignatureMatrix)
        assert isinstance(exec_time, float)

        self.version = version
        self.timestamp = timestamp
        self.user_config = deepcopy(user_config)
        self.final_config = final_config

        self.expression_hash = expression_hash
        self.gene_sets_hash = gene_sets_hash
        self.ontology_hash = ontology_hash

        self.genes = list(genes)
        self.samples = list(samples)
        self.W = W
        self.Y = Y
        self.sig_matrix = sig_matrix
        self.exec_time = exec_time

        # make sure shapes match up
        assert W.shape[0] == len(self.genes)
        assert Y.shape[0] == len(self.samples)
        assert W.shape[1] == Y.shape[1]

    def __repr__(self):
        return '<GOPCARun instance (version="%s", timestamp="%s", hash="%s">' \
                % (self.version, self.timestamp, self.hash)

    def __str__(self):
        return '<GOPCARun instance with %d signatures>' % self.sig_matrix.q

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return repr(self) == repr(other)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def hash(self):
        data_str = ';'.join([
            str(repr(var)) for var in
            [self.version, self.timestamp,
             self.user_config, self.final_config,
             self.expression_hash, self.gene_sets_hash, self.ontology_hash,
             self.genes, self.samples, self.sig_matrix, self.exec_time]
        ])
        data_str += ';'
        data = data_str.encode('UTF-8') + \
            b';'.join([a.tobytes() for a in [
                self.Y, self.W
            ]])
        return str(hashlib.md5(data).hexdigest())

    @classmethod
    def write_pickle(cls, path):
        """Save the current object to a pickle file.

        Parameters
        ----------
        path: str
            The output file.

        Returns
        -------
        None
        """
        with open(path, 'wb') as ofh:
            pickle.dump(self, ofh, pickle.HIGHEST_PROTOCOL)

    def read_pickle(self, path):
        """Read a run from a pickle file.

        Parameters
        ----------
        path: str
            The pickle file.

        Returns
        -------
        None
        """
        with open(path, 'rb') as fh:
            run = pickle.load(fh)
        assert isinstance(run, cls)
        return run