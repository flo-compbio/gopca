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

import logging
from copy import deepcopy
import cPickle as pickle

import numpy as np

from gopca import GOPCAInput, GOPCASignature

logger = logging.getLogger(__name__)

class GOPCAOutput(object):
    """Class representing a complete set of GO-PCA output data.

    Parameters
    ----------
    input_: `go_pca.GOPCAInput`
        The GO-PCA input data.
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

    def __init__(self,input_,genes,samples,W,Y,signatures,S):

        # W = PCA loading matrix
        # Y = PCA score matrix
        # S = GO-PCA signature matrix

        # get logger
        self._logger = logging.getLogger(__name__)

        # make sure input is valid
        if not input_.valid:
            input_.validate()

        if input_.hash is None:
            input_.calculate_hash()

        # checks
        assert isinstance(input_,GOPCAInput)
        assert isinstance(genes,list) or isinstance(genes,tuple)
        assert isinstance(samples,list) or isinstance(samples,tuple)
        assert isinstance(W,np.ndarray)
        assert isinstance(Y,np.ndarray)
        assert isinstance(signatures,list) or isinstance(signatures,tuple)
        for s in signatures:
            assert isinstance(s,GOPCASignature)
        assert isinstance(S,np.ndarray)

        assert W.shape[0] == len(genes)
        assert Y.shape[0] == len(samples)
        assert W.shape[1] == Y.shape[1]
        assert S.shape[0] == len(signatures)

        # the following does not have to hold true
        #assert W.shape[1] == input_.n_components
        #assert Y.shape[1] == input_.n_components
        #assert S.shape[1] == len(samples)

        # initialization
        self.input = deepcopy(input_)
        self.genes = tuple(genes)
        self.samples = tuple(samples)
        self.W = W.copy()
        self.Y = Y.copy()
        self.signatures = tuple(signatures)
        self.S = S.copy()

    ### magic functions
    def __repr__(self):
        param_str = '%d genes; %d samples; %d PCs; %d signatures' \
                %(self.p, self.n, self.D, self.q)
        return '<GOPCAOutput object (%s); hash = %s>' \
                %(param_str, self.__get_hash())

    def __str__(self):
        return '<GOPCAOutput object (%d signatures, %d samples)>' \
                %(self.q, self.n)

    def __eq__(self,other):
        if type(self) is not type(other):
            return False
        elif self.__get_hash() == other.__get_hash():
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.__get_hash())

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
    def __get_hash(self):
        """Calculates a MD5 hash based on GO-PCA input and output data.

        Parameters
        ----------
        None

        Returns
        -------
        str
            MD5 hash as a hex string.
        """
        hashes = [self.input.hash, hash(self.genes), hash(self.samples)]
        hashes.extend([self._hash_ndarray(a) for a in
                [self.W, self.Y, self.S]])
        hashes.extend([hash(sig) for sig in self.signatures])

        hash_str = ','.join([str(h) for h in hashes])
        h = haslib.md5(hash_str).hexdigest()
        return h
    ### end private members
  
    ### protected members
    def _hash_ndarray(self,a):
        """Calculate a hash for a NumPy array."""
        before = a.flags.writable
        a.flags.writable = False
        h = hash(a.data)
        a.flags.writable = before
        return a

    # logging convenience functions
    def _debug(self,s,*args):
        self._logger.debug(s,*args)

    def _info(self,s,*args):
        self._logger.info(s,*args)

    def _warning(self,s,*args):
        self._logger.warning(s,*args)

    def _error(self,s,*args):
        self._logger.error(s,*args)


    ### public members
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

    def get_hash(self):
        """ Calculates a MD5 hash based on GO-PCA input and output data.
        
        See documentation for `__get_hash`.
        """
        return self.__get_hash()

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
        self._info('Saving GO-PCA output to file "%s"...', path)
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
        self._info('Reading GO-PCA output from file "%s"...', path)
        output = None
        with open(path,'rb') as fh:
            output = pickle.load(fh)
        return output
