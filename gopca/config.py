# Copyright (c) 2016 Florian Wagner
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

"""Module containing the `GOPCAConfig` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging
import hashlib
import copy

from genometools.basic import GeneSetCollection
from genometools.ontology import GeneOntology
from . import GOPCAParams

logger = logging.getLogger(__name__)


class GOPCAConfig(object):
    """A GO-PCA configuration.

    This class represents a complete set of inputs to the GO-PCA algorithm,
    excluding the expression data. Specifically, the inputs consist of
    parameter settings, represented by the `GOPCAParams` class, as well as a
    list of gene sets (optionally including Gene Ontology data; represented by
    the `genometools.ontology.GeneOntology` class).

    Furthermore, two versions of the parameter settings are stored: The
    first version (`user_params`) holds the values exactly as specified by the
    user. Some of those values can have special meaning, such as ``-1`` for the
    XL-mHG ``L`` parameter (``mHG_L``), which will be automatically converted
    to ``p/8`` (where ``p`` is the number of genes). The second version
    (simply called `params`) holds the values actually used by the algorithm.
    So in the previous example, the ``mHG_L`` parameter stored in this version
    would equal ``p/8`` (rounded down to the nearest integer).

    Parameters
    ----------
    user_params: `GOPCAParams`
        See :attr:`user_params` attribute.
    gene_sets: `genometools.basic.GeneSetCollection`
        See :attr:`gene_set_coll` attribute.
    gene_ontology: `genometools.ontology.GeneOntology`, optional
        See :attr:`gene_ontology` attribute.
    params: `GOPCAParams`, optional
        See :attr:`params` attribute.

    Attributes
    ----------
    user_params: `GOPCAParams`
        The GO-PCA parameter settings, as specified by the user.
    gene_sets: `genometools.basic.GeneSetCollection`
        The list of gene sets to be used by GO-PCA.
    gene_ontology: `genometools.ontology.GeneOntology` or None
        The Gene Ontology (only if gene sets are based on GO annotations.)
    params: `GOPCAParams` or None
        The final GO-PCA parameter settings, after resolving parameter values
        with special meanings (see above).
    """
    def __init__(self, user_params, gene_sets,
                 gene_ontology=None, params=None):

        # store configuration
        assert isinstance(user_params, GOPCAParams)
        assert isinstance(gene_sets, GeneSetCollection)

        if gene_ontology is not None:
            assert isinstance(gene_ontology, GeneOntology)

        if params is not None:
            assert isinstance(params, GOPCAParams)

        self.user_params = user_params
        self.params = params
        self.gene_sets = gene_sets
        self.gene_ontology = gene_ontology

    def __repr__(self):
        return '<%s instance (hash="%s")>' % \
               (self.__class__.__name__, self.hash)

    def __str__(self):
        return '<%s instance (hash="%s")>' % \
               (self.__class__.__name__, self.hash)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return self.__dict__ == other.__dict__
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def hash(self):
        data_str = ';'.join(
            repr(v) for v in [self.user_params, self.params,
                              self.gene_sets, self.gene_ontology])
        data = data_str.encode('UTF-8')
        return str(hashlib.md5(data).hexdigest())

    def finalize_params(self, num_genes):
        """Replace parameters set to special values with final value.

        Parameters
        ----------
        num_genes: int
            The number of genes in the analysis.
        """
        params = copy.deepcopy(self.user_params)

        # determine mHG_L, if -1 or 0
        if params.mHG_L == -1:
            # -1 = determine L automatically => p / 8
            params.set_param('mHG_L', int(num_genes / 8.0))
        elif params.mHG_L == 0:
            # 0 = "disable" effect of L => set it to the number of genes
            params.set_param('mHG_L', num_genes)
            
        self.params = params