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

"""Module containing the `GOPCASignature` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import re
import logging
import copy
from collections import OrderedDict
import hashlib

import pandas as pd
import numpy as np

from genometools.expression import ExpMatrix, ExpProfile
from genometools.expression.visualize import ExpHeatmap
from genometools.expression import cluster
from genometools.enrichment import RankBasedGSEResult

# from .signature_matrix import GOPCASignatureMatrix

logger = logging.getLogger(__name__)


class GOPCASignature(object):
    """A GO-PCA signature.

    The goal of the GO-PCA algorithm is to define gene "signatures" that are
    likely to represent biologically relevant similarities and differences
    among samples.

    A GO-PCA signature consists of a set genes and their expression profiles.
    Genes in a signature are related to each other in two ways:

    1. All signature genes are members of a specific gene set.
       Gene sets are supplied to GO-PCA by the user and correspond to groups of
       genes that are known to be related to each other in some way. For
       example, when GO-PCA is run with gene sets derived from Gene Ontology
       (GO) annotations, all genes in a gene set are known to be annotated
       with the same GO term, indicating a functional relationship among them.

    2. The genes have been found to be strongly correlated with each other
       in the sense that they all contribute strongly to the same principal
       component (PC) of the expression matrix.

    Parameters
    ----------
    pc: int
        See :attr:`pc` attribute.
    gse_result: `genometools.enrichment.RankBasedGSEResult`
        See :attr:`gse_result` attribute.
    matrix: `genometools.expression.ExpMatrix`
        See :attr:`matrix` attribute.

    Attributes
    ----------
    pc: int
        The principal component (PC) that the signature was derived from
        (starting at 1), with the sign of the integer indicating the way in
        which genes were ranked based on their PC loadings. If the sign is
        positive, then the signature was derived based on an ascending order.
        Conversely, if the sign is negative, then the signature was dervied
        based on a descending ranking.
    gse_result: `RankBasedGSEResult`
        The result of the XL-mHG test that was conducted after ranking the
        genes based on their principal component loadings.
    matrix: `genometools.expression.ExpMatrix`
        Gene-by-sample matrix containing the original expression values of
        all signature genes.
    """
    _abbrev = [('positive ', 'pos. '), ('negative ', 'neg. '),
               ('interferon-', 'IFN-'), ('proliferation', 'prolif.'),
               ('signaling', 'signal.')]
    """Abbreviations used in generating signature labels."""

    def __init__(self, pc, gse_result, matrix):

        assert isinstance(pc, int)
        assert isinstance(gse_result, RankBasedGSEResult)
        assert isinstance(matrix, ExpMatrix)

        self._pc = pc
        self._gse_result = gse_result
        self._matrix = matrix

    def __repr__(self):
        return '<%s instance (pc=%d, k=%d; pval=%.1e; hash="%s")>' \
                % (self.__class__.__name__,
                   self._pc, self.k, self.pval, self.hash)

    def __str__(self):
        return '<%s "%s" (pval=%.1e, escore=%.1f)>' \
                % (self.__class__.__name__,
                   self.label, self.pval, self.escore)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return repr(self) == repr(other)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not (self == other)

    @property
    def _data(self):
        data_str = ';'.join([
            str(repr(var)) for var in
            [self._pc, self._gse_result, self._matrix]
        ])
        data = data_str.encode('UTF-8')
        return data

    def __hash__(self):
        return hash(self._data)

    @property
    def pc(self):
        return self._pc

    @property
    def gse_result(self):
        return copy.deepcopy(self._gse_result)

    @property
    def matrix(self):
        return copy.deepcopy(self._matrix)

    @property
    def hash(self):
        return str(hashlib.md5(self._data).hexdigest())

    @property
    def k(self):
        """ The number of genes in the signature. """
        return self._matrix.p

    @property
    def n(self):
        """ The number of samples. """
        return self._matrix.n

    @property
    def genes(self):
        """The genes in the signature."""
        return self._matrix.genes

    @property
    def samples(self):
        """The sample labels."""
        return self._matrix.samples

    @property
    def X(self):
        """The expression array."""
        return self._matrix.X

    @property
    def expression(self):
        return self.get_expression()

    @property
    def s(self):
        """The signature expression vector."""
        return self.expression.values

    @property
    def gene_set(self):
        """ The gene set that the signature is based on."""
        return self._gse_result.gene_set

    @property
    def gene_set_id(self):
        """ The ID of the gene set that the signature is based on."""
        return self._gse_result.gene_set.id

    @property
    def pval(self):
        """The p-value of the enrichment test."""
        return self._gse_result.pval

    @property
    def escore(self):
        """The E-escore of the enrichment test."""
        return self._gse_result.escore


    @property
    def mHG_cutoff(self):
        """ The cutoff at which the XL-mHG test statistic was attained. """
        return self._gse_result.cutoff

    @property
    def mHG_k(self):
        """ The number of genes within the gene set above the mHG cutoff. """
        return self._gse_result.k

    @property
    def mHG_K(self):
        """ The total number of genes in the gene set. """
        return self._gse_result.K

    @property
    def mHG_N(self):
        """ The total number of genes in the analysis. """
        return self._gse_result.N

    @property
    def escore_pval_thresh(self):
        return self._gse_result.escore_pval_thresh

    @property
    def label(self):
        return self.get_label(include_id=False)

    @property
    def median_correlation(self):
        C = np.corrcoef(self.X)
        ind = np.triu_indices(self.k, k=1)
        return np.median(C[ind])

    def get_expression(self, standardize=False, center=True, use_median=True):
        #label = self.get_label(max_name_length=max_name_length,
        #                       include_id=include_id)

        matrix = self.matrix # creates a copy
        if standardize:
            matrix.standardize_genes(inplace=True)
        elif center:
            matrix.center_genes(use_median=use_median, inplace=True)

        return ExpProfile(matrix.mean(axis=0), label=self)

    def get_ordered_dict(self):
        elements = OrderedDict([
            ['label', ['Label', r'%s']],
            ['pc', ['PC', r'%d']],
            ['gene_set_id', ['Gene set ID', r'%s']],
            ['k', ['k', r'%d']],
            ['K', ['K', r'%d']],
            ['pval', ['P-value', r'%.1e']],
            ['escore', ['E-score (psi=%.1e)' % self.escore_pval_thresh,
                        r'%.1f']],
            ['median_correlation', ['Median correlation', r'%.2f']],
            ['genes', ['Genes', r'%s']]
        ])
        od = OrderedDict([v[0], v[1] % getattr(self, k)]
                         for k, v in elements.items())
        return od

    def get_label(self, max_name_length=0, include_stats=True,
                  include_id=True, include_pval=False,
                  include_coll=True):
        """Generate a signature label."""

        assert isinstance(max_name_length, int)

        gene_set = self.gene_set
        name = gene_set.name

        # make sure name does not exceed max. length
        for abb in self._abbrev:
            name = re.sub(abb[0], abb[1], name)
        if 0 < max_name_length < len(name):
            name = name[:(max_name_length - 3)] + '...'

        label = name
        if include_coll and gene_set.collection is not None:
            label = '%s: %s' % (gene_set.collection, label)

        if include_id:
            label = label + ' (%s)' % gene_set.id

        stats_str = ''
        if include_stats:
            e_str = ''
            p_str = ''
            if include_pval:
                p_str = ', p=%.1e' % self.pval
                if self.escore is not None:
                    e_str = ', e=%.1f' % self.escore
            stats_str = ' [%d:%d/%d%s%s]' \
                        % (self._pc, self.k, self.mHG_K, e_str, p_str)

        label = label + stats_str
        return label

    def get_heatmap(
            self, sig_matrix=None,
            standardize=False, center=True, use_median=True,
            include_id=False,
            cluster_genes=True,
            gene_cluster_metric='correlation',
            cluster_samples=True,
            sample_cluster_metric='euclidean',
            cluster_method='average',
            sig_matrix_kw=None,
            **kwargs):

        """Generate a plotly heatmap showing the expr. of signature genes."""
        # TODO: Finish docstring

        if sig_matrix_kw is None:
            sig_matrix_kw = {}

        assert isinstance(cluster_genes, bool)
        assert isinstance(cluster_samples, bool)
        assert isinstance(gene_cluster_metric, str)
        assert isinstance(sample_cluster_metric, str)
        assert isinstance(cluster_method, str)
        assert isinstance(sig_matrix_kw, dict)

        # we cannot simply import GOPCASignatureMatrix because it creates a
        # circular dependency (how to fix this?)
        # if sig_matrix is not None:
        #    assert isinstance(sig_matrix, GOPCASignatureMatrix)

        matrix = self.matrix  # this creates a copy
        if standardize:
            matrix.standardize_genes(inplace=True)
        elif center:
            matrix.center_genes(use_median=use_median, inplace=True)

        # clustering
        if sig_matrix is not None:
            # user has provided a GOPCASignatureMatrix instance
            # make sure its samples match the signature's
            logger.info('Ordering samples to match order in signature matrix.')
            assert set(sig_matrix.samples) == set(self.samples.values)

            # get signature matrix expression values
            sig_matrix_expr = sig_matrix.get_expression_matrix(**sig_matrix_kw)

            # re-arrange samples according to clustering of signature matrix
            matrix = matrix.loc[:, sig_matrix_expr.samples]

        elif cluster_samples:
            # cluster samples
            matrix = cluster.cluster_samples(
                matrix, metric=sample_cluster_metric, method=cluster_method
            )

        # clustering
        if cluster_genes:
            # cluster signatures
            matrix = cluster.cluster_genes(
                matrix, metric=gene_cluster_metric, method=cluster_method
            )

        cb_label = kwargs.pop('colorbar_label', 'Expression')

        # add a "Signature"-labeled row to the top,
        # which represents the signature expression vector
        title = self.get_label(include_id=include_id)
        mean = np.mean(matrix.X, axis=0)
        header_row = ExpMatrix(genes=['<b>Signature</b>'],
                               samples=matrix.samples,
                               X=np.atleast_2d(mean))
        combined_matrix = pd.concat([header_row, matrix], axis=0)

        heatmap = ExpHeatmap(combined_matrix)

        return heatmap


