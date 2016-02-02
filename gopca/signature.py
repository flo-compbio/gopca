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

"""Module containing the `GOPCASignature` class.

"""

import re
import logging
import copy
from collections import OrderedDict

import numpy as np

logger = logging.getLogger(__name__)

class GOPCASignature(object):
    """A GO-PCA signature.

    The goal of the GO-PCA algorithm is to define gene "signatures" that are
    likely to represent biologically relevant similarities and differences
    among samples.

    Given a gene expression matrix ``X``, a GO-PCA signature consists of a set
    of genes and their expression profiles. Genes in a signature are related to
    each other in two ways:

    First, all genes are members of a single gene set.
    Gene sets are supplied to GO-PCA by the user and correspond to groups of
    genes that are known to be related to each other in some way. For example,
    when GO-PCA is run with gene sets derived from Gene Ontology (GO)
    annotations, all genes in a gene set are known to be annotated with the
    same GO term, indicating a functional relationship among them.

    Second, the genes have been found to be strongly correlated with each other,
    in the sense that they all contribute strongly to the same principal
    component (PC) of the expression matrix.

    Parameters
    ----------
    genes: list or tuple of str
        See :attr:`genes` attribute.
    X: ndarray
        See :attr:`X` attribute.
    pc: int
        See :attr:`pc` attribute.
    enr: `enrichment.GSEResult`
        See :attr:`enr` attribute.

    Attributes
    ----------
    genes: tuple of str
        The list of genes in the signature. The ordering of the genes must
        correspond to the ordering of the rows in ``S``.
    X: `numpy.ndarray`
        A matrix containing the expression profiles of the ``genes``. Each gene
        corresponds to one row in the matrix, so ``E.shape`` should be
        ``(p,n)``, where ``p`` is the number of genes, and ``n`` is the number
        of samples.
    pc: int
        The principal component (PC) that the signature was derived from
        (starting at 1), with the sign of the integer indicating the way in
        which genes were ranked based on their PC loadings. If the sign is
        positive, then the signature was derived based on an ascending order.
        Conversely, if the sign is negative, then the signature was dervied
        based on a descending ranking.
    enr: `enrichment.GSEResult`
        The result of the XL-mHG test that was conducted after ranking the
        genes based on their principal component loadings.
    """

    _abbrev = [('positive ', 'pos. '), ('negative ', 'neg. '),
            ('interferon-', 'IFN-'), ('proliferation', 'prolif.'),
            ('signaling', 'signal.')]
    """Abbreviations used in generating signature labels."""

    def __init__(self, genes, X, pc, enr):
        self.genes = tuple(genes) # genes in the signature (!= self.enr.genes)

        self.X = X.copy()
        self.X.flags.writeable = False

        self.pc = pc
        self.enr = copy.deepcopy(enr)

    def __repr__(self):
        return '<%s "%s" (k=%d; p=%.1e; e=%.1f; hash=%d)>' \
                %(self.__class__.__name__, self.label,
                self.k, self.pval, self.escore, hash(self))

    def __str__(self):
        return '<%s "%s" (%d genes; p-value %.1e / E-score %.1fx)>' \
                %(self.__class__.__name__, self.label,
                self.k, self.pval, self.escore)

    def __hash__(self):
        data = []
        data.append(self.genes)
        data.append(self.X.data)
        data.append(self.pc)
        data.append(self.enr)
        return hash(tuple(data))
        
    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __ne__(self, other):
        return not (self == other)

    def __setstate__(self, d):
        self.__dict__ = d
        self.X.flags.writeable = False

    @property
    def gene_set(self):
        """ The gene set that the signature is based on. """
        return self.enr.gene_set

    @property
    def gene_set_id(self):
        return self.gene_set.id

    @property
    def pval(self):
        """ The enrichment p-value of the GO term that the signature is based on. """
        return self.enr.pval

    @property
    def escore(self):
        return self.enr.escore
    
    @property
    def k(self):
        """ The number of genes in the signature. """
        return len(self.genes)

    @property
    def K(self):
        """ The number of genes annotated with the GO term whose enrichment led to the generation of the signature. """
        return self.enr.K

    @property
    def n(self):
        """ The threshold used to generate the signature. """
        return self.enr.n

    @property
    def N(self):
        """ The total number of genes in the data. """
        return self.enr.N

    @property
    def label(self):
        return self.get_label(include_id = False)

    @property
    def median_correlation(self):
        C = np.corrcoef(self.X)
        ind = np.triu_indices(self.k,k=1)
        return np.median(C[ind])

    @property
    def escore_pval_thresh(self):
        return self.enr.escore_pval_thresh

    @property
    def gene_list(self):
        return ','.join(sorted(self.genes))

    def get_ordered_dict(self):
        elements = OrderedDict([
                ['label',['Label',r'%s']],
                ['pc',['PC',r'%d']],
                ['gene_set_id',['Gene set ID',r'%s']],
                ['k',['k',r'%d']],
                ['K',['K',r'%d']],
                ['pval',['P-value',r'%.1e']],
                ['escore',['E-score (psi=%.1e)' %(self.escore_pval_thresh),r'%.1f']],
                ['median_correlation',['Median correlation',r'%.2f']],
                ['gene_list',['Genes',r'%s']]
        ])
        od = OrderedDict([v[0],v[1] % (getattr(self,k))] for k,v in elements.iteritems())
        return od

    def get_label(self, max_name_length = 0, include_stats = True,
            include_id = True, include_pval = False,
            include_coll = True):
        """Generate a signature label."""
        enr = self.enr

        gene_set = self.gene_set
        name = gene_set.name

        # make sure name does not exceed max. length
        for abb in self._abbrev:
            name = re.sub(abb[0], abb[1], name)
        if max_name_length > 0 and len(name) > max_name_length:
            name = name[:(max_name_length - 3)] + '...'

        label = name
        if include_coll:
            label = '%s: %s' %(gene_set.collection, label)

        if include_id:
            label = label + ' (%s)' %(gene_set.id)

        stats_str = ''
        if include_stats:

            e_str = ''
            p_str = ''
            if include_pval:
                p_str = ', p=%.1e' %(self.pval)
                if self.escore is not None:
                    e_str = ', e=%.1f' %(self.escore)

            stats_str = ' [%d:%d/%d%s%s]' \
                    %(self.pc,self.k,self.K,e_str,p_str)

        label = label + stats_str
        return label


