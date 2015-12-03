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

import re
from collections import OrderedDict

import numpy as np

"""Module containing the `GOPCASignature` class.

"""

class GOPCASignature(object):
    """Class representing a GO-PCA signature.

    The goal of the GO-PCA algorithm is to define gene "signatures" that are
    likely to represent biologically relevant similarities and differences
    among the samples.

    Given an gene expression matrix ``E``, a GO-PCA signature consists of a set
    of genes and their expression profiles. The genes in a signature have been
    found by GO-PCA to be related to each other in two ways:
    
    First, all genes in the signature were found to be annotated with the same
    GO term. In other words, there existed at least one GO term ``T`` such that
    all genes in the signature were annotated with T. (The set of genes that is
    considered annotated with a particular term can be determined based on
    data from the `UniProt-GOA`__ database.)
    
    __ uniprot_goa_
    
    Secondly, the genes in the signature were found to all contribute to the
    same principal component (PC) of ``E``, albeit to various extents. As a
    consequence, their expression profiles exhibit a certain degree of
    correlation.

    .. _uniprot_goa: http://www.ebi.ac.uk/GOA

    Parameters
    ----------
    genes: list or tuple of str
        See :attr:`genes` attribute.
    E: ndarray of floats
        See :attr:`E` attribute.
    pc: int
        See :attr:`pc` attribute.
    enr: GOTermEnrichment object
        See :attr:`enr` attribute.

    Attributes
    ----------
    genes: list or tuple of str
        The list of genes in the signatures. The ordering of the genes must
        correspond to the ordering of the rows in ``E``.
    E: ndarray of floats
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
    enr: GOTermEnrichment object
        The result of the XL-mHG test that was conducted after ranking the
        genes based on their principal component loadings. See the
        `GOTermEnrichment` class for details.
    """

    _abbrev = [('positive ','pos. '),('negative ','neg. '),('interferon-','IFN-'),('proliferation','prolif.'),('signaling','signal.')]

    def __init__(self,genes,E,pc,enr,label=None):
        self.genes = tuple(genes) # genes in the signature (NOT equal to self.enr.genes, which contains the gene names corresponding to all the 1's)
        self.E = E # expression of the genes in the signture, with ordering matching that of self.genes
        self.pc = pc # principal component (sign indicates whether ordering was ascending or descending)
        self.enr = enr # GO enrichment this signature is based on

    def __repr__(self):
        return '<GOPCASignature: label="%s", pc=%d, es=%.1f; %s>' \
                %(self.label,self.pc,self.escore,repr(self.enr))

    def __str__(self):
        return '<GO-PCA Signature "%s" (PC %d / E-score %.1fx / %s)>' \
                %(self.label,self.pc,self.escore, str(self.enr))

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self,other):
        if type(self) is not type(other):
            return False
        elif repr(self) == repr(other):
            return True
        else:
            return False

    @property
    def term(self):
        """ The GO term that the signature is based on. """
        return self.enr.term

    @property
    def term_id(self):
        return self.enr.term[0]

    @property
    def term_name(self):
        return self.enr.term[3]

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
        return self.enr.ranks[self.k-1]

    @property
    def N(self):
        """ The total number of genes in the data. """
        return self.enr.N

    @property
    def label(self):
        return self.get_label(include_id=False)

    @property
    def median_correlation(self):
        C = np.corrcoef(self.E)
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
                ['term_id',['GO Term ID',r'%s']],
                ['k',['k',r'%d']],['K',['K',r'%d']],
                ['pval',['P-value',r'%.1e']],
                ['escore',['E-score (psi=%.1e)' %(self.escore_pval_thresh),r'%.1f']],
                ['median_correlation',['Median Correlation',r'%.2f']],
                ['gene_list',['Genes',r'%s']]
        ])
        od = OrderedDict([v[0],v[1] % (getattr(self,k))] for k,v in elements.iteritems())
        return od

    def get_label(self,max_name_length=0,include_stats=True,include_id=True,include_pval=False,include_collection=True):
        enr = self.enr

        term = enr.term
        term_name = term[3]
        for abb in self._abbrev:
            term_name = re.sub(abb[0],abb[1],term_name)
        if max_name_length > 0 and len(term_name) > max_name_length:
            term_name = term_name[:(max_name_length-3)] + '...'

        term_str = term_name
        if include_collection:
            term_str = '%s: %s' %(term[2],term_str)

        if include_id:
            term_str = term_str + ' (%s)' %(term[0])

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

        return term_str + stats_str


