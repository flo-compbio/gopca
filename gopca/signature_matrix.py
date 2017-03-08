#  (c) 2015, 2016 Florian Wagner
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

"""Module containing the `GOPCASignatureMatrix` class.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import logging
import hashlib
import copy
from collections import Iterable

import six
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

from genometools.expression import ExpProfile, ExpMatrix
from genometools.expression import cluster
from genometools.expression.visualize import ExpHeatmap, HeatmapGeneAnnotation

# from .config import GOPCAParams
from . import GOPCASignature
# from gopca import util

if six.PY2:
    import cPickle as pickle
else:
    import pickle

logger = logging.getLogger(__name__)


class GOPCASignatureMatrix(ExpMatrix):
    """A GO-PCA signature matrix (the result of a GO-PCA run).

    """
    def __init__(self, *args, **kwargs):
        return ExpMatrix.__init__(self, *args, **kwargs)


    @classmethod
    def from_signatures(
            cls,
            signatures,
            standardize=False,
            center=True,
            use_median=True,
            cluster_signatures=True,
            signature_cluster_metric='correlation',
            cluster_samples=True,
            sample_cluster_metric='euclidean',
            cluster_method='average'
        ):
        """Generate a GO-PCA signature matrix from individual signatures.

        The GO-PCA signature matrix contains the expression levels of all
        signatures (rows) generated, across all samples (columns) in the
        analysis. See the documentation of the `GOPCASignature` class for
        details on how signature expression levels are calculated.

        Parameters
        ----------
        signatures: Iterable of `GOPCASignature`
        The signatures generated.
        """
        # TODO: finish docstring
        assert isinstance(signatures, Iterable)
        assert isinstance(standardize, bool)
        assert isinstance(center, bool)
        assert isinstance(use_median, bool)
        assert isinstance(cluster_signatures, bool)
        assert isinstance(cluster_samples, bool)

        ### generate the expression matrix
        matrix = ExpMatrix(pd.concat(
            [sig.get_expression(standardize=standardize, center=center,
                                use_median=use_median)
             for sig in signatures],
            axis=1
        ).T)
        matrix.genes.name = 'Signatures'
        matrix.samples.name = 'Samples'

        if matrix.p == 1:
            cluster_signatures = False
            cluster_samples = False

        ### clustering
        if cluster_signatures:
            # cluster signatures
            matrix = cluster.cluster_genes(
                matrix, metric=signature_cluster_metric, method=cluster_method
            )

        order_samples = None
        if cluster_samples:
            # cluster samples
            matrix = cluster.cluster_samples(
                matrix, metric=sample_cluster_metric, method=cluster_method
            )

        return cls(matrix)

    # magic functions
    def __repr__(self):
        return ('<GOPCASignatureMatrix instance (q=%d, hash="%s")>'
                % (self.q, self.hash))

    def __str__(self):
        return ('<GOPCASignatureMatrix instance with %d signatures (n=%d)>'
                % (self.q, self.n))

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
    def _constructor(self):
        return GOPCASignatureMatrix

    @property
    def _constructor_sliced(self):
        return ExpProfile

    @property
    def hash(self):
        """An MD5 hash string for the signature."""
        data_str = ';'.join([
            str(repr(var)) for var in
            [self.samples, self.signatures]
        ])
        data = data_str.encode('UTF-8')
        return str(hashlib.md5(data).hexdigest())

    @property
    def q(self):
        """The number of signatures in the matrix."""
        return self.p

    @property
    def n(self):
        """The number of samples in the matrix."""
        return len(self.samples)

    @property
    def signatures(self):
        """This returns the list of signatures (as a `pandas.Index` object)."""
        return self.genes

    def get_signature(self, name, pc=None, i=None):
        """Look up a signature by name, PC, and index.

        """
        # TODO: Finish docstring.
        # find all signatures that start with name and (optionally)
        # were generated by a certain PC

        # type checks
        assert isinstance(name, (str, _oldstr))
        if pc is not None:
            assert isinstance(pc, int) and pc >= 1
        if i is not None:
            assert isinstance(i, int)

        if i is None:
            i = 0

        found = []
        for sig in self.signatures:
            if sig.gene_set.name.startswith(name):
                if pc is None or abs(sig.pc) == pc:
                    found.append(sig)

        if not found:
            raise ValueError('Signature not found!')

        if i >= len(found):
            raise ValueError('Index "%d" out of bounds.' % i)

        return found[i]

    def filter_collection_signatures(self, corr_thresh=0.9, source=None):
        """Filter signatures by collection."""

        assert isinstance(corr_thresh, float)
        if source is not None:
            assert isinstance(source, (str, _oldstr))

        sig_matrix = copy.deepcopy(self)
        signatures = sig_matrix.signatures

        # sort signatures first by PC, then by E-score
        sig_abs_pcs = np.absolute(np.int64([sig.pc for sig in signatures]))
        sig_escore = np.float64([sig.escore for sig in signatures])
        a = np.lexsort([-sig_escore, sig_abs_pcs])

        # filtering
        q = len(signatures)
        sel = np.ones(q, dtype=np.bool_)
        for pos, i in enumerate(a):
            sig = signatures[i]

            if not sel[i]:
                # signature is already excluded => skip
                continue

            if source is not None and sig.gene_set.source != source:
                # signature doesn't have the selected source => ignore
                continue

            src = sig.gene_set.source
            coll = sig.gene_set.collection

            for i2 in a[(pos+1):]:
                other = signatures[i2]

                if other.gene_set.source != src:
                    # other signature does not have the same source => ignore
                    continue
                elif other.gene_set.collection != coll:
                    # other signature does not have the same collection
                    # => ignore
                    continue

                # otherwise, exclude the other signature if it's too highly
                # correlated
                r = pearsonr(sig.expression, other.expression)[0]
                if r >= corr_thresh:
                    # signature is too highly correlated => exclude
                    logger.debug('Excluding signature "%s" due to correlation'
                                 'with "%s".', other.label, sig.label)
                    sel[i2] = False

        sel = np.nonzero(sel)[0]
        sig_matrix = GOPCASignatureMatrix.from_signatures(
            sig_matrix.signatures[i] for i in sel
        )
        return sig_matrix

    @property
    def signature_labels(self):
        """The list of signature labels."""
        sig_labels = self.get_signature_labels()
        return sig_labels


    def get_signature_labels(self, **kwargs):
        """Generate a list of GO-PCA signature labels (convenience function).

        Returns
        -------
        list of str
            List of signature labels.
        """
        sig_labels = [sig.get_label(**kwargs) for sig in self.signatures]
        return sig_labels

    def get_heatmap(self, max_name_length=40, include_id=False,
                    highlight_sig=None, highlight_source=None,
                    matrix_kw=None,
                    colorbar_label=('Signature expression<br>'
                                    '(log<sub>2</sub>-scale)'),
                    annotation_transparency=0.8):
        """Generate an `ExpHeatMap` instance."""

        if matrix_kw is None:
            matrix_kw = {}

        if highlight_sig is None:
            highlight_sig = {}

        if highlight_source is None:
            highlight_source = {}

        assert isinstance(highlight_sig, dict)
        assert isinstance(highlight_source, dict)
        assert isinstance(matrix_kw, dict)
        if colorbar_label is not None:
            assert isinstance(colorbar_label, (str, _oldstr))

        # generate expresssion matrix
        matrix = self.copy()

        # generate signature labels
        #signatures = matrix.index.values
        sig_labels = self.get_signature_labels(
            max_name_length=max_name_length,
            include_id=include_id)

        # add signature annotations
        sig_annotations = []

        # annotations for individual signatures
        for sig, color in highlight_sig.items():
            try:
                # pd.Index.get_loc() does not work with objects (bug?), so
                # we have to do it the slow way using np.nonzero
                i = np.nonzero(self.signatures == sig)[0][0]
                sig_annotations.append(
                    HeatmapGeneAnnotation(sig_labels[i], color,
                                          label=sig_labels[i])
                )
            except KeyError as e:
                raise e
                # raise ValueError('%s not found in signature matrix.'
                #                  % str(sig))

        # annotations for all signatures from certain gene set sources
        for src, color in highlight_source.items():
            for sig, label in zip(matrix.index, sig_labels):
                if sig.gene_set.source == src:
                    sig_annotations.append(
                        HeatmapGeneAnnotation(
                            label, color, label=label,
                            transparency=annotation_transparency)
                    )

        # replace signatures with labels in expression matrix
        matrix.index = sig_labels
        matrix.index.name = 'Signatures'

        # generate ExpHeatMap
        heatmap = ExpHeatmap(matrix, gene_annotations=sig_annotations,
                             colorbar_label=colorbar_label)

        return heatmap


    def get_figure(self, heatmap_kw=None, **kwargs):
        """Generate a plotly figure showing the signature matrix as a heatmap.

        This is a shortcut for
        ``SignatureMatrix.get_heatmap(...).get_figure(...)``.

        See :func:`ExpHeatmap.get_figure` for keyword arguments.

        Parameters
        ----------
        heatmap_kw: dict or None
            If not None, dictionary containing keyword arguments to be passed
            to the `ExpHeatmap` constructor.

        Returns
        -------
        `plotly.graph_objs.Figure`
            The plotly figure.
        """
        if heatmap_kw is not None:
            assert isinstance(heatmap_kw, dict)

        if heatmap_kw is None:
            heatmap_kw = {}

        width = kwargs.pop('width', 1350)
        height = kwargs.pop('height', 800)

        margin_left = kwargs.pop('margin_left', 350)
        margin_bottom = kwargs.pop('margin_bottom', 100)

        emin = kwargs.pop('emin', -3.0)
        emax = kwargs.pop('emax', 3.0)

        font_size = kwargs.pop('font_size', 10)
        title_font_size = kwargs.pop('title_font_size', font_size)

        return self.\
            get_heatmap(**heatmap_kw).\
            get_figure(width=width, height=height,
                       margin_left=margin_left, margin_bottom=margin_bottom,
                       emin=emin, emax=emax,
                       font_size=font_size, title_font_size=title_font_size,
                       **kwargs)


    def filter_signatures(self, corr_thresh, inplace=False):
        """Remove "redundant" signatures."""

        # checks
        assert isinstance(corr_thresh, (float, int))
        assert 0 < corr_thresh <= 1.0

        matrix = self
        if not inplace:
            matrix = matrix.copy()

        if corr_thresh == 1.0:
            # no filtering
            return matrix

        signatures = matrix.signatures
        matrix = matrix.expression_matrix

        # sort signatures first by PC, then by E-score
        sig_abs_pcs = np.absolute(np.int64([sig.pc for sig in signatures]))
        sig_escore = np.float64([sig.escore for sig in signatures])
        a = np.lexsort([-sig_escore, sig_abs_pcs])

        # filtering
        S = matrix.values
        q, n = matrix.shape
        sel = np.ones(q, dtype=np.bool_)
        for i in a:

            if not sel[i]:
                # already excluded
                continue

            for i2, sig in enumerate(signatures):
                if i == i2 or not sel[i2]:
                    continue
                # assert np.corrcoef(np.vstack([S[i,:],S[i2,:]])).shape == (2,2)
                if np.corrcoef(np.vstack([S[i, :], S[i2, :]]))[0, 1] >= \
                        corr_thresh:
                    logger.info(
                        'Excluding signature "%s" due to correlation with '
                        '"%s".',
                        signatures[i2].label, signatures[i].label)
                    sel[i2] = False

        sel = np.nonzero(sel)[0]
        matrix.signatures = [matrix.signatures[i] for i in sel]
        return matrix

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
