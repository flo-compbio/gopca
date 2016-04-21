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

"""Module containing the GOPCAPlotter class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging
from collections import OrderedDict

import pandas as pd
import numpy as np

from genometools.expression import ExpMatrix
from genometools.expression import plot as eplt

from . import GOPCAResult, GOPCASignature

logger = logging.getLogger(__name__)


class GOPCAPlotter(object):
    """A plotter for GO-PCA results."""

    def __init__(self, result):
        assert isinstance(result, GOPCAResult)

        self.result = result

    def plot_signature_matrix(
            self, sig_max_name_length=50, sig_include_id=False,
            emin=-3.0, emax=3.0,
            margin_left=150, margin_bottom=50,
            show_sample_labels=False, matrix_kw=None, **kwargs):
        """Plot a GO-PCA signature matrix as a heat map."""

        if matrix_kw is not None:
            assert isinstance(matrix_kw, (dict, OrderedDict))

        if matrix_kw is None:
            matrix_kw = {}

        matrix_kw['max_name_length'] = sig_max_name_length
        matrix_kw['include_id'] = sig_include_id

        result = self.result

        # generate signature matrix
        S, a_signatures, a_samples = result.get_signature_matrix(**matrix_kw)

        # plot heat map
        fig = eplt.get_heatmap(
            S, yaxis_label='Signatures',
            colorbar_label='Standardized Expression',
            emin=emin, emax=emax,
            show_sample_labels=show_sample_labels,
            margin_left=margin_left, margin_bottom=margin_bottom,
            **kwargs
        )

        return fig, S, a_signatures, a_samples

    def plot_signature(
            self, name=None, pc=None, i=None, sig=None,
            standardize=False, center=True, cluster_samples=True,
            include_id=False,
            emin=None, emax=None,
            margin_left=70, margin_bottom=50, margin_top=50,
            show_sample_labels=False, matrix_kw=None, **kwargs):

        if sig is not None:
            assert isinstance(sig, GOPCASignature)
        else:
            assert name is not None, 'No signature and no signature name given!'
            assert isinstance(name, str)

        if sig is None:
            # look up signature
            sig = self.result.get_signature(name, pc, i)

        if matrix_kw is None:
            matrix_kw = {}

        # get signature expression matrix
        E, a_genes, _ = sig.get_expression_matrix(
            standardize=standardize,
            center=center,
            **matrix_kw)

        cb_label = 'Centered Expression'
        if standardize:
            cb_label = 'Standardized Expression'
    
        if cluster_samples:
            # use entire signature matrix as the basis for clustering
            S, _, a_samples = self.result.get_signature_matrix()
            E = E.iloc[:, a_samples]

        # add a "Signature" row to the top
        # (representing the signature expression)
        title = sig.get_label(include_id=include_id)
        mean = np.mean(E.X, axis=0)
        header_row = ExpMatrix(genes=['Signature'], samples=E.samples,
                               X=np.atleast_2d(mean))
        E = pd.concat([header_row, E], axis=0)

        # plot heat map
        fig = eplt.get_heatmap(
            E, title=title, yaxis_label='Genes',
            colorbar_label=cb_label,
            emin=emin, emax=emax,
            show_sample_labels=show_sample_labels,
            margin_left=margin_left,
            margin_bottom=margin_bottom, margin_top=margin_top,
            **kwargs
        )

        return fig
