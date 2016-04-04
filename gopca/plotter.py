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
            self, sig_max_name_length = 50, sig_include_id = False,
            emin = -3.0, emax = 3.0, margin_left = 150, margin_bottom = 50,
            show_sample_labels = False, matrix_kw = None, **kwargs):
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

    """
    def plot_signature(
            self, sig, include_id = False,
            cmap = None, vmin = None, vmax = None,
            width = 800, height = 400,
            font = None, font_size = None,
            show_sample_labels = False, sig_kw = None):
        #Plot the gene expression matrix of a GO-PCA signature as a heat map.

        assert isinstance(sig, GOPCASignature)

        if sig_kw is None:
            sig_kw = {}

        E, _, _ = sig.get_expression_matrix(**sig_kw) # clustering etc.

        # add a "Signature" row to the top
        # (representing the signature expression)
        title = sig.get_label(include_id = include_id)
        mean = np.mean(E.X, axis = 0)
        header_row = ExpMatrix(genes = ['Signature'], samples = E.samples, X = np.atleast_2d(mean))
        E = pd.concat([header_row, E], axis = 0)

        hm = eplt.plot_heatmap(
                E, cmap = cmap, title = title, width = width, height = height
        )
        # extract expression matrix for signature genes
        
        return hm
    """
