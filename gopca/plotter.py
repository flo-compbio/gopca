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

import logging
from collections import OrderedDict

import pandas as pd
import numpy as np

from bokeh.models import HoverTool

from genometools.expression import ExpMatrix
from genometools.expression import plot as eplt

from . import GOPCAResult, GOPCASignature

logger = logging.getLogger(__name__)

class GOPCAPlotter(object):
    """A plotter for GO-PCA results."""

    def __init__(self):
        pass

    def plot_signature_matrix(
            self, result, sig_max_name_length = 50, sig_include_id = False,
            width = 800, height = 400,
            cmap = None, vmin = None, vmax = None,
            font = None, font_size = None,
            show_sample_labels = False, matrix_kw = None):
        """Plot a GO-PCA signature matrix as a heat map."""

        assert isinstance(result, GOPCAResult)

        if matrix_kw is None:
            matrix_kw = {}

        assert isinstance(matrix_kw, (dict, OrderedDict))

        matrix_kw['max_name_length'] = sig_max_name_length
        matrix_kw['include_id'] = sig_include_id

        # generate signature matrix
        S, _, _ = result.get_signature_matrix(**matrix_kw)

        # plot heat map
        hm = eplt.plot_heatmap(
                S, cmap = cmap, yaxis_label = 'Signatures', 
                width = width, height = height,
                vmin = -3.0, vmax = 3.0, font_size = '8pt'
        )

        return hm

    def plot_signature(
            self, sig, include_id = False,
            cmap = None, vmin = None, vmax = None,
            width = 800, height = 400,
            font = None, font_size = None,
            show_sample_labels = False, sig_kw = None):
        """Plot the gene expression matrix of a GO-PCA signature as a heat map."""

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

