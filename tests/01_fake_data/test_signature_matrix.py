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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# from builtins import *
# from builtins import open
from builtins import str as text

from copy import deepcopy

import pytest

from plotly import graph_objs as go

from genometools.expression.visualize import ExpHeatmap
from gopca import GOPCASignatureMatrix

def test_basic(my_sig_matrix):
    assert isinstance(my_sig_matrix, GOPCASignatureMatrix)
    assert isinstance(repr(my_sig_matrix), str)
    assert isinstance(str(my_sig_matrix), str)
    assert isinstance(text(my_sig_matrix), text)
    assert isinstance(my_sig_matrix.hash, text)

    other = deepcopy(my_sig_matrix)
    assert other is not my_sig_matrix
    assert other == my_sig_matrix
    other = GOPCASignatureMatrix.from_signatures(
        my_sig_matrix.signatures.tolist()[:-1],
        cluster_signatures=False)
    # signature matrix with only one signature
    assert other != my_sig_matrix


def test_heatmap(my_sig_matrix, my_signature):
    sig_matrix = my_sig_matrix
    sig1 = sig_matrix.get_signature(my_signature.gene_set.name)

    hl_col = 'blue'

    heatmap = sig_matrix.get_heatmap(
        highlight_sig={sig1: hl_col},
        colorbar_label=r'Median-centered expression (log<sub>2</sub>-RPKM)',
        matrix_kw={'cluster_samples': False},
    )
    assert isinstance(heatmap, ExpHeatmap)

    fig = heatmap.get_figure(
        width=1350, height=800,
        margin_left=350,
        margin_bottom=100,
        font_size=12,
        show_sample_labels=True,
    )
    assert isinstance(fig, go.graph_objs.Figure)
