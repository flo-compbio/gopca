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

import pytest

from plotly import graph_objs as go

from genometools.expression.visualize import ExpHeatmap

def test_heatmap(my_sig_matrix):
    sig1 = my_sig_matrix.get_signature('regulation of os')
    sig2 = my_sig_matrix.get_signature('body morpho')
    sig3 = my_sig_matrix.get_signature('cell fate')
    sig4 = my_sig_matrix.get_signature('DNA amplification')
    sig5 = my_sig_matrix.get_signature('Notch')
    sig6 = my_sig_matrix.get_signature('P granule')
    sig7 = my_sig_matrix.get_signature('photoreceptor')

    hl_col = 'blue'

    heatmap = my_sig_matrix.get_heatmap(
        highlight_sig={sig1: hl_col, sig2: hl_col, sig3: hl_col, sig4: hl_col,
                       sig5: hl_col, sig6: hl_col, sig7: hl_col},
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
    #print(type(fig))
    #iplot(fig)
