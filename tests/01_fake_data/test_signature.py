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
from builtins import str as text
# from builtins import int as newint

from copy import deepcopy

import pytest
import numpy as np
from plotly import graph_objs as go

from genometools.basic import GeneSet
from genometools.expression import ExpMatrix
from genometools.enrichment import RankBasedGSEResult
from genometools.expression.visualize import ExpHeatmap

from gopca import GOPCASignature


def test_basic(my_signature):
    assert isinstance(my_signature, GOPCASignature)
    assert isinstance(repr(my_signature), str)
    assert isinstance(str(my_signature), str)
    assert isinstance(text(my_signature), text)
    assert isinstance(my_signature.hash, text)

    other = deepcopy(my_signature)
    assert other is not my_signature
    assert other == my_signature
    other.pc = abs(my_signature.pc) + 1
    assert other != my_signature

def test_label(my_signature):
    assert isinstance(my_signature.label, text)
    label = my_signature.get_label(max_name_length=5, include_stats=True,
                                   include_id=True, include_pval=True,
                                   include_coll=True)
    assert isinstance(label, text)

def test_heatmap(my_signature):
    heatmap = my_signature.get_heatmap(
        cluster_samples=False,
        colorbar_label=r'Median-centered expression (log<sub>2</sub>-RPKM)',
    )
    assert isinstance(heatmap, ExpHeatmap)

    fig = heatmap.get_figure(
        emin=-3, emax=3,
        height=800, font_size=12, title_font_size=18, show_sample_labels=True,
        margin_bottom=100,
    )
    assert isinstance(fig, go.graph_objs.Figure)