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
from builtins import str as text

from copy import deepcopy

import pytest

from gopca import GOPCAParams

import six
if six.PY2:
    import cPickle as pickle
else:
    import pickle


def test_basic(my_params):
    assert isinstance(my_params, GOPCAParams)
    assert isinstance(repr(my_params), str)
    assert isinstance(str(my_params), str)
    assert isinstance(text(my_params), text)
    assert isinstance(my_params.hash, text)

    other = deepcopy(my_params)
    assert other is not my_params
    assert other == my_params
    other.mHG_min = 7
    assert other != my_params

    params = my_params.params
    assert params is not my_params.params
    assert params == my_params.params
    assert params == my_params.param_defaults
    assert isinstance(my_params.param_names, list)
    for n in my_params.param_names:
        assert isinstance(n, text)
    assert isinstance(my_params.param_strings, list)
    for s in my_params.param_strings:
        assert isinstance(s, text)


def test_check(my_params):
    assert my_params.check_params()


def test_pickle(my_params, tmpdir):
    path = text(tmpdir.join('gopca_config.pickle'))
    with open(path, 'wb') as ofh:
        pickle.dump(my_params, ofh, pickle.HIGHEST_PROTOCOL)
    with open(path, 'rb') as fh:
        config = pickle.load(fh)
    assert config.hash == my_params.hash


def test_access(my_params):
    """Tests the other access methods, besides ``config[param]``."""
    other = deepcopy(my_params)
    for n in sorted(my_params.param_names):
        assert n in my_params
        ref = my_params[n]
        assert getattr(my_params, n) == ref
        assert my_params.get_default(n) is not None


def test_assign(my_params):
    other = deepcopy(my_params)
    for n in sorted(my_params.param_names):
        other.set_param(n, my_params[n])
    assert other == my_params
    other.set_params(my_params.params)
    assert other == my_params


def test_set_param(my_params):
    other = deepcopy(my_params)
    other.set_param('pval_thresh', 2e-6)
    assert other != my_params
    other.set_params(my_params.params)
    assert other == my_params
    other['pval_thresh'] = 2e-6
    assert other != my_params
    other.reset_params()
    assert other == my_params


@pytest.mark.xfail(reason=(
    'Fails in Python 2.7 due to '
    'https://github.com/PythonCharmers/python-future/issues/118'))
def test_read_write(tmpdir, my_params):
    other = deepcopy(my_params)
    other.set_param('pval_thresh', 2e-6)
    output_file = tmpdir.join('config.ini').strpath
    other.write_ini(output_file)
    other2 = GOPCAParams.read_ini(output_file)
    assert other == other2


def test_error(my_params):
    with pytest.raises(AttributeError):
        test = my_params['hello world']
    with pytest.raises(AttributeError):
        test = my_params.get_param_default('hello world')