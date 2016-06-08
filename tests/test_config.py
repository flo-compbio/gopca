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

from gopca import GOPCAConfig

import six
if six.PY2:
    import cPickle as pickle
else:
    import pickle


@pytest.fixture
def my_config():
    config = GOPCAConfig()
    return config


def test_basic(my_config):
    assert isinstance(my_config, GOPCAConfig)
    assert isinstance(repr(my_config), str)
    assert isinstance(str(my_config), str)
    assert isinstance(text(my_config), text)
    assert isinstance(my_config.hash, text)

    other = deepcopy(my_config)
    assert other is not my_config
    assert other == my_config
    other.sel_var_genes = 1001
    assert other != my_config

    params = my_config.params
    assert params is not my_config.params
    assert params == my_config.params
    assert params == my_config.param_defaults
    assert isinstance(my_config.param_names, list)
    for n in my_config.param_names:
        assert isinstance(n, text)
    assert isinstance(my_config.param_strings, list)
    for s in my_config.param_strings:
        assert isinstance(s, text)


def test_check(my_config):
    assert my_config.check_params()

def test_pickle(my_config, tmpdir):
    path = text(tmpdir.join('gopca_config.pickle'))
    with open(path, 'wb') as ofh:
        pickle.dump(my_config, ofh, pickle.HIGHEST_PROTOCOL)
    with open(path, 'rb') as fh:
        config = pickle.load(fh)
    assert config.hash == my_config.hash

def test_access(my_config):
    """Tests the other access methods, besides ``config[param]``."""
    other = deepcopy(my_config)
    for n in sorted(my_config.param_names):
        assert n in my_config
        ref = my_config[n]
        assert getattr(my_config, n) == ref
        assert my_config.get_default(n) is not None


def test_assign(my_config):
    other = deepcopy(my_config)
    for n in sorted(my_config.param_names):
        other.set_param(n, my_config[n])
    assert other == my_config
    other.set_params(my_config.params)
    assert other == my_config


def test_set_param(my_config):
    other = deepcopy(my_config)
    other.set_param('pval_thresh', 2e-6)
    assert other != my_config
    other.set_params(my_config.params)
    assert other == my_config
    other['pval_thresh'] = 2e-6
    assert other != my_config
    other.reset_params()
    assert other == my_config


@pytest.mark.xfail(
    reason= ('Fails in Python 2.7 due to '
             'https://github.com/PythonCharmers/python-future/issues/118')
    )
def test_read_write(tmpdir, my_config):
    other = deepcopy(my_config)
    other.set_param('pval_thresh', 2e-6)
    output_file = tmpdir.join('config.ini').strpath
    other.write_ini(output_file)
    other2 = GOPCAConfig.read_ini(output_file)
    assert other == other2


def test_error(my_config):
    with pytest.raises(AttributeError):
        test = my_config['hello world']
    with pytest.raises(AttributeError):
        test = my_config.get_param_default('hello world')