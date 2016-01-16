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

"""Module containing the `GOPCARun` class.
"""

import logging
from copy import deepcopy
import cPickle as pickle

from gopca import GOPCAConfig, GOPCAResult

logger = logging.getLogger(__name__)

class GOPCARun(object):
    """A GO-PCA run, consisting of configuration, run, and result data.

    Parameters
    ----------
    version: str
        The GO-PCA version.
    user_config: `gopca.GOPCAConfig`
        The configuration data provided by the user.
    timestamp: str
        The timestamp.
    result: `GOPCAResult`
        The result (signatures).
    exec_time: float
        The execution time (in seconds).
    """

    def __init__(self, version, user_config, timestamp, result, exec_time):
        # checks
        assert isinstance(version, (str, unicode))
        assert isinstance(user_config, GOPCAConfig)
        assert isinstance(timestamp, (str, unicode))
        assert isinstance(result, GOPCAResult)
        assert isinstance(exec_time, float)
        # initialization
        self.version = version
        self.user_config = deepcopy(user_config)
        self.timestamp = timestamp
        self.result = result
        self.exec_time = exec_time

    ### magic functions
    def __repr__(self):
        return '<GOPCARun (version=%s; timestamp=%s; exec_time=%.1f; hash=%d)>' \
                %(self.version, self.timestamp, hash(self))

    def __str__(self):
        return '<GOPCARun (version %s; %d signatures; %s)>' \
                %(self.result.q, self.timestamp)

    def __eq__(self,other):
        if type(self) is not type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __hash__(self):
        data = []
        data.append(self.user_config)
        data.append(self.version)
        data.append(self.timestamp)
        data.append(self.result)
        data.append(self.exec_time)
        return hash(tuple(data))
    ### end magic functions

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
        with open(path, 'wb') as ofh:
            pickle.dump(self, ofh, pickle.HIGHEST_PROTOCOL)
