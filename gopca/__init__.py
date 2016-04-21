from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import pkg_resources

from .config import GOPCAConfig
from .go_pca import GOPCA
from .run import GOPCARun
from .result import GOPCAResult
from .signature import GOPCASignature
from .plotter import GOPCAPlotter

__version__ = str(pkg_resources.require('gopca')[0].version)

__all__ = ['GOPCAConfig', 'GOPCA', 'GOPCARun', 'GOPCAPlotter',
           'GOPCAResult', 'GOPCASignature']
