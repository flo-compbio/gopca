# from __future__ import (absolute_import, division,
#                        print_function)
from builtins import str as text

import pkg_resources

from .params import GOPCAParams
from .config import GOPCAConfig
from .signature import GOPCASignature
from .signature_matrix import GOPCASignatureMatrix
from .run import GOPCARun
from .gopca import GOPCA

__version__ = text(pkg_resources.require('gopca')[0].version)

__all__ = ['GOPCAParams', 'GOPCAConfig', 'GOPCA', 'GOPCARun',
           'GOPCASignatureMatrix', 'GOPCASignature']
