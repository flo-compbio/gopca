import pkg_resources

from gopca.input import GOPCAInput
from gopca.signature import GOPCASignature
from gopca.output import GOPCAOutput
from gopca.go_pca import GOPCA

__version__ = pkg_resources.require('gopca')[0].version

__all__ = ['GOPCAInput','GOPCA','GOPCAOutput']
