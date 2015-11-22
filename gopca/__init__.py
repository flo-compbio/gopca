import pkg_resources

from go_pca_objects import GOPCAConfig, GOPCA

__version__ = pkg_resources.require('gopca')[0].version

__all__ = ['GOPCAConfig','GOPCA']
