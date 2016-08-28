Installing GO-PCA
=================

GO-PCA is available in two versions: The latest official release (``latest``),
and the current development version (``develop``). The official release is more
stable, but the development version might contain some features that haven't
made it into an official release yet.

Installing the latest release
-----------------------------

The latest GO-PCA release is `hosted on PyPI`__, the Python package Index.
The recommended procedure for installing GO-PCA differs slightly, depending on
whether you use a Windows, Linux, or Mac OS X operating system:

__ pypi_

- `Installation on Ubuntu Linux <install_ubuntu>`
- `Installation on Windows <install_windows>`
- `Installation on Mac OS X <install_mac>`


.. _pypi: https://pypi.python.org/pypi/gopca

Installing the development version
----------------------------------

The recommended installation procedures for the development version of GO-PCA
are identical to those for the latest release, except for the last step.
Instead of running ``pip install gopca``, run the following::

    pip install git+git://github.com/flo-compbio/gopca.git@develop


.. toctree::
    :hidden:
    :maxdepth: 2

    install_windows
    install_ubuntu