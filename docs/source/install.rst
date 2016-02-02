Installation
============

GO-PCA is available in two versions: The latest official release (``latest``), and the current development version (``develop``). The official release is more stable, but the development version might contain some features that haven't made it into an official release yet.

Installing the latest GO-PCA release
------------------------------------

The latest GO-PCA release can be `found on PyPI`__, the Python package Index.

__ pypi_

The recommended installation procedure for GO-PCA differs slightly, depending
on whether you use a Windows, Linux, or OS X operating system:

- `Installation on Windows <install_windows>`
- `Installation on Ubuntu Linux <install_ubuntu>`
- ``Installation on OS X`` (to-do! Install Anaconda, as for Windows. => Which C compiler?)


.. _pypi: https://pypi.python.org/pypi/gopca

Installing the development version of GO-PCA
--------------------------------------------

The recommended installation procedures for the development version of GO-PCA
are identical to those for the latest release, except for the last step.
Instead of installing GO-PCA from PyPI, check out or download the `"develop"
branch`__ of the GO-PCA GitHub repository, and then run (on Ubuntu):

.. code-block:: bash

    $(env) cd gopca
    $(env) pip install -e .


__ develop_

.. _develop: https://github.com/flo-compbio/gopca/tree/develop
