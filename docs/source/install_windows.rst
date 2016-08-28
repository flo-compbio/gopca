Installation Instructions for Windows
=====================================

The easiest way to install GO-PCA on Windows consists of a two-step procedure:
First, install `Anaconda`__, a Python distribution provided by
Continuum Analytics. Second, install GO-PCA with `pip`__.

Here is a step-by-step guide:

1. `Download`__ and install Anaconda (~350 MB download). The Python 3.5 version
   of Anaconda is recommended, but 2.7 works, too. During the installation,
   make sure you leave the options "Add Anaconda to my PATH environment
   variable" and "Register Anaconda as my default Python 3.5" (or Python 2.7)
   enabled.

2. Open a command prompt (Start->Run->"Cmd").

3. Type ``pip install gopca`` and hit Enter.

__ anaconda_
__ pip_
__ download_anaconda_


Installing a specific version of GO-PCA
---------------------------------------

To install a specific version of GO-PCA, e.g., "2.0.0", run::

    pip install gopca==2.0.0


Why is Anaconda required?
-------------------------

Many third-party Python packages do not require a large, third-party Python
distribution like Anaconda for their installation on Windows (i.e., they can
be installed with the `official Python distribution`__). However, GO-PCA
depends on the *SciPy* library, which is difficult to install with ``pip`` on
Windows (the reasons for this are beyond the scope of this manual). Anaconda
comes with SciPy pre-installed, which means that pip does not need to
install the package when installing GO-PCA. This makes the installation of
GO-PCA much more straightforward.

__ python_download_


Installation instructions for advanced Python users
---------------------------------------------------

Advanced Python users are likely familiar with the concept of Python
"environments", and know how to use `conda`__ for creating and managing
environments. When installing GO-PCA into a separate conda environment on
Windows, this environment should have SciPy installed in it (e.g., use
something like ``conda create -n gopca_env python=3.5 scipy`` to create the
environment). Then, inside the environment, install GO-PCA using ``pip``, as
described above.

__ conda_

.. _anaconda: https://www.continuum.io/

.. _download_anaconda: https://www.continuum.io/downloads#windows

.. _python_download: https://www.python.org/downloads/

.. _pip: https://pip.pypa.io/en/stable/

.. _conda: http://conda.pydata.org/docs/