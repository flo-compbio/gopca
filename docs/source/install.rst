Installation
============

List of Dependencies
--------------------

GO-PCA is implemented in Python 2.7, and directly depends on the following third-party Python packages:

- Cython
- NumPy
- scikit-learn (which in turn depends on SciPy)
- matplotlib
- sphinx (for documentation)

Both Ubuntu Linux and Windows offer convenient options for automatically satisfying these dependencies. GO-PCA further depends on three Python packages developed by me:

- genometools
- goparser
- xlmhg

These packages are available from the `Python Package Index <https://pypi.python.org/pypi>`_, and can be installed using pip, the Python package manager.

The following sections provide detailed instructions for installing GO-PCA and its dependencies on both Ubuntu Linux and Microsoft Windows.

Installation Instructions for Ubuntu Linux
-------------------------------------------

In a terminal window:

1. Make sure you have Python 2.7 installed:

	.. code-block:: bash

		$ python -V
		Python 2.7.6

	Ubuntu 14.04 (trusty) currently has Python version 2.7.6, but any Python 2.7 release should work.

2. Install GO-PCA dependencies: (pip, Cython, NumPy, scikit-learn, Matplotlib, and sphinx for building this documentation locally):
	
	.. code-block:: bash
	
		$ sudo apt-get install python-pip cython python-numpy python-matplotlib python-scikits-learn
	
2. Make sure you have pip installed:
	
	.. code-block:: bash
	
		$ pip -V
		pip 7.1.0 from [path] (python 2.7)
	
3. Install GO-PCA:
	
	.. code-block:: bash
	
		$ pip install gopca

Note: GO-PCA directly depends on the following third-party Python packages:

- Cython
- NumPy
- scikit-learn (which in turn depends on SciPy)
- matplotlib
- sphinx (for documentation)

.. pip will attempt to download and install the latest versions of these packages automatically from the `PyPI, Python Package Index <https://pypi.python.org>`_, but SciPy in particular has additional dependencies that 


Installation Instructions for Microsoft Windows
-----------------------------------------------

1. `Download <http://continuum.io/downloads>`_ and install the "Anaconda" Python distribution (if you haven't already).

2. This 

On the command line (Run-> "cmd"):
1. Download Microsoft Visual C++ Compiler for Python 2.7
