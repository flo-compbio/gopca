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
	
		$ sudo apt-get install python-pip cython python-numpy python-matplotlib python-scikits-learn ipython ipython-notebook

3. Make sure the Ubuntu package python-sphinx is *not* installed:

	.. code-block:: bash
	
		$ sudo apt-get remove python-sphinx

	(The reason the package needs to be uninstalled is that it is an older version that conflicts with the version required by GO-PCA.)

4. Install the xlmhg python package using pip:

	.. code-block:: bash
	
		$ sudo pip install xlmhg

5. Install GO-PCA using pip (including all remaining dependencies):
	
	.. code-block:: bash
	
		$ sudo pip install gopca


Installation Instructions for Microsoft Windows
-----------------------------------------------

1. `Download <http://continuum.io/downloads>`_ and install the "Anaconda" Python distribution (if you haven't already). This distribution includes all third-party packages required by GO-PCA.

2. `Download <http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_ and install the Microsoft Visual C++ Compiler for Python 2.7. This compiler is required for installing the xlmhg python package (Step 3).

Then, on the command line (Start -> Run -> "cmd"):

3. Install the xlmhg python package using pip:
	
	.. code-block:: bash
	
		> pip install xlmhg

4. Install GO-PCA using pip (including all remaining dependencies):

	.. code-block:: bash
	
		> pip install gopca
