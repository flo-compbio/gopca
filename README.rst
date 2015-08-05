GO-PCA
======

GO-PCA is an unsupervised method for the exploratory analysis of heterogeneous transcriptomic data. It uses prior knowledge, in the form of `Gene Ontology (GO) <http://geneontology.org/>`_ anntoations, in combination with `principal component analysis (PCA) <https://en.wikipedia.org/wiki/Principal_component_analysis>`_, in order to generate *signatures*, i.e., small sets of genes that are both strongly correlated and closely functionally related. It then produces a *signature matrix*, which contains the expression profiles of all signatures across all samples, based on the untransformed data.

Installation
------------

GO-PCA is `available on PyPI <https://pypi.python.org/pypi/gopca>`_, the Python Package Index, and can be installed using ``pip``, the Python package manager:

.. code-block:: bash

	$ pip install gopca


Installation of Prerequisites
-----------------------------

GO-PCA depends on the following third-party Python packages:

- NumPy
- SciPy
- Scikit-Learn
- Cython

To-Do!

How to Cite GO-PCA
------------------

If you use GO-PCA in your research, please cite `Wagner (2015) <http://dx.doi.org/10.1101/018705>`_.

Documentation
-------------

In the works!

Copyright and License
---------------------

Copyright (c) 2015 Florian Wagner

::

  GO-PCA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License, Version 3,
  as published by the Free Software Foundation.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
