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

Feedback
--------

Please feel free to email me with bug reports, questions, comments, and suggestions: florian (dot) wagner (at) duke (dot) edu.

