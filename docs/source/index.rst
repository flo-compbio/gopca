.. GO-PCA documentation master file, created by
   sphinx-quickstart on Fri Sep 11 11:35:24 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======
GO-PCA
======

.. "**Exploring Gene Expression Data Using Prior Knowledge**

.. image:: /_static/gopca.png
	:width: 500 px
	:align: center

`GO-PCA <https://github.com/flo-compbio/gopca>`_ is an unsupervised method to **explore gene expression data using prior knowledge**. It combines `principal component analysis (PCA) <https://en.wikipedia.org/wiki/Principal_component_analysis>`_ with `nonparametric GO enrichment analysis <http://dx.doi.org/10.1186/1471-2105-10-48>`_, in order to generate *signatures*, i.e., small sets of genes that are both strongly correlated and closely functionally related. It then produces a *signature matrix*, which contains the expression profiles of all signatures across all samples, based on the untransformed data.

Demo
----

Check out a `demonstration of GO-PCA <http://nbviewer.ipython.org/github/flo-compbio/gopca/blob/master/notebooks/GO-PCA_Demo.ipynb>`_ based on the `DMAP dataset <http://dx.doi.org/10.1016/j.cell.2011.01.004>`_.

Getting Started
---------------

GO-PCA is available as Python package from PyPI and has been tested under both Ubuntu Linux and Windows. Check out the :doc:`Installation section</install>`.


.. "Table of Contents
	-----------------

.. toctree::
	:maxdepth: 2
	:hidden:

	install
	running
	visual_scripts
	process_scripts
	modules


..	intro
	tutorial
	autodoc


.. "Indices and tables
	==================
	
	* :ref:`genindex`
	* :ref:`search`
	
	.. "* :ref:`modindex`
