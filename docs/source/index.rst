GO-PCA |version|
================

.. "**Exploring Gene Expression Data Using Prior Knowledge**

.. "image:: /_static/gopca.png
    :width: 500 px
    :align: right

.. figure:: /_static/gopca.png
    :width: 400 px
    :align: right
    
    Detail from a GO-PCA *signature plot*.

GO-PCA (`Wagner, 2015`__) is an unsupervised method to **explore gene
expression data using prior knowledge**. This is the documentation for
the `Python implementation of GO-PCA`__, which is **free and open-source
software** (see `License <license>`).

__ go_pca_paper_
__ go_pca_

Briefly, GO-PCA combines `principal component analysis (PCA)`__  with
`nonparametric GO enrichment analysis`__ in order to defined **signatures**,
i.e., small sets of genes that are both strongly correlated and closely
functionally related. It then plots the expression profiles of all
signatures in a heat map. This visualization, referred to as the **signature
matrix**, aims to provide a systematic and easily interpretable view
of biologically relevant expression patterns in the data. These and other
GO-PCA outputs can serve as a powerful starting point for data analysis and
hypothesis generation. The `GO-PCA paper`__ contains a detailed description of
the method.

__ pca_
__ go_enrich_
__ go_pca_paper_

.. _go_pca: https://github.com/flo-compbio/gopca
.. _go_pca_paper: https://dx.doi.org/10.1371/journal.pone.0143196
.. _pca: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _go_enrich: https://dx.doi.org/10.1186/1471-2105-10-48

Table of contents
-----------------

.. toctree::
    :maxdepth: 2

    Home <self>
    GO-PCA Demos <demos>
    Getting Started <getting_started>
    License <license>

..  intro
    tutorial
    Modules <modules>
    autodoc


.. "Indices and tables
    ==================
    
    * :ref:`genindex`
    * :ref:`search`
    
    .. "* :ref:`modindex`
