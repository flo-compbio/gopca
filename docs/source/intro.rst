Introduction
============

This section provides a brief summary of *GO-PCA*, a list of key features of
the Python implementation, and a link to Demo notebooks.

What is *GO-PCA*?
-----------------

*GO-PCA* is an unsupervised method to explore gene expression data using prior
knowledge. Briefly, *GO-PCA* combines `principal component analysis (PCA)`__
with `nonparametric GO enrichment analysis`__ in order to define
**signatures**, i.e., small sets of genes that are both strongly correlated and
closely functionally related.

__ pca_
__ go_enrich_

The expression profiles of all signatures generated can be conveniently
visualized as a heat map. This visualization, referred to as the
**signature matrix**, aims to provide a systematic and easily interpretable
view of biologically relevant expression patterns in the data. Together with
other *GO-PCA* visualizations, it can serve as a powerful starting point for
exploratory data analysis and hypothesis generation. The method is described in
detail in an `open-access research article`__.

__ go_pca_paper_

.. _pca: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _go_enrich: https://dx.doi.org/10.1186/1471-2105-10-48
.. _go_pca_paper: https://dx.doi.org/10.1371/journal.pone.0143196


Key features
------------

GO-PCA is implemented in Python 2.7, a `high-level programming language`__ that
is widely used in both scientific and non-scientific settings. The key features
of GO-PCA are:

- Accessibility and transparency: GO-PCA is `free and open-source software`__.
- Cross-platform compatibility: GO-PCA can be easily
  `installed <install>` on Windows, OS X, and Linux, and runs under both
  Python 2.7.x and 3.5.x.
- Simple command-line interface: GO-PCA can be
  run directly from the command-line (`go-pca.py`), and command-line
  scripts can be used to generate output files containing the signatures
  created in tab-separated text (\*.tsv) or Excel spreadsheet (\*.xlsx) format.
- Powerful Python API (documentation forthcoming): The GO-PCA Python API
  can be used to create high-quality figures displaying the signature matrix
  or individual matrices in detail. This API in turn relies on the powerful
  and open-source `plotly`__ plotting library.
- Speed: GO-PCA takes about 60 seconds to run on the ``DMAP`` dataset,
  consisting  of ~8,000 genes and ~200 samples. The most computationally
  intensive part of GO-PCA (GO enrichment analysis using the XL-mHG test)
  is implemented in `Cython`__, a Python extension which produces efficient
  C code.
- Reproducibility: GO-PCA is a deterministic algorithm, and supports the
  calculation of `MD5 hash values`__ for all input, configuration, and output
  data. These values make it easy to establish e.g. whether two GO-PCA runs
  used identical parameter settings.
- Extensibility: GO-PCA's code is modular and well-documented, making it
  straightforward to implement modifications, new features and extensions.

__ python_
__ foss_
__ plotly_
__ cython_
__ md5_

.. _python: https://www.python.org/
.. _foss: https://en.wikipedia.org/wiki/Free_and_open-source_software
.. _plotly: https://plot.ly/
.. _cython: http://cython.org/A
.. _md5: https://en.wikipedia.org/wiki/MD5


Demos
-----

Demos of GO-PCA in action can be found in a `separate GitHub repository`__.
Note: These demos were created using an older version of this package
that relied on the `matplotlib`__ library for plotting.

__ demos_
__ matplotlib_

.. _demos: https://github.com/flo-compbio/gopca-demos

.. _matplotlib: http://matplotlib.org/
