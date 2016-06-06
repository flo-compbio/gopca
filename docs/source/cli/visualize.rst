Visualizing GO-PCA results
==========================

This section documents all commands for visualizing GO-PCA results, most
notably the *signature matrix* as well as individual signatures. All commands
not related to visualizations are `documented in the previous section
<process>`.

.. note::

  All visualization commands produce PNG images as output files, except for
  ``plot_all_signatures.py``,  which produces a PDF file with multiple pages
  (one for each signature).


.. contents:: Contents
    :depth: 2
    :local:
    :backlinks: none


.. _plot_signature_matrix:

Plotting the signature matrix: ``gopca_plot_signature_matrix.py``
-----------------------------------------------------------------

This command plots the GO-PCA signature matrix in PNG format. This command is
typically the first command to run after the GO-PCA run has finished.

.. argparse::
   :ref: gopca.cli.plot_signature_matrix.get_argument_parser
   :prog: gopca_plot_signature_matrix.py


Plotting an individual signature in detail: ``gopca_plot_signature.py``
-----------------------------------------------------------------------

This command generates a plot that provides detailed view of the expression
pattern of each gene in a specific signature.

.. argparse::
   :ref: gopca.cli.plot_signature.get_argument_parser
   :prog: gopca_plot_signature.py

Plotting all signatures in detail: ``gopca_plot_all_signature.py``
------------------------------------------------------------------

This command works like ``gopca_plot_signature.py``, applied to every signature
in the GO-PCA result. The output is a PDF that contains one signature plot per
page (the number of pages equals the number of signatures).

.. argparse::
   :ref: gopca.cli.plot_all_signatures.get_argument_parser
   :prog: gopca_plot_all_signatures.py

Ploting the term-by-PC matrix: ``gopca_plot_term_by_pc_matrix.py``
------------------------------------------------------------------

This command plots a heatmap of a matrix showing the GO enrichment p-values,
for each GO term that used to generate a signature, for each principal
component (PC) tested. This visualization provides additional transparency in
terms of what PCs appear associated with which GO terms.

This visualization was originally proposed by Dr. Meromit Singer from Dr. Aviv
Regev's group at the Broad Institute.

.. argparse::
   :ref: gopca.cli.plot_term_by_pc_matrix.get_argument_parser
   :prog: gopca_plot_term_by_pc_matrix.py


