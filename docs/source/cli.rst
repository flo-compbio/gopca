Command-line interface
======================

The GO-PCA command-line interface (CLI) consists of individual scripts that
can be used to process and visualize the results of a GO-PCA run.

.. contents:: Contents
    :depth: 2
    :local:
    :backlinks: none


.. _extract_go_gene_sets:

Generate custom GO-derived gene sets: ``gopca_extract_go_gene_sets.py``
-----------------------------------------------------------------------

Generating custom GO-derived gene sets for use with GO-PCA is a two-step
process: First, the script ``ensembl_extract_protein_coding_genes.py`` from
the `genometools` package has to be used to create a tab-delimited text file
with a list of protein-coding genes. The input for this script is an Ensembl
GTF file (see the "Gene sets" column on Ensembl's `FTP Download`__ page):

__ ensembl_download_


.. code-block:: bash
    
   ensembl_extract_protein_coding_genes.py -a [gtf_file] -o [output_file]


The output file can then be used as the "gene file" (``-g``) for the script
``gopca_extract_go_gene_sets.py``.

.. argparse::
   :ref: gopca.cli.extract_go_gene_sets.get_argument_parser
   :prog: gopca_extract_go_gene_sets.py


.. _ensembl_download: http://www.ensembl.org/info/data/ftp/index.html


.. _go_pca:

Running GO-PCA: ``go-pca.py``
-----------------------------

.. ".. code-block:: bash
    
    go-pca.py -g [gene_file] -a [annotation_file] -t [ontology_file] -e [expression_file] -o [output_file]

``go-pca.py`` is the command to run *GO-PCA*. All parameters can either be
spcefied directly on the command line, or in a separate configuration file,
using the ``-c`` option.

.. note::

  The configuration file is expected to follow the Windows "INI-style" format,
  with a single "[GO-PCA]" section, followed by "parameter=value" entries. 
  If a configuration file is given, and a parameter is set both in the
  configuration file and on the command line, the command line setting takes
  precedence.

The only required parameters are:::

 -e  (The expression file.)
 -s  (The gene set file.)
 -o  (The output file.)

However, if the expression matrix is not pre-filtered to only contain expressed
genes, it is also highly advisable to specify the ``-G`` option.

.. argparse::
   :ref: gopca.cli.go_pca.get_argument_parser
   :prog: go-pca.py


Inspecting the results: ``gopca_print_info.py``
-----------------------------------------------

In order to simply get a summary of the results contained in a particular
GO-PCA result file, the ``gopca_print_info.py`` command can be used. It prints
things like the number of principal components analyzed, the number of
signatures generated etc.

.. Hide
   
   It also outputs a list of all parameter settings
   used, as well as the names and MD5 hashsums of all input files.

.. argparse::
   :ref: gopca.cli.print_info.get_argument_parser
   :prog: gopca_print_info.py


.. _extract_signatures:


Extracting the signature matrix (as tab-delimited text file): ``gopca_extract_signature_matrix.py``
---------------------------------------------------------------------------------------------------

This command generates a tab-delimited text file which contains a matrix with
the signature expression values for each signature and each sample. (This is
the data visualized by the ``gopca_plot_signature_matrix.py`` command).

.. argparse::
   :ref: gopca.cli.extract_signature_matrix.get_argument_parser
   :prog: gopca_extract_signature_matrix.py


Plotting the signature matrix as a heatmap: ``gopca_plot_signature_matrix.py``
------------------------------------------------------------------------------

This command generates an interactive plot (embedded into an HTML file) of the
GO-PCA signature matrix, visualized as a heatmap.

The HTML file also allows exporting the figure to the PNG format.

.. argparse::
   :ref: gopca.cli.plot_signature_matrix.get_argument_parser
   :prog: gopca_plot_signature_matrix.py


Extracting the signatures (as tab-delimited text file): ``gopca_extract_signatures.py``
---------------------------------------------------------------------------------------

This command generates a tab-delimited text file in which each row corresponds
to a signature. The columns contain detailed information for each signature,
e.g., the gene set enrichment it was based on, and the list of genes contained in it.

.. argparse::
   :ref: gopca.cli.extract_signatures.get_argument_parser
   :prog: gopca_extract_signatures.py


Extracting the signatures (as Excel spreadsheet): ``gopca_extract_signatures_excel.py``
---------------------------------------------------------------------------------------

This command generates a file with the same information as
``gopca_extract_signatures.py``, but in the form of an Excel spreadsheet.

.. argparse::
   :ref: gopca.cli.extract_signatures_excel.get_argument_parser
   :prog: gopca_extract_signatures_excel.py

.. Hide
   
    Converting the results to MATLAB format: ``gopca_convert_to_matlab.py``
    -----------------------------------------------------------------------
   
    This command converts a GO-PCA result file to MATLAB format, using scipy's
    `io.savemat` function command from the `scipy` package.
   
    .. argparse::
    :ref: gopca.cli.convert_to_matlab.get_argument_parser
    :prog: gopca_convert_to_matlab.py

.. Hide
   
    Filtering the signatures: ``gopca_filter_signatures.py``
    --------------------------------------------------------
   
    GO-PCA tends to generate some highly correlated signatures that represent the
    same underlying signal. To some extent, this redundancy is intentional, as the
    different signature labels offer users alternative interpretations for the
    biological relevance of the underlying signal. However, sometimes these
    redundant signatures get in the way to result in an excessively long (tall)
    signature matrix that is difficult to read. In these cases, the
    ``gopca_filter_signature.py`` command can generate a reduced set of signatures
    so that their pair-wise correlation coefficients do not exceed a certain value.
    This can effectively remove highly correlated signatures.
   
    .. argparse::
    :ref: gopca.cli.filter_signatures.get_argument_parser
    :prog: gopca_filter_signatures.py

.. Hide
   
    Combining the signatures from two or more GO-PCA runs: ``gopca_combine_signatures.py``
    --------------------------------------------------------------------------------------

    This command does exactly what the name implies: It combines the signatures
    contained in two or more individual GO-PCA result files into a single, new
    result file.

    .. argparse::
    :ref: gopca.cli.combine_signatures.get_argument_parser
    :prog: gopca_combine_signatures.py
