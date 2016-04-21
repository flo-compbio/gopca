Processing GO-PCA results
=========================

This section documents all commands for processing GO-PCA results,
excluding all visualization commands, which are
`documented in the next section <visualize>`.

.. contents:: Contents
    :depth: 2
    :local:
    :backlinks: none

Inspecting the results: ``gopca_print_info.py``
-----------------------------------------------

In order to simply get a summary of the results contained in a particular
GO-PCA result file, the ``gopca_print_info.py`` command can be used. It prints
things like the number of principal components analyzed, the number of
signatures generated etc. It also outputs a list of all parameter settings
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


Converting the results to MATLAB format: ``gopca_convert_to_matlab.py``
-----------------------------------------------------------------------

This command converts a GO-PCA result file to MATLAB format, using scipy's
`io.savemat` function command from the `scipy` package.

.. argparse::
   :ref: gopca.cli.convert_to_matlab.get_argument_parser
   :prog: gopca_convert_to_matlab.py


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

Combining the signatures from two or more GO-PCA runs: ``gopca_combine_signatures.py`` 
--------------------------------------------------------------------------------------

This command does exactly what the name implies: It combines the signatures
contained in two or more individual GO-PCA result files into a single, new
result file.

.. argparse::
   :ref: gopca.cli.combine_signatures.get_argument_parser
   :prog: gopca_combine_signatures.py
