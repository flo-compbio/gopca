Installation Instructions for Microsoft Windows
===============================================

1. `Download`__ and install the "Anaconda" Python distribution (if you haven't already). This distribution includes all third-party packages required by GO-PCA.

__ anaconda_
.. _anaconda: http://continuum.io/downloads

2. `Download <http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_ and install the Microsoft Visual C++ Compiler for Python 2.7. This compiler is required for installing the xlmhg python package (Step 3).

Then, on the command line (Start -> Run -> "cmd"):

3. Install the xlmhg python package using pip:
    
    .. code-block:: bat
    
        > pip install xlmhg

4. Install GO-PCA using pip (including all remaining dependencies):
    
    .. code-block:: bat
    
        > pip install gopca
