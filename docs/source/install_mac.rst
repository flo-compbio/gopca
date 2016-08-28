Installation Instructions for Mac OS X
======================================

1. `Download`__ and install Python 3.5.

2. In a terminal window, run:

    .. code-block:: bash

        $ pip3 install gopca

  This will install the ``go-pca.py`` script into the following directory::

    /Library/Frameworks/Python.framework/Versions/3.5/bin/


3. To test if the installation was successful, try running ``go-pca.py``:

    .. code-block:: bash

        $ /Library/Frameworks/Python.framework/Versions/3.5/bin/go-pca.py --version

4. (optional) If you want to avoid having to type the full path everytime
   you want to run GO-PCA, add the following two lines to your
   ``~/.bash_profile`` file::

     PATH="/Library/Frameworks/Python.framework/Versions/3.5/bin:${PATH}"
     export PATH


__ download_python_

.. _download_python: https://www.python.org/downloads/