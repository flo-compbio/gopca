Installation Instructions for Ubuntu Linux
==========================================

In a terminal window:

1. Make sure you have Python 2.7 installed:
    
    .. code-block:: bash
        
        $ python2.7 -V
        Python 2.7.6
    

    Ubuntu 14.04 (trusty) currently has Python version 2.7.6, but any Python 2.7 release should work.

2. Install GO-PCA dependencies: (pip, Cython, NumPy, scikit-learn, Matplotlib, and sphinx for building this documentation locally):
    
    .. code-block:: bash
    
        $ sudo apt-get install python-pip cython python-numpy python-matplotlib python-scikits-learn ipython ipython-notebook

.. "3. Make sure the Ubuntu package python-sphinx is *not* installed:
    
    .. code-block:: bash
    
        $ sudo apt-get remove python-sphinx
    
    (The reason the package needs to be uninstalled is that it is an older version that conflicts with the version required by GO-PCA.)

4. Install the xlmhg python package using pip:

    .. code-block:: bash
    
        $ sudo pip install xlmhg

5. Install GO-PCA using pip (including all remaining dependencies):
    
    .. code-block:: bash
    
        $ sudo pip install gopca
