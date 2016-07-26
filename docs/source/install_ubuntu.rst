Installation Instructions for Ubuntu Linux
==========================================

In a terminal window:

1. Make sure you have Python 2.7.x or Python 3.5.x installed:
    
    .. code-block:: bash
        
        $ python -V
        Python 3.5.2
    

    Note: Ubuntu 14.04 (trusty) comes with Python version 2.7.6, and
    Ubuntu 16.06 (xenial) comes with Python 3.5.1, so in both cases, you're
    already good to go.

2. Install GO-PCA using `pip`__, the Python package manager
    
    .. code-block:: bash
    
        $ sudo pip install gopca

   Note: If you need to be able to install gopca without admin privileges,
   look into `conda`__ (which you can get by installnig `Miniconda`__) or
   `virtualenv`__ (installable using pip). Both tools allow the creation of
   local Python environments that allow packages to be installed without
   admin privileges.

__ conda_
__ miniconda_
__ virtualenv_

.. _conda: : http://conda.pydata.org/docs/

.. _miniconda: : http://conda.pydata.org/miniconda.html

.. _virtualenv: https://virtualenv.pypa.io/en/stable/