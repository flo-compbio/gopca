Installation Instructions for Ubuntu Linux
==========================================

In a terminal window:

1. Make sure you have Python 2.7.x or Python 3.5.x installed:
    
    .. code-block:: bash
        
        $ python -V
        Python 3.5.1

    Note: Ubuntu 12.04 (precise) and Ubuntu 14.04 (trusty) both come with
    Python 2.7.x, and Ubuntu 16.04 (xenial) comes with Python 3.5.1.
    Therefore, if you're running any of those Ubuntu versions, you should
    already be good to go.

2. Install GO-PCA using `pip`__, the Python package manager
    
    .. code-block:: bash
    
        $ sudo pip install gopca

    Note: If you need to be able to install gopca without admin privileges,
    look into `conda`__ (which you can get by installing `Miniconda`__) or
    `virtualenv`__ (installable using pip). Both tools allow the creation of
    local Python environments that allow packages to be installed without
    admin privileges.


__ pip_
__ conda_
__ miniconda_
__ virtualenv_

.. _pip: https://pip.pypa.io/en/stable/

.. _conda: : http://conda.pydata.org/docs/

.. _miniconda: : http://conda.pydata.org/miniconda.html

.. _virtualenv: https://virtualenv.pypa.io/en/stable/