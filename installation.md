--- 
layout: page 
title: "Installation Instructions" 
description: "Installation page" 
--- 
{% include JB/setup %}

pyClusterROI requires [Python](http://www.python.org/) along with the
[NiBabel](http://nipy.sourceforge.net/nibabel/), [NumPy](http://www.numpy.org/)
and [SciPy](http://www.scipy.org/) Python libraries. 

###Python 

The easiest method for installing python on Microsoft Windows or Mac OS
X (and Linux) is to use a Python distribution such as [Enthought
Canopy](https://www.enthought.com),
[Anaconda](https://store.continuum.io/cshop/anaconda/), or
[ActivePython](http://www.activestate.com/activepython).

Python comes preinstalled on many Linux distributions. If it is not available,
it can usually be installed using either:

    sudo apt-get install python

on Debian or Ubuntu systems, or

    sudo yum install python

on Redhat, CentOs, and Fedora systems.

###SciPy and NumPy

Both SciPy and NumPy are bundled with Enthought Canopy, Anaconda, and
ActivePython.  Other methods for installing NumPy and SciPy can be found on the
official [SciPy installation page](http://www.scipy.org/install.html).

On Debian and Ubuntu systems they can be installed using:

    sudo apt-get install python-numpy python-scipy

or on Redhat, CentOs and Fedora systems using:

    sudo yum install numpy scipy

###NiBabel

NiBabel can be installed using the package managers that are bundled with the
Enthought Canopy, Anaconda, or ActivePython. On Linux based distributions
it can be installed from pypi using setuptools.

On Debian and Ubuntu systems this is accomplished by:

    sudo apt-get install python-setuptools
    sudo easy_install nibabel

or on Redhat, Fedora, and CentOs by:

    sudo yum install python-setuptools
    sudo easy_install nibabel

###pyClusterROI

Once all of the pre-requisites are installed, pyClusterROI can be installed by
simply downloading the zip file and uncompressing it into a convenient
directory. The code can be either executed directly from its containing
directory or by adding the directory to the python path. This can be
accomplished by adding the following line to the script that you would like to
call the pyClusterROI libraries from:

    import sys
    sys.path.append("/path/to/pyClusterROI")
    
