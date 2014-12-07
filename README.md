PySeidon
================

### Project description ###
* This project aims to meet multiple objectives of the EcoEII consortium
  through the setting of a dedicated server and the development of Python
  based packages. This project can be seen as two folded. On the one 
  hand, it aims to enhance data accessibility for all the partners of 
  the EcoEII consortium thanks to simple client protocols. On the other 
  hand, it aims to develop standardised numerical toolbox gathering 
  specific analysis functions for measured and simulated data (FVCOM model)
  to the EcoEII partners.
* Additionally, this project was the ideal opportunity to transport various
  scripts and packages accumulated over the years into Python. These scripts
  and packages have been extensively used by the tidal energy community for
  more than a decade. The 'Contributors' section of this document is a 
  mere attempt to acknowledge the work of those who participated more or
  less indirectly to the development of this tool box. We are consciously
  standing on the shoulders of a multitude giants...so please forgive us
  if we forgot one of them.  
* The present package is still a work in progress, so the more feedback,
  the better

### Installation ###
Hydrodynamic model:
* This package has been primarily developed and designed for post-processing FVCOM outputs. One can download FVCOM from [here](http://fvcom.smast.umassd.edu/fvcom/) 

Requirements:
* This package ha been designed for Python 2.7: One can download Python from [here](http://www.python.org/download)
* It recommended to install Anaconda beforehand: One can download Anaconda from [here](http://continuum.io/downloads#all)

Dependencies:

Althought they should be automatically resolved during the installation, this package relies on the following dependencies:
* setuptools: One can download setuptools from [here](https://pypi.python.org/pypi/setuptools#installation-instructions)
* UTide: One can download UTide from [here](https://github.com/wesleybowman/UTide)
* Pydap: One can download Pydap from [here](http://www.pydap.org/)
* NetworkX: One can download NetworkX from [here](http://networkx.github.io/documentation/latest/install.html)
* Pandas: One can download Pandas from [here](http://pandas.pydata.org/pandas-docs/stable/install.html)
* Seaborn: One can download Seaborn from [here](http://web.stanford.edu/~mwaskom/software/seaborn/installing.html)

Installation:
* Step 1a: Download PySeidon package, save it on your machine and Unzip
* Step 1b: or clone the repository
* Step 2: from a shell, change directory to PySeidon-master folder
* Step 3: from the shell, as superuser/admin, type `python setup.py install`
  or `python setup.py install --user`
* Finally, in order to test the installation, type `from pyseidon import *` in Ipython shell.

Up-dating:
* The code will evolve and improve with time. To up-date, simply go through
  the installation procedure again.

Recommendations:
* The tutorials and package functioning have been designed for use in IPython shell: One can download IPython from [here](http://ipython.org/)


### Contribution guidelines ###
* [Tutorial 0: First steps](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/master/PySeidon_tuto_0.ipynb)
* [Tutorial 1: FVCOM class](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/master/PySeidon_tuto_1.ipynb)

### Contacts ###
* Project Leader: [Richard Karsten](richard.karsten@acadiau.ca)
* Repository Admin & Software Development Manager: [Thomas Roc](thomas.roc@acadiau.ca)
* Main Developers: [Wesley Bowman](https://github.com/wesleybowman), [Thomas Roc](thomas.roc@acadiau.ca), [Jonathan Smith](https://github.com/LaVieEnRoux)

### Contributors ###
Dr. Richard Karsten, [Aidan Bharath](https://github.com/Aidan-Bharath), Mitchell O'Flaherty-Sproul, Robie Hennigar, Dr. Joel Culina, Justine McMillan, Dr. Brian Polagye, [Dr. Kristen Thyng](https://github.com/kthyng)...

### Legal Information ###
* Authorship attributed to Wesley Bowman, Thomas Roc and Jonathan Smith
* Copyright (c) 2014 [EcoEnergyII](http://tidalenergy.acadiau.ca/EcoEII.html)
* Licensed under an Affero GPL style license v3.0 (see License_PySeidon.txt)
