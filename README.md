PySeidon_dvt
================
###Warning###
####This is the development branch of PySeidon and it is subject to daily changes!###

### Project description ###
* This project aims to meet multiple objectives of the [EcoEnergyII](http://tidalenergy.acadiau.ca/EcoEII.html) consortium
  through the setting of a dedicated server and the development of Python
  based packages. This project can be seen as two folded. On the one 
  hand, it aims to enhance data accessibility for all the partners of 
  the [EcoEII](http://tidalenergy.acadiau.ca/EcoEII.html) consortium thanks to simple client protocols. On the other 
  hand, it aims to develop a standardised numerical toolbox gathering 
  specific analysis functions for measured and simulated data (FVCOM model)
  to the [EcoEII](http://tidalenergy.acadiau.ca/EcoEII.html) partners.
* Additionally, this project was the ideal opportunity to transport various
  scripts and packages accumulated over the years into Python. These scripts
  and packages have been extensively used by the tidal energy community for
  more than a decade. The 'Contributors' section of this document is a 
  mere attempt to acknowledge the work of those who participated directly or
  indirectly to the development of this tool box. We are consciously
  standing on the shoulders of a multitude of giants...so please forgive us
  if we forgot one of them.  
* The present package is still a work in progress, so the more feedback,
  the better

### Installation ###
Hydrodynamic model:
* This package has been primarily developed and designed for post-processing FVCOM outputs. One can download FVCOM from [here](http://fvcom.smast.umassd.edu/fvcom/) 

Requirements:
* This package has been designed for Python 2.7: one can download Python from [here](http://www.python.org/download)
* It is also recommended to install Anaconda beforehand: one can download Anaconda from [here](http://continuum.io/downloads#all)
* The HDF5 library is also needed for this package to work: one can download the HDF5 library from [here](https://www.hdfgroup.org/HDF5/)

Dependencies:
Althought they should be automatically resolved during the installation, this package relies on the following dependencies:
* setuptools: One can download setuptools from [here](https://pypi.python.org/pypi/setuptools#installation-instructions)
* UTide: One can download UTide from [here](https://github.com/wesleybowman/UTide)
* Pydap: One can download Pydap from [here](http://www.pydap.org/)
* NetworkX: One can download NetworkX from [here](http://networkx.github.io/documentation/latest/install.html)
* Pandas: One can download Pandas from [here](http://pandas.pydata.org/pandas-docs/stable/install.html)
* Seaborn: One can download Seaborn from [here](http://web.stanford.edu/~mwaskom/software/seaborn/installing.html)
* netCDF4: One can download netCDF4 from [here](https://pypi.python.org/pypi/netCDF4/0.8.2)
* gdal: One can download gdal from [here](https://pypi.python.org/pypi/GDAL/)

Installation:
* Step 1a: Download PySeidon package, save it on your machine and Unzip
* Step 1b: or clone the repository
* Step 2: from a shell, change directory to PySeidon-master folder
* Step 3: from the shell, as superuser/admin, type `python setup.py install`
  or `python setup.py install --user`
* Step 4: choose to automatically resolve (y) or not (n) the dependencies
* Finally, in order to test the installation, type `from pyseidon_dvt import *` in Ipython shell.

Up-dating:
* The code will evolve and improve with time. To up-date, simply "git pull" or download the package
  and go through the installation procedure again.

Recommendations:
* The tutorials and package functioning have been designed for use in IPython shell: One can download IPython from [here](http://ipython.org/)

### Documentation ###
Package's documentation can be found [here](http://grumpynounours.github.io/PySeidon/index.html)

### Contribution guidelines ###
* [Tutorial 0: First steps](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/master/PySeidon_tuto_0.ipynb)
* [Tutorial 1: FVCOM class](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/master/PySeidon_tuto_1.ipynb)
* [Tutorial 2: Station class](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/development/PySeidon_tuto_2.ipynb)
* [Tutorial 3: ADCP class](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/development/PySeidon_tuto_3.ipynb)
* [Tutorial 4: TideGauge class](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/development/PySeidon_tuto_4.ipynb)
* [Tutorial 5: Drifter class](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/development/PySeidon_tuto_5.ipynb)
* [Tutorial 6: Validation class](http://nbviewer.ipython.org/github/GrumpyNounours/PySeidon/blob/development/PySeidon_tuto_6.ipynb)

### Contacts ###
* Project Leader: [Richard Karsten](richard.karsten@acadiau.ca)
* Repository Admin & Software Development Manager: [Thomas Roc](thomas.roc@acadiau.ca)
* Main Developers: [Thomas Roc](thomas.roc@acadiau.ca), [Jonathan Smith](https://github.com/LaVieEnRoux), [Wesley Bowman](https://github.com/wesleybowman), [Kody Crowell](https://github.com/TheKingInYellow)

### Contributors ###
Dr. Richard Karsten, [Aidan Bharath](https://github.com/Aidan-Bharath), Mitchell O'Flaherty-Sproul, Robie Hennigar, [Robert Covill](http://tekmap.ns.ca/), Dr. Joel Culina, Justine McMillan, Dr. Brian Polagye, [Dr. Kristen Thyng](https://github.com/kthyng)...

### Legal Information ###
* Original authorship attributed to Thomas Roc, Wesley Bowman and Jonathan Smith
* Copyright (c) 2014 [EcoEnergyII](http://tidalenergy.acadiau.ca/EcoEII.html)
* Licensed under an Affero GPL style license v3.0 (see License_PySeidon.txt)
