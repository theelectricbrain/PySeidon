#!/usr/bin/python2.7
# encoding: utf-8
from setuptools import setup, find_packages
from numpy.distutils.misc_util import Configuration

def readme():
    with open('README.md') as f:
        return f.read()

option = raw_input("Resolve dependencies (y) or (n): ")
option = option.lower()

if option=='n':
    setup(name='PySeidon',
      version='1.5',
      description='Suite of tools for tidal-energy and FVCOM-user communities',
      long_description=readme(),
      url='https://github.com/GrumpyNounours/PySeidon',
      author='Wesley Bowman, Thomas Roc, Jon Smith',
      author_email='thomas.roc@acadiau.ca,wesley.bowman23@gmail.com,'+\
                   'lavieenroux20@gmail.com',
      maintainer='Thomas Roc',
      license='GNU Affero GPL v3.0',
      packages=find_packages(),
      package_dir={'PySeidon' :'pyseidon'},
      zip_safe=False)
else:
    setup(name='PySeidon',
      version='2.0',
      description='Suite of tools for tidal-energy and FVCOM-user communities',
      long_description=readme(),
      url='https://github.com/GrumpyNounours/PySeidon',
      author='Wesley Bowman, Thomas Roc, Jon Smith',
      author_email='thomas.roc@acadiau.ca,wesley.bowman23@gmail.com,'+\
                   'lavieenroux20@gmail.com',
      maintainer='Thomas Roc',
      license='GNU Affero GPL v3.0',
      packages=find_packages(),
      package_dir={'PySeidon' :'pyseidon'},
      install_requires=['setuptools', 'utide', 'numpy', 'pandas', 'pydap', 'pydap',
                        'networkx', 'seaborn', 'scipy','matplotlib', 'h5py', 'numexpr',
                        'datetime', 'netCDF4'],
      zip_safe=False)

