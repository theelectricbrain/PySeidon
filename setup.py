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
    setup(name='PySeidon_dvt',
      version='2.1',
      description='Suite of tools for tidal-energy and FVCOM-user communities',
      long_description=readme(),
      url='https://github.com/GrumpyNounours/PySeidon',
      author='Thomas Roc, Jon Smith, Wesley Bowman',
      author_email='thomas.roc@acadiau.ca,wesley.bowman23@gmail.com,'+\
                   'lavieenroux20@gmail.com',
      maintainer='Thomas Roc',
      license='GNU Affero GPL v3.0',
      packages=find_packages(),
      package_dir={'PySeidon_dvt' :'pyseidon_dvt'},
      zip_safe=False)
else:
    setup(name='PySeidon_dvt',
      version='2.1',
      description='Suite of tools for tidal-energy and FVCOM-user communities',
      long_description=readme(),
      url='https://github.com/GrumpyNounours/PySeidon',
      author='Thomas Roc, Jon Smith, Wesley Bowman',
      author_email='thomas.roc@acadiau.ca,wesley.bowman23@gmail.com,'+\
                   'lavieenroux20@gmail.com',
      maintainer='Thomas Roc',
      license='GNU Affero GPL v3.0',
      packages=find_packages(),
      package_dir={'PySeidon_dvt' :'pyseidon_dvt'},
      install_requires=['setuptools', 'utide', 'numpy', 'pandas', 'pydap', 'pydap',
                        'networkx', 'seaborn', 'scipy','matplotlib', 'h5py', 'numexpr',
                        'datetime', 'netCDF4', 'gdal', 'reportlab'],
      zip_safe=False)

