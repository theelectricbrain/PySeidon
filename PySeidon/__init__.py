#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import os
import sys

sys.path.append('./PySeidon/pyseidon/fvcomClass/')
sys.path.append('./PySeidon/pyseidon/adcpClass/')
sys.path.append('./PySeidon/pyseidon/drifterClass/')
sys.path.append('./PySeidon/pyseidon/stationClass/')
sys.path.append('./PySeidon/pyseidon/tidegaugeClass/')
sys.path.append('./PySeidon/pyseidon/validationClass/')
sys.path.append('./PySeidon/pyseidon/utilities/')

#Local import
from fvcomClass import *
from adcpClass import *
from drifterClass import *
from tidegaugeClass import *
from validationClass import *
from stationClass import *

__version__ = '0.0'
__authors__ = ['Wesley Bowman, Thomas Roc, Jonathan Smith']
__licence__ = 'GNU Affero GPL v3.0'
__copyright__ = 'Copyright (c) 2014 EcoEnergyII'

