#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import os
import sys

sys.path.append('./PySeidon/fvcomClass/')
sys.path.append('./PySeidon/adcpClass/')
sys.path.append('./PySeidon/drifterClass/')
sys.path.append('./PySeidon/stationClass/')
sys.path.append('./PySeidon/tidegaugeClass/')
sys.path.append('./PySeidon/validationClass/')
sys.path.append('./PySeidon/utilities/')

#Local import
from fvcomClass import *
from adcpClass import *
#from drifterClass import *
from tidegaugeClass import *
from validationClass import *
from stationClass import *

__version__ = '0.0'
__authors__ = ['Thomas Roc, Wesley Bowman, Jonathan Smith']
__licence__ = 'GNU Affero GPL v3.0'
__copyright__ = 'Copyright (c) 2014 EcoEnergyII'

