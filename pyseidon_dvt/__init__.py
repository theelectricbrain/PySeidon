#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import os
import sys

local = os.path.dirname(__file__)
sys.path.append(os.path.join(local,'fvcomClass'))
sys.path.append(os.path.join(local,'adcpClass'))
sys.path.append(os.path.join(local,'drifterClass'))
sys.path.append(os.path.join(local,'stationClass'))
sys.path.append(os.path.join(local,'tidegaugeClass'))
sys.path.append(os.path.join(local,'validationClass'))
sys.path.append(os.path.join(local,'utilities'))

#Local import
from utilities import *
from adcpClass import *
from drifterClass import *
from tidegaugeClass import *
from stationClass import *
from fvcomClass import *
from validationClass import *

# Custom error
from pyseidon_dvt.utilities.pyseidon_error import PyseidonError

#Permission info for OpenDap server
#print "OpenDap server connexion info:"

__version__ = '2.0'
__all__ = ["FVCOM", "ADCP", "Drifter", "TideGauge",\
           "Validation", "Station", "utilities", "PyseidonError"]
__authors__ = ['Wesley Bowman, Thomas Roc, Jonathan Smith']
__licence__ = 'GNU Affero GPL v3.0'
__copyright__ = 'Copyright (c) 2014 EcoEnergyII'

