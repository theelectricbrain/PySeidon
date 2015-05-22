#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from numpy.ma import MaskError
import h5py
from miscellaneous import *

class _load_drifter:
    """
'Variables' subset in Tidegauge class contains the following numpy arrays:
    """
    def __init__(self,cls, debug=False):
        if debug:
            print 'Loading variables...'
        self.matlabTime = cls.Data['velocity'].vel_time[:]
        #Sorting values with increasing time step
        sortedInd = self.matlabTime.argsort()
        self.matlabTime.sort()
        self.lat = cls.Data['velocity'].vel_lat[sortedInd]
        self.lon = cls.Data['velocity'].vel_lon[sortedInd]
        self.u = cls.Data['velocity'].u[sortedInd]
        self.v = cls.Data['velocity'].v[sortedInd]       
        if debug:
            print '...Passed'
