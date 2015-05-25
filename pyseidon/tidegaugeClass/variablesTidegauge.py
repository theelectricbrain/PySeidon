#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np

class _load_tidegauge:
    """
'Variables' subset in TideGauge class contains the following numpy arrays:

                       _lat = latitude, float, decimal degrees
                      |_lon = lontitude, float, decimal degrees
 TideGauge.Variables._|_RBR = Raw recording and sampling features 
                      |_matlabTime = matlab time, 1D array
                      |_el = sea surface elevation (m) timeserie, 1D array 
    """
    def __init__(self, cls, debug=False):
        if debug: print 'Loading variables...'
        self.RBR = cls['RBR']
        data = self.RBR.data
        self.matlabTime = self.RBR.date_num_Z
        self.lat = self.RBR.lat
        self.lon = self.RBR.lon
        self.el = data - np.mean(data)

        if debug: print '...Done'

