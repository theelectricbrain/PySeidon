#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from pyseidon_dvt.utilities.miscellaneous import mattime_to_datetime

class _load_tidegauge:
    """
    **'Variables' subset in TideGauge class**

    It contains the following numpy arrays: ::

                           _lat = latitude, float, decimal degrees
                          |_lon = lontitude, float, decimal degrees
     TideGauge.Variables._|_RBR = Raw recording and sampling features
                          |_matlabTime = matlab time, 1D array
                          |_el = sea surface elevation (m) timeserie, 1D array
    """
    def __init__(self, cls, History, debug=False):
        if debug: print 'Loading variables...'

        # Pointer to History
        setattr(self, '_History', History)
        
        self.RBR = cls['RBR']
        self.data = self.RBR.data
        self.matlabTime = self.RBR.date_num_Z
        self.lat = self.RBR.lat
        self.lon = self.RBR.lon
        self.el = self.data - np.mean(self.data)

        #-Append message to History field
        start = mattime_to_datetime(self.matlabTime[0])
        end = mattime_to_datetime(self.matlabTime[-1])
        text = 'Temporal domain from ' + str(start) +\
                ' to ' + str(end)
        self._History.append(text)

        if debug: print '...Done'

