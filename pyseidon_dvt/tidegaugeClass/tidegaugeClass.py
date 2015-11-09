#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

import scipy.io as sio

#Local import
from variablesTidegauge import _load_tidegauge
from functionsTidegauge import *
from plotsTidegauge import *

class TideGauge:
    """
    **A class/structure for tide gauge data**

    Functionality structured as follows: ::

                    _Data. = raw matlab file data
                   |_Variables. = useable tide gauge variables and quantities
        TideGauge._|_History = Quality Control metadata
                   |_Utils. = set of useful functions
                   |_Plots. = plotting functions

    Inputs:
      - Only takes a file name as input, ex: testTG=TideGauge('./path_to_matlab_file/filename')

    *Notes*
      - Only handle fully processed tide gauge matlab data at the mo.
      Throughout the package, the following conventions apply:
      - Coordinates = decimal degrees East and North
      - Directions = in degrees, between -180 and 180 deg., i.e. 0=East, 90=North,
                     +/-180=West, -90=South
      - Depth = 0m is the free surface and depth is negative
    """
    def __init__(self, filename, debug=False):
        self._debug = debug
        if debug: print '-Debug mode on-'
        if debug: print 'Loading...'
        #Metadata
        self._origin_file = filename
        self.History = ['Created from ' + filename]

        self.Data = sio.loadmat(filename,
                               struct_as_record=False, squeeze_me=True)
        self.Variables = _load_tidegauge(self.Data, self.History, debug=self._debug)

        self.Plots = PlotsTidegauge(self.Variables, debug=self._debug)

        self.Utils = FunctionsTidegauge(self.Variables,
                                        self.Plots,
                                        self.History,
                                        debug=self._debug)

        ##Re-assignement of utility functions as methods
        self.dump_profile_data = self.Plots._dump_profile_data_as_csv

