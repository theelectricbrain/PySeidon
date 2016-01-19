#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division
import numpy as np
import scipy.io as sio
import h5py

#Local import
from variablesAdcp import _load_adcp
from functionsAdcp import *
from plotsAdcp import *


class ADCP:
    """
    **A class/structure for ADCP data**

    Functionality structured as follows: ::

                   _Data. = raw matlab file data
                  |_Variables. = useable adcp variables and quantities
                  |_History = Quality Control metadata
        testAdcp._|_Utils. = set of useful functions
                  |_Plots. = plotting functions
                  |_method_1
                  | ...      = methods and analysis techniques intrinsic to ADCPs
                  |_method_n

    Inputs:
      - Only takes a file name as input, ex: testAdcp=ADCP('./path_to_matlab_file/filename')

    *Notes*
      Only handle fully processed ADCP matlab data previously quality-controlled as well
      as formatted through "EnsembleData_FlowFile" matlab script at the mo.

      Throughout the package, the following conventions apply:
      - Coordinates = decimal degrees East and North
      - Directions = in degrees, between -180 and 180 deg., i.e. 0=East, 90=North,
                     +/-180=West, -90=South
      - Depth = 0m is the free surface and depth is negative
    """

    def __init__(self, filename, debug=False):
        """ Initialize ADCP class."""
        self._debug = debug
        self._origin_file = filename
        if debug:
            print '-Debug mode on-' 
        #TR_comments: find a way to dissociate raw and processed data
        self.History = ['Created from ' + filename]
        #TR_comments: *_Raw and *_10minavg open with h5py whereas *_davgBS
        try:
            self.Data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        except NotImplementedError:
            self.Data = h5py.File(filename, 'r')
        #TR_comments: Initialize class structure
        self.Variables = _load_adcp(self, self.History, debug=self._debug)
        self.Plots = PlotsAdcp(self.Variables, debug=self._debug)
        self.Utils = FunctionsAdcp(self.Variables,
                                   self.Plots,
                                   self.History,
                                   debug=self._debug) 

        ##Re-assignement of utility functions as methods
        self.dump_profile_data = self.Plots._dump_profile_data_as_csv

        return