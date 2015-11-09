#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

import scipy.io as sio
import h5py

#Local import
from variablesDrifter import _load_drifter
# from functionsDrifter import FunctionsDrifter
from plotsDrifter import *


class Drifter:
    """
    **A class/structure for Drifter data**

    Functionality structured as follows: ::

                   _Data. = raw matlab file data
                  |_Variables. = useable drifter variables and quantities
                  |_History = Quality Control metadata
        testAdcp._|_Utils. = set of useful functions
                  |_Plots. = plotting functions
                  |_method_1
                  | ...      = methods and analysis techniques intrinsic to drifters
                  |_method_n

    Inputs:
      - Only takes a file name as input, ex: testDrifter=Drifter('./path_to_matlab_file/filename')

    Notes:
      Only handle fully processed drifter matlab data previously quality-controlled at the mo.

      Throughout the package, the following conventions apply:
      - Coordinates = decimal degrees East and North
      - Directions = in degrees, between -180 and 180 deg., i.e. 0=East, 90=North,
                     +/-180=West, -90=South
      - Depth = 0m is the free surface and depth is negative
    """
    def __init__(self, filename, debug=False):
        """ Initialize Drifter class.
            Notes: only handle processed Drifter matlab data at the mo."""
        self._debug = debug
        self._origin_file = filename
        if debug:
            print '-Debug mode on-'
        #Load data from Matlab file 
        #TR_comments: find a way to dissociate raw and processed data
        #TR_comments: *_Raw and *_10minavg open with h5py whereas *_davgBS
        try:
            self.Data = sio.loadmat(filename,struct_as_record=False, squeeze_me=True)
        except NotImplementedError:
            self.Data = h5py.File(filename)

        #Store info in "History" field
        self.History = ['Created from ' + filename]

        # KC comment: for some reason, some drifter MATLAB files 
        # have 'Comments' as a key in the variables structure,
        # while others have 'comments' as a key.
        if debug: print '-adding comments to History-'

        if 'Comments' in self.Data:
            for comment in self.Data['Comments'][:]:
                self.History.append(str(comment))
        elif 'comments' in self.Data:
            for comment in self.Data['comments'][:]:
                self.History.append(str(comment))
        elif debug:
            print '-no comments found-'


        #Initialize class structure
        self.Variables = _load_drifter(self, self.History, debug=self._debug)
        self.Plots = PlotsDrifter(self.Variables, debug=self._debug)
        #self.Utils = FunctionsAdcp(self.Variables,
        #                           self.Plots,
        #                           self.History,
        #                           debug=self._debug) 

        ##Re-assignement of utility functions as methods
        self.dump_data_as_csv = self.Plots._dump_data_as_csv    

        return