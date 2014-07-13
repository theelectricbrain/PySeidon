#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np
import sys
#TR_Alternative: import scipy.io as sio
import h5py

#Add local path to utilities
sys.path.append('../utilities/')

#Local import
from variablesAdcp import _load_adcp
from functionsAdcp import *
from plotsAdcp import *


class ADCP:
    ''' 
Description:
-----------
  A class/structure for ADCP data.
  Functionality structured as follows:
               _Data. = raw matlab file data
              |_Variables. = useable adcp variables and quantities
              |_QC = Quality Control metadata
    testAdcp._|_Utils. = set of useful functions
              |_Plots. = plotting functions
              |_method_1
              | ...      = methods and analysis techniques intrinsic to ADCPs
              |_method_n

Inputs:
------
  Only takes a file name as input, ex: testAdcp=ADCP('./path_to_matlab_file/filename')
    '''

    def __init__(self, filename, debug=False):
        ''' Initialize ADCP class.
            Notes: only handle raw ADCP matlab data at the mo.'''    
        self._debug = debug
        if debug:
            print '-Debug mode on-' 

        self.QC = ['Raw data']
        self.Data = h5py.File(filename)
        #TR_comments: *_Raw and *_10minavg open with h5py whereas *_davgBS
        #TR_Alternative: self.Data = sio.loadmat(filename,struct_as_record=False, squeeze_me=True)
        self.Variables = _load_adcp(self.Data['adcp'], debug=self._debug)
        self.Utils = FunctionsAdcp(self)
        self.Plots = PlotsAdcp(self)   


if __name__ == '__main__':
    filename = '../test_files/adcp/GP-120726-BPd_raw.mat'
    data = ADCP(filename, debug=True)
