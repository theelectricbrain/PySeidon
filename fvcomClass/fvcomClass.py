#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
import numpy as np
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr
import netCDF4 as nc
#WB_Alternative: import scipy.io as sio

#Add local path to utilities
sys.path.append('../utilities/')

#Utility import
from shortest_element_path import shortest_element_path

#Local import
from variablesFvcom import _load_var, _load_grid
from functionsFvcom import *
from functionsFvcomThreeD import *
from plotsFvcom import *

class FVCOM:
    '''
Description:
----------
  A class/structure for FVCOM data.
  Functionality structured as follows:
                _Data. = raw netcdf file data
               |_Variables. = fvcom variables and quantities
               |_Grid. = fvcom grid data
               |_QC = Quality Control metadata
    testFvcom._|_Utils. = set of useful functions
               |_Plots. = plotting functions
               |_method_1
               | ...      = methods and analysis techniques intrinsic to fvcom runs
               |_method_n

Inputs:
------
  Takes a file name as input, ex: testFvcom=FVCOM('./path_to_FVOM_output_file/filename')

Options:
-------
    Can be defined for a region only, i.e. a bounding box, as such:
        ax = [min(lon_coord), max(lon_coord), min(lat_coord), max(lat_coord)]
    Can be defined for a time period only, as such:
        tx = ['2012.11.07','2012.11.09'], string of year.month.date
Notes:
-----
    As of right now, only takes a filename as input. It will then load in the
    data (except for timeseries, since loading in the whole time series can be
    too large)
    '''

    def __init__(self, filename, ax=[], tx=[], debug=False):
        ''' Initialize FVCOM class.'''
        self._debug = debug
        if debug:
            print '-Debug mode on-'
        #TR_comments: Add input check and alternative (extract from server)
        #WB_Alternative: self.Data = sio.netcdf.netcdf_file(filename, 'r')
        #WB_comments: scipy has causes some errors, and even though can be
        #             faster, can be unreliable
        self.Data = nc.Dataset(filename, 'r')

        #Metadata
        if hasattr(self.Data, 'QC'):
            self.QC = self.Data.QC
        else:
            text = 'Created from ' + filename
            self.QC = [text]
            # Custom fields
            self.Grid = _load_grid(self.Data, ax, self.QC,
                                   debug=self._debug)
            self.Variables = _load_var(self, tx, self.QC,
                                       debug=self._debug)
            self.Plots = PlotsFvcom(self)
            self.Utils = FunctionsFvcom(self)
            if self.Variables._3D:
                self.UtilsThreeD = FunctionsFvcomThreeD(self)
            if ax:
                self.Grid._ax = ax
            else:
                self.Grid._ax = [min(self.Grid.lon), max(self.Grid.lon),
                                 min(self.Grid.lat), max(self.Grid.lat)]

    def harmonics(self, ind, twodim=True, **kwarg):
        '''Use/Inputs/Outputs of this method has to be clarified !!!'''
        #TR_comments: Add debug flag in Utide: debug=self._debug
        if twodim:
            self.coef = ut_solv(self.Variables.matlabTime,
                                self.Variables.ua[:, ind],
                                self.Variables.va[:, ind],
                                self.Variables.lat[ind],
                                debug=self._debug, **kwarg)
            self.QC.append('ut_solv done for velocity')

        else:
            self.coef = ut_solv(self.Variables.matlabTime,
                                self.Variables.ua[:, ind], [],
                                self.Variables.lat[ind], **kwarg)
            self.QC.append('ut_solv done for elevation')

    def reconstr(self, time):
        '''Use/Inputs/Outputs of this method has to be clarified !!!'''
        #TR_comments: Add debug flag in Utide: debug=self._debug
        if self.coef['aux']['opt']['twodim']:
            self.U, self.V = ut_reconstr(time, self.coef)
            self.QC.append('ut_reconstr done for velocity')
        else:
            self.ts_recon, _ = ut_reconstr(time, self.coef)
            self.QC.append('ut_reconstr done for elevation')

#Test section when running in shell >> python fvcomClass.py
if __name__ == '__main__':

    filename = './test_file/dn_coarse_0001.nc'
    #filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
    test = FVCOM(filename)
    test.harmonics(0, cnstit='auto', notrend=True, nodiagn=True)
    #WB_COMMENTS: fixed matlabttime to matlabtime
    test.reconstr(test.Variables.matlabTime)

    #WB_COMMENTS: This doesn't work with your TR variable convention
    # t = shortest_element_path(test.latc,test.lonc,test.lat,test.lon,test.nv,test.h)
    t = shortest_element_path(test.Grid.latc, test.Grid.lonc,
                              test.Grid.lat, test.Grid.lon,
                              test.Grid.nv, test.Grid.h)

    elements, _ = t.getTargets([[41420, 39763], [48484, 53441],
                                [27241, 24226], [21706, 17458]])



    t.graphGrid()
