#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
import numpy as np
import sys
from utide import ut_solv, ut_reconstr
import scipy.io as sio
import netCDF4 as nc
from variables import _load_var, _load_grid

#Add local path to utilities
sys.path.append('../utilities/')

#Local import
from shortest_element_path import shortest_element_path
from functions import *
from plots import *

class FVCOM:
    '''
    A class for FVCOM data.
    As of right now, only takes a filename as input. It will then load in the
    data (except for timeseries, since loading in the whole time series can be
    too large)

    ax can be defined as a region, i.e. a bounding box.
    An example:
        ax = [min(lon_coord), max(lon_coord), min(lat_coord), max(lat_coord)]
    '''

    def __init__(self, filename, ax=[]):
        self.QC = ['raw data']

        # Add input check and alternative (extract from server)
        self.Data = nc.Dataset(filename, 'r')
        #self.Data = sio.netcdf.netcdf_file(filename, 'r')

        self.Variables = _load_var(self.Data)
        self.Grid = _load_grid(self.Data)

        # Invisible parameter
        if ax:
            self._ax = ax
        else:
            self._ax = [min(self.Variables.lon), max(self.Variables.lon), min(self.Variables.lat), max(self.Variables.lat)]

        self.Utils = Functions(self)
        self.Plots = Plots(self)

    def harmonics(self, ind, twodim=True, **kwarg):

        if twodim:
            self.coef = ut_solv(self.Variables.matlabTime, self.Variables.ua[:, ind],
                                self.Variables.va[:, ind], self.Variables.lat[ind], **kwarg)
            self.QC.append('ut_solv done for velocity')

        else:	
            self.coef = ut_solv(self.Variables.matlabTime, self.Variables.ua[:, ind], [],
                                self.Variables.lat[ind], **kwarg)
            self.QC.append('ut_solv done for elevation')

    def reconstr(self, time):
        if self.coef['aux']['opt']['twodim']:
            self.U, self.V = ut_reconstr(time, self.coef)
            self.QC.append('ut_reconstr done for velocity')
        else:
            self.ts_recon, _ = ut_reconstr(time, self.coef)
            self.QC.append('ut_reconstr done for elevation')

if __name__ == '__main__':

    filename = './test_file/dn_coarse_0001.nc'
    test = FVCOM(filename)
    test.harmonics(0, cnstit='auto', notrend=True, nodiagn=True)
    test.reconstr(test.time)

    t = shortest_element_path(test.latc,test.lonc,test.lat,test.lon,test.nv,test.h)
    elements, _ = t.getTargets([[41420,39763],[48484,53441],
                                [27241,24226],[21706,17458]])



    # t.graphGrid()
