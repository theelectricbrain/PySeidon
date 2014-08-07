#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
from __future__ import division
import numpy as np
import sys
from utide import ut_solv, ut_reconstr
#import netCDF4 as nc
from scipy.io import netcdf
import cPickle as pkl
#import pickle as pkl
#WB_Alternative: import scipy.io as sio

#Add local path to utilities
sys.path.append('../utilities/')

#Utility import
from shortest_element_path import shortest_element_path
from object_from_dict import ObjectFromDict

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
    testFvcom._|_Utils2D. = set of useful functions
               |_Utils3D. = set of useful functions
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
        ax = [minimun longitude, maximun longitude,
              minimun latitude, maximum latitude]
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

        #Loading pickle file
        if filename.endswith('.p'):
            f = open(filename, "rb")
            data = pkl.load(f)
            self.QC = data['QC']
            self.Grid = ObjectFromDict(data['Grid'])
            self.Variables = ObjectFromDict(data['Variables'])
            try:
                #self.Data = nc.Dataset(data['Origin'], 'r')
                self.Data = netcdf.netcdf_file(data['Origin'], 'r',mmap=True)
            except: #TR: need to precise the type of error here
                print "the original *.nc file has not been found"
                pass

        #Loading netcdf file         
        elif filename.endswith('.nc'):
            #WB_Alternative: self.Data = sio.netcdf.netcdf_file(filename, 'r')
            #WB_comments: scipy has causes some errors, and even though can be
            #             faster, can be unreliable
            #self.Data = nc.Dataset(filename, 'r')
            self.Data = netcdf.netcdf_file(filename, 'r',mmap=True)
            text = 'Created from ' + filename
            self._origin_file = filename
            #Metadata
            self.QC = [text]
            # Calling sub-class
            print "Initialisation..."
            #print "This might take some time..."
            try:
                self.Grid = _load_grid(self.Data,
                                       ax,
                                       self.QC,
                                       debug=self._debug)
                self.Variables = _load_var(self.Data,
                                           self.Grid,
                                           tx,
                                           self.QC,
                                           debug=self._debug)
            except MemoryError:
                print '---Data too large for machine memory---'
                print 'Tip: use ax or tx during class initialisation'
                print '---  to use partial data'
                raise

        elif filename.endswith('.mat'):
            print "Cannot handle *.mat file yet"

        else:
            print "---Wrong file format---"
            sys.exit()

        self.Plots = PlotsFvcom(self.Variables,
                                self.Grid,
                                self._debug)
        self.Util2D = FunctionsFvcom(self.Variables,
                                     self.Grid,
                                     self.Plots,
                                     self.QC,
                                     self._debug)

        if self.Variables._3D:
            self.Util3D = FunctionsFvcomThreeD(self.Variables,
                                               self.Grid,
                                               self.Plots,
                                               self.Util2D,
                                               self.QC,
                                               self._debug)
            self.Plots.vertical_slice = self.Util3D._vertical_slice

    def Save_as(self, filename, fileformat='pickle', debug=False):
        """
        Save the current FVCOM structure as:
           - *.p, i.e. python file
           - *.mat, i.e. Matlab file
        Inputs:
        ------
          - filename = name of the file to be saved, string
        Keywords:
        --------
          - fileformat = format of the file to be saved, i.e. 'pickle' or 'matlab'
        """
        debug = debug or self._debug
        if debug:
            print 'Saving file...'
        #TR: to be developed
        if fileformat=='pickle':
            filename = filename + ".p"
            f = open(filename, "wb")
            data = {}
            data['Origin'] = self._origin_file
            data['QC'] = self.QC
            data['Grid'] = self.Grid.__dict__
            data['Variables'] = self.Variables.__dict__
            #TR: Force caching Variables otherwise error during loading
            #    with 'netcdf4.Variable' type (see above)
            for key in data['Variables']:
                if type(data['Variables'][key]).__name__=='Variable':
                    if debug:
                        print "Force caching for " + key
                    data['Variables'][key] = data['Variables'][key][:]
            #Save in pickle file
            if debug:
                print 'Dumping in pickle file...'
            try:    
                pkl.dump(data, f, protocol=pkl.HIGHEST_PROTOCOL)
            except SystemError:
                print '---Data too large for machine memory---'
                print 'Tip: use ax or tx during class initialisation'
                print '---  to use partial data'
                raise
           
            f.close()
        elif fileformat=='matlab':
            filename = filename + ".mat"
            print "Functionality not yet implemented"
        else:
            print "---Wrong file format---"

    def Harmonic_analysis(self, ind, twodim=True, **kwarg):
        '''
        Description:
        ----------
        Harmonic_analysis calls ut_solv. Depending on whether the user wants velocity
        or elevation, it will call the correct version of ut_solv based on the
        twodim option.

        Inputs:
        ------
        - Takes an index so that ut_solv can be run on the data and twodim.
        - twodim=True means that ut_solv will be done for velocity, and
        - twodim=False means that ut_solv will be done for elevation.

        Options:
        -------
        - Options are the same as for ut_solv, which are shown below with
          their default values:
            conf_int=True; cnstit='auto'; notrend=0; prefilt=[]; nodsatlint=0;
            nodsatnone=0; gwchlint=0; gwchnone=0; infer=[]; inferaprx=0;
            rmin=1; method='cauchy'; tunrdn=1; linci=0; white=0; nrlzn=200;
            lsfrqosmp=1; nodiagn=0; diagnplots=0; diagnminsnr=2;
            ordercnstit=[]; runtimedisp='yyy'

        Notes:
        -----
        For more detailed information about ut_solv, please see
        https://github.com/wesleybowman/UTide

        '''
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
                                self.Variables.el[:, ind], [],
                                self.Variables.lat[ind], **kwarg)
            self.QC.append('ut_solv done for elevation')

    def Harmonic_reconstruction(self, time):
        '''
        Description:
        ----------
        Harmonic_reconstruction calls ut_reconstr. This function assumes harmonics (ut_solv)
        has already been executed. If it has not, it will inform the user of
        the error and ask them to run harmonics. It asks the user to run it
        since it needs an index at which to run, and there isn't a default
        index.

        Inputs:
        ------
        Takes a time series for ut_reconstr to do the reconstruction to.

        Options:
        -------
        Options are the same as for ut_reconstr, which are shown below with
        their default values:
        cnstit = []
        minsnr = 2
        minpe = 0

        Notes:
        -----
        For more detailed information about ut_reconstr, please see
        https://github.com/wesleybowman/UTide

        '''
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
    t = shortest_element_path(test.Grid.latc, test.Grid.lonc,
                              test.Grid.lat, test.Grid.lon,
                              test.Grid.trinodes, test.Grid.h)

    elements, _ = t.getTargets([[41420, 39763], [48484, 53441],
                                [27241, 24226], [21706, 17458]])



    t.graphGrid()
