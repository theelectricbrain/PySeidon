#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
from __future__ import division
import numpy as np
import sys
from utide import ut_solv, ut_reconstr
#TR comment: 2 alternatives
#import netCDF4 as nc
from scipy.io import netcdf
from scipy.io import savemat
from scipy.io import loadmat
from pydap.client import open_url
import cPickle as pkl
import copy

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
           |_History = Quality Control metadata
    FVCOM._|_Utils2D. = set of useful functions and methods for 2D and 3D runs
           |_Utils3D. = set of useful functions and methods for 3D runs
           |_Plots. = plotting functions
           |_Save_as = "save as" methods

Inputs:
------
  - filename = path to file, string, 
               ex: testFvcom=FVCOM('./path_to_FVOM_output_file/filename').
               Note that the file can be a pickle file (i.e. *.p)
               or a netcdf file (i.e. *.nc).
               Additionally, either a file path or a OpenDap url could be used. 

Options:
-------
  - ax = defines for a specific spatial region to work with, as such:
             ax = [minimun longitude, maximun longitude,
                 minimun latitude, maximum latitude]
         or use one of the following pre-defined region:
             ax = 'GP', 'PP', 'DG' or 'MP'
         Note that this option permits to extract partial data from the overall file
         and therefore reduce memory and cpu use.

  - tx = defines for a specific temporal period to work with, as such:
             tx = ['2012-11-07T12:00:00','2012.11.09T12:00:00'],
         string of 'yyyy-mm-ddThh:mm:ss'.
         Note that this option permits to extract partial data from the overall file
         and therefore reduce memory and cpu use.

Notes:
-----
  Throughout the package, the following conventions apply:
  - Date = string of 'yyyy-mm-ddThh:mm:ss'
  - Coordinates = decimal degrees East and North
  - Directions = in degrees, between -180 and 180 deg., i.e. 0=East, 90=North,
                 +/-180=West, -90=South
  - Depth = 0m is the free surface and depth is negative
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
            self._origin_file = data['Origin']
            self.History = data['History']
            if debug: print "Turn keys into attributs"
            self.Grid = ObjectFromDict(data['Grid'])
            self.Variables = ObjectFromDict(data['Variables'])
            try:
                if self._origin_file.startswith('http'):
                    #Look for file through OpenDAP server
                    print "Retrieving data through OpenDap server..."
                    self.Data = open_url(data['Origin'])
                    #Create fake attribut to be consistent with the rest of the code
                    self.Data.variables = self.Data
                else:
                    #WB_Alternative: self.Data = sio.netcdf.netcdf_file(filename, 'r')
                    #WB_comments: scipy has causes some errors, and even though can be
                    #             faster, can be unreliable
                    #self.Data = nc.Dataset(data['Origin'], 'r')
                    self.Data = netcdf.netcdf_file(data['Origin'], 'r',mmap=True)
            except: #TR: need to precise the type of error here
                print "the original *.nc file has not been found"
                pass

        #Loading netcdf file         
        elif filename.endswith('.nc'):
            if filename.startswith('http'):
                #Look for file through OpenDAP server
                print "Retrieving data through OpenDap server..."
                self.Data = open_url(filename)
                #Create fake attribut to be consistent with the rest of the code
                self.Data.variables = self.Data
            else:
                #Look for file locally
                print "Retrieving data from " + filename + " ..."
                #WB_Alternative: self.Data = sio.netcdf.netcdf_file(filename, 'r')
                #WB_comments: scipy has causes some errors, and even though can be
                #             faster, can be unreliable
                #self.Data = nc.Dataset(filename, 'r')
                self.Data = netcdf.netcdf_file(filename, 'r',mmap=True)
            text = 'Created from ' + filename
            self._origin_file = filename
            #Metadata
            self.History = [text]
            # Calling sub-class
            print "Initialisation..."
            #print "This might take some time..."
            try:
                self.Grid = _load_grid(self.Data,
                                       ax,
                                       self.History,
                                       debug=self._debug)
                self.Variables = _load_var(self.Data,
                                           self.Grid,
                                           tx,
                                           self.History,
                                           debug=self._debug)
            except MemoryError:
                print '---Data too large for machine memory---'
                print 'Tip: use ax or tx during class initialisation'
                print '---  to use partial data'
                raise

        elif filename.endswith('.mat'):
            print "---Functionality not yet implemented---"
            sys.exit()
        else:
            print "---Wrong file format---"
            sys.exit()

        self.Plots = PlotsFvcom(self.Variables,
                                self.Grid,
                                self._debug)
        self.Util2D = FunctionsFvcom(self.Variables,
                                     self.Grid,
                                     self.Plots,
                                     self.History,
                                     self._debug)

        if self.Variables._3D:
            self.Util3D = FunctionsFvcomThreeD(self.Variables,
                                               self.Grid,
                                               self.Plots,
                                               self.Util2D,
                                               self.History,
                                               self._debug)
            self.Plots.vertical_slice = self.Util3D._vertical_slice

    #Special methods
    def __add__(self, FvcomClass, debug=False):
        """
        This special method permits to stack variables
        of 2 FVCOM objects through a simple addition:
          fvcom1 += fvcom2

        Notes:
        -----
          - fvcom1 and fvcom2 have to cover the exact
            same spatial domain
          - last time step of fvcom1 must be <= to the 
            first time step of fvcom2 
        """
        debug = debug or self._debug
        #Define bounding box
        if debug:
            print "Computing bounding box..."
        if self.Grid._ax == []:
            lon = self.Grid.lon[:]
            lat = self.Grid.lat[:]
            self.Grid._ax = [lon.min(), lon.max(),
                             lat.min(), lat.max()]
        if FvcomClass.Grid._ax == []:
            lon = FvcomClass.Grid.lon[:]
            lat = FvcomClass.Grid.lat[:]
            FvcomClass.Grid._ax = [lon.min(), lon.max(),
                                   lat.min(), lat.max()]
        #series of test before stacking
        if not (self.Grid._ax == FvcomClass.Grid._ax):
            print "---Spatial regions do not match---"
            sys.exit()
        elif not ((self.Grid.nele == FvcomClass.Grid.nele) and
                  (self.Grid.nnode == FvcomClass.Grid.nnode) and
                  (self.Variables._3D == FvcomClass.Variables._3D)):
            print "---Data dimensions do not match---"
            sys.exit()
        else:
            if not (self.Variables.julianTime[-1]<=
                    FvcomClass.Variables.julianTime[0]):
                print "---Data not consecutive in time---"
                sys.exit()
            #Copy self to newself
            newself = copy.copy(self)
            #TR comment: it still points toward self and modifies it
            #            so cannot do fvcom3 = fvcom1 + fvcom2
            if debug:
                print 'Stacking variables...'
            #keyword list for hstack
            kwl=['matlabTime', 'julianTime']
            for key in kwl:
                tmpN = getattr(newself.Variables, key)
                tmpO = getattr(FvcomClass.Variables, key)
                setattr(newself.Variables, key,
                np.hstack((tmpN[:], tmpO[:])))

            #keyword list for vstack
            kwl=['u', 'v', 'w', 'ua', 'va', 'el', 'tke', 'gls']
            for key in kwl:
                try:
                    tmpN = getattr(newself.Variables, key)
                    tmpO = getattr(FvcomClass.Variables, key)
                    setattr(newself.Variables, key,
                    np.vstack((tmpN[:], tmpO[:])))          
                except AttributeError:
                    continue
            #New time dimension
            newself.Grid.ntime = newself.Grid.ntime + FvcomClass.Grid.ntime
            #Append to new object history
            text = 'Data from ' + FvcomClass.History[0].split('/')[-1] \
                 + ' has been stacked'
            newself.History.append(text)

        return newself  
   
    #Methods
    def Save_as(self, filename, fileformat='pickle', debug=False):
        """
        This method saves the current FVCOM structure as:
           - *.p, i.e. python file
           - *.mat, i.e. Matlab file

        Inputs:
        ------
          - filename = path + name of the file to be saved, string

        Keywords:
        --------
          - fileformat = format of the file to be saved, i.e. 'pickle' or 'matlab'
        """
        debug = debug or self._debug
        if debug:
            print 'Saving file...'
        #Define bounding box
        if debug:
            print "Computing bounding box..."
        if self.Grid._ax == []:
            lon = self.Grid.lon[:]
            lat = self.Grid.lat[:]
            self.Grid._ax = [lon.min(), lon.max(),
                             lat.min(), lat.max()]
        #Save as different formats
        if fileformat=='pickle':
            filename = filename + ".p"
            f = open(filename, "wb")
            data = {}
            data['Origin'] = self._origin_file
            data['History'] = self.History
            data['Grid'] = self.Grid.__dict__
            data['Variables'] = self.Variables.__dict__
            #TR: Force caching Variables otherwise error during loading
            #    with 'netcdf4.Variable' type (see above)
            for key in data['Variables']:
                listkeys=['Variable', 'ArrayProxy', 'BaseType'] 
                if any([type(data['Variables'][key]).__name__==x for x in listkeys]):
                    if debug:
                        print "Force caching for " + key
                    data['Variables'][key] = data['Variables'][key][:]
            #Unpickleable objects
            data['Grid'].pop("triangle", None)
            #TR: Force caching Variables otherwise error during loading
            #    with 'netcdf4.Variable' type (see above)
            for key in data['Grid']:
                listkeys=['Variable', 'ArrayProxy', 'BaseType'] 
                if any([type(data['Grid'][key]).__name__==x for x in listkeys]):
                    if debug:
                        print "Force caching for " + key
                    data['Grid'][key] = data['Grid'][key][:]
            #Save in pickle file
            if debug:
                print 'Dumping in pickle file...'
            try:    
                pkl.dump(data, f, protocol=pkl.HIGHEST_PROTOCOL)
            except SystemError:
                try:
                    print "---Very large data, this may take a while---"
                    pkl.dump(data, f)
                except SystemError:                 
                    print "---Data too large for machine memory---"
                    print "Tip: use ax or tx during class initialisation"
                    print "---  to use partial data"
                    sys.exit()
           
            f.close()
        elif fileformat=='matlab':
            filename = filename + ".mat"
            #TR comment: based on MitchellO'Flaherty-Sproul's code
            dtype = float
            data = {}
            Grd = {}
            Var = {}
            data['Origin'] = self._origin_file
            data['History'] = self.History
            Grd = self.Grid.__dict__
            Var = self.Variables.__dict__
            #TR: Force caching Variables otherwise error during loading
            #    with 'netcdf4.Variable' type (see above)
            for key in Var:
                listkeys=['Variable', 'ArrayProxy', 'BaseType'] 
                if any([type(Var[key]).__name__==x for x in listkeys]):
                    if debug:
                        print "Force caching for " + key
                    Var[key] = Var[key][:]
                #keyV = key + '-var'
                #data[keyV] = Var[key]
                data[key] = Var[key]
            #Unpickleable objects
            Grd.pop("triangle", None)
            for key in Grd:
                listkeys=['Variable', 'ArrayProxy', 'BaseType'] 
                if any([type(Grd[key]).__name__==x for x in listkeys]):
                    if debug:
                        print "Force caching for " + key
                    Grd[key] = Grd[key][:]
                #keyG = key + '-grd'
                #data[keyG] = Grd[key]
                data[key] = Grd[key]

            #Save in mat file file
            if debug:
                print 'Dumping in matlab file...'
            savemat(filename, data, oned_as='column')       
        else:
            print "---Wrong file format---"                            

#Test section when running in shell >> python fvcomClass.py
#if __name__ == '__main__':

    #filename = './test_file/dn_coarse_0001.nc'
    #test = FVCOM(filename)
    #test.harmonics(0, cnstit='auto', notrend=True, nodiagn=True)
    #WB_COMMENTS: fixed matlabttime to matlabtime
    #test.reconstr(test.Variables.matlabTime)
