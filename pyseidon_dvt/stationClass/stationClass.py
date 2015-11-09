#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

import numpy as np
import netCDF4 as nc
from scipy.io import netcdf
from scipy.io import savemat
from pydap.client import open_url
import cPickle as pkl
import copy

#Utility import
from pyseidon_dvt.utilities.object_from_dict import ObjectFromDict
from pyseidon_dvt.utilities.miscellaneous import findFiles

# Custom error
from pyseidon_dvt.utilities.pyseidon_error import PyseidonError

#Local import
from variablesStation import _load_var, _load_grid
from functionsStation import *
from functionsStationThreeD import *
from plotsStation import *

class Station:
    """
    **A class/structure for Station data**

    Functionality structured as follows: ::

                  _Data. = raw netcdf file data
                 |_Variables. = fvcom station variables and quantities
                 |_Grid. = fvcom station grid data
                 |_History = Quality Control metadata
        Station._|_Utils2D. = set of useful functions for 2D and 3D station
                 |_Utils3D. = set of useful functions for 3D station
                 |_Plots. = plotting functions
                 |_Harmonic_analysis = harmonic analysis based UTide package
                 |_Harmonic_reconstruction = harmonic reconstruction based UTide package

    Inputs:
      - filename = path to netcdf file or folder, string,
                   ex: testFvcom=Station('./path_to_FVOM_output_file/filename')
                       testFvcom=Station('./path_to_FVOM_output_file/folder/')

        Note that if the path point to a folder all the similar netCDF station files
        will be stack together.
        Note that the file can be a pickle file (i.e. *.p) or a netcdf file (i.e. *.nc).

    Options:
      - elements = indices to extract, list of integers

    *Notes*
      Throughout the package, the following conventions apply:
      - Date = string of 'yyyy-mm-dd hh:mm:ss'
      - Coordinates = decimal degrees East and North
      - Directions = in degrees, between -180 and 180 deg., i.e. 0=East, 90=North,
                     +/-180=West, -90=South
      - Depth = 0m is the free surface and depth is negative

    """
    def __init__(self, filename, elements=slice(None), debug=False):
        #Class attributs
        self._debug = debug
        self._isMulti(filename)
        if not self._multi:
            self._load(filename, elements, debug=debug )
            self.Plots = PlotsStation(self.Variables,
                                      self.Grid,
                                      self._debug)
            self.Util2D = FunctionsStation(self.Variables,
                                           self.Grid,
                                           self.Plots,
                                           self.History,
                                           self._debug)
            if self.Variables._3D:
                self.Util3D = FunctionsStationThreeD(
                                       self.Variables,
                                       self.Grid,
                                       self.Plots,
                                       self.History,
                                       self._debug) 
        else:
            print "---Finding matching files---"
            self._matches = findFiles(filename, 'STATION')
            filename = self._matches.pop(0)
            self._load(filename, elements, debug=debug )
            self.Plots = PlotsStation(self.Variables,
                                      self.Grid,
                                      self._debug)
            self.Util2D = FunctionsStation(self.Variables,
                                           self.Grid,
                                           self.Plots,
                                           self.History,
                                           self._debug)
            if self.Variables._3D:
                self.Util3D = FunctionsStationThreeD(
                                       self.Variables,
                                       self.Grid,
                                       self.Plots,
                                       self.History,
                                       self._debug) 
            for entry in self._matches:
                #Define new 
                text = 'Created from ' + entry
                tmp = {}
                tmp['Data'] = self._load_nc(entry)
                tmp['History'] = [text]
                tmp['Grid'] = _load_grid(tmp['Data'], elements, [], debug=self._debug)
                tmp['Variables'] = _load_var(tmp['Data'], elements, tmp['Grid'], [],
                                             debug=self._debug)
                tmp = ObjectFromDict(tmp)
                self = self.__add__(tmp)

        ##Re-assignement of utility functions as methods
        self.dump_profile_data = self.Plots._dump_profile_data_as_csv

        return

    def _isMulti(self, filename):
        """Tells if filename point to a file or a folder"""
        split = filename.split('/')
        if split[-1]:
            self._multi = False
        else:
            self._multi = True

    def _load(self, filename, elements, debug=False):
        """Loads data from *.nc, *.p and OpenDap url"""
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
                    self.Data = self._load_nc(data['Origin'])
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
                self.Data = self._load_nc(filename)
            #Metadata
            text = 'Created from ' + filename
            self._origin_file = filename
            self.History = [text]
            # Calling sub-class
            print "Initialisation..."
            try:
                self.Grid = _load_grid(self.Data,
                                       elements,
                                       self.History,
                                       debug=self._debug)
                self.Variables = _load_var(self.Data,
                                           elements,
                                           self.Grid,
                                           self.History,
                                           debug=self._debug)

            except MemoryError:
                print '---Data too large for machine memory---'
                print 'Tip: use ax or tx during class initialisation'
                print '---  to use partial data'
                raise

        elif filename.endswith('.mat'):
            raise PyseidonError("---Functionality not yet implemented---")
        else:
            raise PyseidonError("---Wrong file format---")

    def _load_nc(self, filename):
        """loads netcdf file"""
        #Look for file locally
        #print "Retrieving data from " + filename + " ..."
        # WB_Alternative: self.Data = sio.netcdf.netcdf_file(filename, 'r')
        # WB_comments: scipy has causes some errors, and even though can be
        #             faster, can be unreliable
        try:
            Data = netcdf.netcdf_file(filename, 'r', mmap=True)
        except ValueError: #TR: quick fix due to mmap
            Data = nc.Dataset(filename, 'r')
        return Data

    #Special methods
    def __add__(self, StationClass, debug=False):
        """
        This special method permit to stack variables
        of 2 Station objects through a simple addition: ::

          station1 += station2

        *Notes*
          - station1 and station2 have to cover the exact
            same spatial domain
          - last time step of station1 must be <= to the 
            first time step of station2 
        """
        debug = debug or self._debug
        if debug: print "Find matching elements..."
        #Find matching elements
        origNele = self.Grid.nele
        origEle = []
        #origName = self.Grid.name
        origX = self.Grid.x[:]
        origY = self.Grid.y[:]
        newNele = StationClass.Grid.nele
        newEle = []
        #newName = StationClass.Grid.name
        newX = StationClass.Grid.x[:]
        newY = StationClass.Grid.y[:]
        for i in range(origNele):
            for j in range(newNele):
                #Match based on names
                #if (all(origName[i,:]==newName[j,:])):
                #    origEle.append(i)
                #    newEle.append(j)
                #Match based on coordinates
                if ((origX[i]==newX[j]) and (origY[i]==newY[j])):
                    origEle.append(i)
                    newEle.append(j)
                
        print len(origEle), " points will be stacked..."

        if len(origEle)==0:
            raise PyseidonError("---No matching element found---")
        elif not (self.Variables._3D == StationClass.Variables._3D):
            raise PyseidonError("---Data dimensions do not match---")
        else:
            if not (self.Variables.julianTime[-1]<=
                    StationClass.Variables.julianTime[0]):
                raise PyseidonError("---Data not consecutive in time---")
            #Copy self to newself
            newself = copy.copy(self)
            #TR comment: it still points toward self and modifies it
            #            so cannot do Station3 = Station1 + Station2
            if debug:
                print 'Stacking variables...'
            #keyword list for hstack
            kwl=['matlabTime', 'julianTime', 'secondTime']
            for key in kwl:
                tmpN = getattr(newself.Variables, key)
                tmpO = getattr(StationClass.Variables, key)
                setattr(newself.Variables, key,
                np.hstack((tmpN[:], tmpO[:])))

            #keyword list for vstack
            kwl=['u', 'v', 'w', 'tke', 'gls', 'ua', 'va','el']
            kwl2D=['ua', 'va','el']
            for key in kwl:
                try:
                    if key in kwl2D:
                        tmpN = getattr(newself.Variables, key)\
                               [:,newEle[:]]
                        tmpO = getattr(StationClass.Variables, key)\
                               [:,origEle[:]]                   
                        setattr(newself.Variables, key,
                        np.vstack((tmpN[:], tmpO[:])))
                        if debug: print "Stacking " + key + "..."
                    else:
                        tmpN = getattr(newself.Variables, key)\
                               [:,:,newEle[:]]
                        tmpO = getattr(StationClass.Variables, key)\
                               [:,:,origEle[:]]                   
                        setattr(newself.Variables, key,
                        np.vstack((tmpN[:], tmpO[:])))  
                        if debug: print "Stacking " + key + "..."    
                except AttributeError:
                    continue
            #New time dimension
            newself.Grid.ntime = newself.Grid.ntime + StationClass.Grid.ntime
            #Keep only matching names
            newself.Grid.name = self.Grid.name[origEle[:]]
            #Append to new object history
            text = 'Data from ' + StationClass.History[0].split('/')[-1] \
                 + ' has been stacked'
            newself.History.append(text)

        return newself  
   
    #Methods
    def Save_as(self, filename, fileformat='pickle', debug=False):
        """
        Save the current Station structure as:
           - *.p, i.e. python file
           - *.mat, i.e. Matlab file

        Inputs:
          - filename = path + name of the file to be saved, string

        Options:
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
                print '---Data too large for machine memory---'
                print 'Tip: use ax or tx during class initialisation'
                print '---  to use partial data'
                raise
           
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

#if __name__ == '__main__':
