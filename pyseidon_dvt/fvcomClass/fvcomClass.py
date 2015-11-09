#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
from __future__ import division
#TR comment: 2 alternatives
import netCDF4 as nc
from scipy.io import netcdf
from pydap.client import open_url
import cPickle as pkl
import pickle as Pkl
import copy
from os.path import isfile
import gc

#Utility import
from pyseidon_dvt.utilities.object_from_dict import ObjectFromDict
from pyseidon_dvt.utilities.pyseidon2pickle import pyseidon_to_pickle
from pyseidon_dvt.utilities.pyseidon2matlab import pyseidon_to_matlab
from pyseidon_dvt.utilities.pyseidon2netcdf import pyseidon_to_netcdf

# Custom error
from pyseidon_error import PyseidonError

#Local import
from variablesFvcom import _load_var, _load_grid
from functionsFvcom import *
from functionsFvcomThreeD import *
from plotsFvcom import *

class FVCOM:
    """
    **A class/structure for FVCOM data**
    Functionality structured as follows: ::

                _Data. = raw netcdf file data
               |_Variables. = fvcom variables and quantities
               |_Grid. = fvcom grid data
               |_History = Quality Control metadata
        FVCOM._|_Utils2D. = set of useful functions and methods for 2D and 3D runs
               |_Utils3D. = set of useful functions and methods for 3D runs
               |_Plots. = plotting functions
               |_Save_as = "save as" methods

    Inputs:
      - filename = path to file, string,
                ex: testFvcom = FVCOM('./path_to_FVOM_output_file/filename')
                Note that the file can be a pickle file (i.e. *.p) or a netcdf file (i.e. *.nc)
                Additionally, either a file path or a OpenDap url could be used

    Options:
      - ax = defines for a specific spatial region to work with, as such:
           ax = [minimun longitude, maximun longitude, minimun latitude, maximum latitude]
           or use one of the following pre-defined region: ax = 'GP', 'PP', 'DG' or 'MP'
           Note that this option permits to extract partial data from the overall file
           and therefore reduce memory and cpu use.

      - tx = defines for a specific temporal period to work with, as such:
           tx = ['2012-11-07T12:00:00','2012.11.09 12:00:00'], string of 'yyyy-mm-dd hh:mm:ss'
           Note that this option permits to extract partial data from the overall file
           and therefore reduce memory and cpu use

    *Notes*
    Throughout the package, the following conventions apply:
      - Date = string of 'yyyy-mm-dd hh:mm:ss'
      - Coordinates = decimal degrees East and North
      - Directions = in degrees, between -180 and 180 deg., i.e. 0=East, 90=North, +/-180=West, -90=South
      - Depth = 0m is the free surface and depth is negative
    """

    def __init__(self, filename, ax=[], tx=[], debug=False):
        """ Initialize FVCOM class."""
        self._debug = debug
        if debug: print '-Debug mode on-'
        #Force garbage collector when fvcom object created
        gc.collect()

        #Loading pickle file
        if filename.endswith('.p'):
            f = open(filename, "rb")
            try:
                data = pkl.load(f)
            except MemoryError:
                try:
                    data = Pkl.load(f)
                except KeyError:
                    data = pkl.load(f,2)
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
                    if isfile(data['Origin']):
                        try:
                            self.Data = netcdf.netcdf_file(data['Origin'], 'r',mmap=True)
                            #due to mmap not coping with big array > 4Gib
                        except (OverflowError, TypeError, ValueError) as e:
                            self.Data = nc.Dataset(data['Origin'], 'r',
                                       format='NETCDF4_CLASSIC')
                    else:
                        print "the original *.nc file has not been found"
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
                try:
                    self.Data = netcdf.netcdf_file(filename, 'r',mmap=True)
                    #due to mmap not coping with big array > 4Gib
                except (OverflowError, TypeError, ValueError) as e:
                    self.Data = nc.Dataset(filename, 'r', format='NETCDF4_CLASSIC')
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
            raise PyseidonError("---Functionality not yet implemented---")
        else:
            raise PyseidonError("---Wrong file format---")

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

        ##Re-assignement of utility functions as methods
        #self.dump_profile_data = self.Plots._dump_profile_data_as_csv
        #self.dump_map_data = self.Plots._dump_map_data_as_csv

        return

    #Special methods
    def __del__(self):
        """making sure that all opened files are closed when deleted or overwritten"""
        #TR: not sure __del__ is the best approach for that
        try:
            if type(self.Data).__name__ == "netcdf_file":
                try:
                    self.Data.close()
                except AttributeError:
                    pass
            elif type(self.Data).__name__ == "Dataset":
                self.Data.close()
            else:
                try:
                    f.close()
                except (NameError,AttributeError) as e:
                    pass
        except AttributeError:
            try:
                f.close()
            except (NameError,AttributeError) as e:
                pass

    def __new__(self):
        """Force garbage collector when new fvcom object created"""
        gc.collect()

    def __add__(self, FvcomClass, debug=False):
        """
        This special method permits to stack variables
        of 2 FVCOM objects through a simple addition: ::
          fvcom1 += fvcom2

        *Notes*
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
            raise PyseidonError("---Spatial regions do not match---")
        elif not ((self.Grid.nele == FvcomClass.Grid.nele) and
                  (self.Grid.nnode == FvcomClass.Grid.nnode) and
                  (self.Variables._3D == FvcomClass.Variables._3D)):
            raise PyseidonError("---Data dimensions do not match---")
        else:
            if not (self.Variables.julianTime[-1]<=
                    FvcomClass.Variables.julianTime[0]):
                raise PyseidonError("---Data not consecutive in time---")
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
    def save_as(self, filename, fileformat='netcdf', exceptions=[], compression=False, debug=False):
        """
        This method saves the current FVCOM structure as:
           - *.nc, i.e. netcdf file
           - *.p, i.e. python file
           - *.mat, i.e. Matlab file

        Inputs:
          - filename = path + name of the file to be saved, string

        Options:
          - fileformat = format of the file to be saved, i.e. 'pickle', .netcdf. or 'matlab'
          - exceptions = list of variables to exclude from output file
                     , list of strings
          - compresion = compresses data with zlib and uses at least 3 significant digits, boolean
            Note: Works only with netcdf format
        """
        debug = debug or self._debug
        if debug:
            print 'Saving file...'
        #Save as different formats
        if fileformat=='pickle':
            pyseidon_to_pickle(self, filename, exceptions=exceptions, debug=debug)
        elif fileformat=='matlab':
            pyseidon_to_matlab(self, filename, exceptions=exceptions, debug=debug)
        elif fileformat=='netcdf':
            pyseidon_to_netcdf(self, filename, exceptions=exceptions, compression=compression, debug=debug)
        else:
            print "---Wrong file format---"

#Test section when running in shell >> python fvcomClass.py
#if __name__ == '__main__':

    #filename = './test_file/dn_coarse_0001.nc'
    #test = FVCOM(filename)
    #test.harmonics(0, cnstit='auto', notrend=True, nodiagn=True)
    #WB_COMMENTS: fixed matlabttime to matlabtime
    #test.reconstr(test.Variables.matlabTime)
