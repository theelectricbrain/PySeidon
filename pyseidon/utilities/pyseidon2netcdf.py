#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
from __future__ import division
import numpy as np
import sys
#TR comment: 2 alternatives
import netCDF4 as nc
#from scipy.io import netcdf

def pyseidon_to_netcdf(fvcom, filename, exceptions=[], compression=False, debug=False):
    """
    Saves fvcom object in a pickle file

    inputs:
      - fvcom = fvcom pyseidon object
      - filename = file name, string
    options:
      - exceptions = list of variables to exclude from output file
                     , list of strings
      - compresion = compresses data with zlib and uses at least 3 significant digits, boolean
        Note: Works only with netcdf format
    """
    #Define bounding box
    if debug: print "Computing bounding box..."
    if compression:
        zlib = True
        least_significant_digit = 3
    else:
        zlib = False
        least_significant_digit = None

    filename = filename + ".nc"
    f = nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC')
    #history attribut
    f.history = fvcom.History[:]

    #create dimensions
    if not fvcom.Variables._3D:
        ##2D dimensions
        dims = {'three':3, 'four':4,
                'nele': fvcom.Grid.nele, 'node': fvcom.Grid.nnode,
                'siglay': 2, 'siglev': 3,
                'time': fvcom.Variables.julianTime.shape[0]}
    else:
        ##3D dimensions
        dims = {'three':3, 'four':4,
                'nele': fvcom.Grid.nele, 'node': fvcom.Grid.nnode,
                'siglay': fvcom.Grid.nlevel, 'siglev': fvcom.Grid.nlevel+1,
                'vertshear':fvcom.Grid.nlevel-1,
                'time': fvcom.Variables.julianTime.shape[0]}
    for key in dims.keys():
        f.createDimension(key, dims[key])

    # list of var
    varname = fvcom.Variables.__dict__.keys()
    gridname = fvcom.Grid.__dict__.keys()
    #   getting rid of the "_var" kind
    iterlist = varname[:]
    for key in iterlist:
        if key[0] == "_": varname.remove(key)
    iterlist = gridname[:]  # getting rid of the "_var" kind
    for key in iterlist:
        if key[0] == "_": gridname.remove(key)
    #   exceptions which need name replacement
    if 'w' in varname: varname[varname.index('w')] = 'ww'

    #load in netcdf file
    if debug: print "Loading variables' matrices in nc file..."
    for var in varname:
        if not var in exceptions: # check if in list of exceptions
            try:
                if hasattr(fvcom.Variables, var):
                    dim = []
                    s = getattr(fvcom.Variables, var).shape
                    for d in s:
                        for key in dims.keys():
                            if dims[key] == d:
                                if len(dim) == len(s):  # in case of similar dims
                                    pass
                                else:
                                    dim.append(key)
                    dim = tuple(dim)
                    tmp_var = f.createVariable(var, 'float', dim,
                                               zlib=zlib, least_significant_digit=least_significant_digit)
                    tmp_var[:] = getattr(fvcom.Variables, var)[:]
            except (AttributeError, IndexError) as e:
                pass
            if debug: print "..."+var+" loaded..."

    if debug: print "Loading grid' matrices in nc file..."
    for grd in gridname:
        if not grd in exceptions: # check if in list of exceptions
            try:
                if hasattr(fvcom.Grid, grd):
                    dim = []
                    s = getattr(fvcom.Grid, grd).shape
                    for d in s:
                        for key in dims.keys():
                            if dims[key] == d:
                                if len(dim) == len(s):  # in case of similar dims
                                    pass
                                else:
                                    dim.append(key)
                    dim = tuple(dim)
                    tmp_var = f.createVariable(grd, 'float', dim)
                    tmp_var[:] = getattr(fvcom.Grid, grd)[:]
            except (AttributeError, IndexError) as e:
                pass
            if debug: print "..."+grd+" loaded..."
    f.close()
    if debug: print "...done"    






