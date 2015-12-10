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

    # Check if netcdffilename has an extension
    if not filename[-3:] == '.nc':
        filename = filename + '.nc'
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
    mstrList = varname + gridname

    #load in netcdf file
    if debug: print "Loading in nc file..."
    for var in mstrList:
        if not var in exceptions: # check if in list of exceptions
            if var in varname:
                data = fvcom.Variables
            else:
                data = fvcom.Grid
                zlib = False
                least_significant_digit = None
            try:
                if hasattr(data, var):
                    dim = []
                    s = getattr(data, var).shape
                    for d in s:
                        flag = 1
                        count = 0
                        while flag:
                            if count == len(dims.keys()):  # when two dimensions are the same by coincidence
                                #dim.append(dim[-1])  # TODO search in the entire list == d rather than last index
                                for k in dim:
                                    if dims[k] == d:
                                        newkey = k
                                dim.append(newkey)
                                flag = 0
                            else:
                                key = dims.keys()[count]
                                if dims[key] == d:
                                    if len(dim) == len(s):  # in case of similar dims
                                        count += 1
                                        pass
                                    else:
                                        if key not in dim:  # make sure dimension doesn't get recorded twice
                                            dim.append(key)
                                            flag = 0
                                            count += 1
                                        else:
                                            count += 1
                                else:
                                    count += 1
                    dim = tuple(dim)
                    #   exceptions which need name replacement
                    if var == 'w':
                        keyAlias = 'ww'
                    elif var == 'el':
                        keyAlias = 'zeta'
                    else:
                        keyAlias = var
                    tmp_var = f.createVariable(keyAlias, 'float', dim,
                                               zlib=zlib, least_significant_digit=least_significant_digit)
                    tmp_var[:] = getattr(data, var)[:]
            except (AttributeError, IndexError) as e:
                pass
            if debug: print "..."+var+" loaded..."
    # clean exit
    f.close()
