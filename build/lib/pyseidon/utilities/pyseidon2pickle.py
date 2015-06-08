#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
from __future__ import division
import numpy as np
import sys
#TR comment: 2 alternatives
#import netCDF4 as nc
from scipy.io import netcdf
import cPickle as pkl

#Local import
from variablesFvcom import _load_var, _load_grid
from functionsFvcom import *
from functionsFvcomThreeD import *
from plotsFvcom import *

def pyseidon_to_pickle(fvcom, filename, debug):
    """
    save fvcom object in a pickle file

    inputs:
      - fvcom = fvcom pyseidon object
      - filename = file name, string
    """
    #Define bounding box
    if debug:
        print "Computing bounding box..."
    if fvcom.Grid._ax == []:
        lon = fvcom.Grid.lon[:]
        lat = fvcom.Grid.lat[:]
        fvcom.Grid._ax = [lon.min(), lon.max(),
                         lat.min(), lat.max()]
    filename = filename + ".p"
    f = open(filename, "wb")
    data = {}
    data['Origin'] = fvcom._origin_file
    data['History'] = fvcom.History
    data['Grid'] = fvcom.Grid.__dict__
    data['Variables'] = fvcom.Variables.__dict__
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
