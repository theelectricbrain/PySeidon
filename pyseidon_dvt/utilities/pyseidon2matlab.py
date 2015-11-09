#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
from __future__ import division
import numpy as np
import sys
from scipy.io import savemat

def pyseidon_to_matlab(fvcom, filename, exceptions=[], debug=False):
    """
    Saves fvcom object in a pickle file

    inputs:
      - fvcom = fvcom pyseidon object
      - filename = file name, string
    options:
      - exceptions = list of variables to exclude from output file
                     , list of strings
    """
    #Define bounding box
    if debug:
        print "Computing bounding box..."
    if fvcom.Grid._ax == []:
        lon = fvcom.Grid.lon[:]
        lat = fvcom.Grid.lat[:]
        fvcom.Grid._ax = [lon.min(), lon.max(),
                         lat.min(), lat.max()]

    filename = filename + ".mat"
    #TR comment: based on MitchellO'Flaherty-Sproul's code
    dtype = float
    data = {}
    Grd = {}
    Var = {}
    data['Origin'] = fvcom._origin_file
    data['History'] = fvcom.History
    Grd = fvcom.Grid.__dict__
    Var = fvcom.Variables.__dict__
    # delete exceptions
    for key in exceptions: Var.pop(key, None)
    for key in exceptions: Grd.pop(key, None)
    #TR: Force caching Variables otherwise error during loading
    #    with 'netcdf4.Variable' type (see above)
    for key in Var:
        listkeys=['Variable', 'ArrayProxy', 'BaseType'] 
        if any([type(Var[key]).__name__==x for x in listkeys]):
            if debug:
                print "Force caching for " + key
            Var[key] = Var[key][:]
        #keyV = key + '-var'
        try:
            data[key] = np.float64(Var[key].copy())
        except AttributeError:
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

        try:
            data[key] = np.float64(Grd[key].copy())
        except AttributeError:
            data[key] = Grd[key]

    #Save in mat file file
    if debug:
        print 'Dumping in matlab file...'
    savemat(filename, data, oned_as='column')