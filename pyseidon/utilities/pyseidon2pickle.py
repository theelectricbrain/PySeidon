#!/usr/bin/python2.7
# encoding: utf-8

#Libs import
from __future__ import division
import cPickle as pkl

#Local import
from functionsFvcomThreeD import *

# Custom error
from pyseidon_error import PyseidonError

def pyseidon_to_pickle(fvcom, filename, exceptions=[], debug=False):
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
    filename = filename + ".p"
    f = open(filename, "wb")
    data = {}
    data['Origin'] = fvcom._origin_file
    data['History'] = fvcom.History
    data['Grid'] = fvcom.Grid.__dict__
    data['Variables'] = fvcom.Variables.__dict__
    # delete exceptions
    for key in exceptions: data['Variables'].pop(key, None)
    for key in exceptions: data['Grid'].pop(key, None)
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
            raise PyseidonError("---Data too large for machine memory---\n"\
                                "Tip: use ax or tx during class initialisation\n"\
                                "     to use partial data")
   
    f.close()
