#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import pandas as pd
import csv
import cPickle as pkl

#import netCDF4 as nc
#Quick fix
import scipy.io.netcdf as nc
from scipy.io import savemat
from scipy.io import loadmat

from datetime import datetime, timedelta
import cPickle as pickle
import sys
import os
from utide import ut_solv
import scipy.io as sio

#Local import
from compareData import *
from valTable import valTable
from variablesValidation import _load_validation
from interpolation_utils import *
from stationClass import Station
from adcpClass import ADCP
from fvcomClass import FVCOM
from tidegaugeClass import TideGauge

class Validation:
    """
    Validation class/structure.
    Functionality structured as follows:
                 _History = Quality Control metadata
    Validation._|_Variables. = observed and simulated variables and quantities
                |_validate. = validation method/function

    Inputs:
    ------
      - observed = any PySeidon measurement object (i.e. ADCP, TideGauge, Drifter,...)
      - simulated = any PySeidon simulation object (i.e. FVCOM or Station)
    """
    def __init__(self, observed, simulated, debug=False, debug_plot=False):
        self._debug = debug
        self._debug_plot = debug_plot
        if debug: print '-Debug mode on-'
        if debug: print 'Loading...'
        #Metadata
        self.History = ['Created from ' + observed._origin_file +\
                        ' and ' + simulated._origin_file]
        self.Variables = _load_validation(observed, simulated, debug=self._debug)
 
    def validate(self, filename=[], depth=[], plot=False, debug=False, debug_plot=False):
        """
        This method computes series of standard validation benchmarks.

        Options:
        ------
          - filename: file name of the .csv file to be saved, string.
          - depth: depth at which the validation will be performed, float.
                   Only applicable for 3D simulations.
          - plot: plot series of valiudation graphs, boolean.

        References:
        ----------
        - NOAA. NOS standards for evaluating operational nowcast and
          forecast hydrodynamic model systems, 2003.

        - K. Gunn, C. Stock-Williams. On validating numerical hydrodynamic
          models of complex tidal flow, International Journal of Marine Energy, 2013

        - N. Georgas, A. Blumberg. Establishing Confidence in Marine Forecast
          Systems: The design and skill assessment of the New York Harbor Observation
          and Prediction System, version 3 (NYHOPS v3), 2009

        - Liu, Y., P. MacCready, B. M. Hickey, E. P. Dever, P. M. Kosro, and
          N. S. Banas (2009), Evaluation of a coastal ocean circulation model for
          the Columbia River plume in summer 2004, J. Geophys. Res., 114
        """
        debug = debug or self._debug
        debug_plot = debug_plot or self._debug_plot    
        #User input
        if filename==[]:
            filename = input('Enter filename (string) for csv file: ')
            filename = str(filename)
        if (depth==[] and self.Variables.sim._3D):
            depth = input('Depth from surface at which the validation will be performed: ')
            depth = float(depth)
            if depth < 0.0: depth = -1.0 * depth
        if depth==[]: depth=5.0

        #initialisation
        vars = []

        if self.Variables.struct['type'] == 'ADCP':
    	    (elev_suite, speed_suite, dir_suite, u_suite, v_suite, 
             vel_suite) = compareUV(self.Variables.struct, self.Variables.sim._3D,
                                    plot=plot, depth=depth,
                                    debug=debug, debug_plot=debug_plot)
            self.Variables.struct['elev_val'] = elev_suite
    	    self.Variables.struct['speed_val'] = speed_suite
    	    self.Variables.struct['dir_val'] = dir_suite
            self.Variables.struct['u_val'] = u_suite
	    self.Variables.struct['v_val'] = v_suite
	    self.Variables.struct['vel_val'] = vel_suite
            #Variable to processed
            vars.append('elev')
            vars.append('speed')
            vars.append('dir')
            vars.append('u')
            vars.append('v')
            vars.append('vel')

        elif self.Variables.struct['type'] == 'TideGauge':
     	    elev_suite_dg = compareTG(self.Variables.struct,
                                      debug=debug, debug_plot=debug_plot)
    	    self.Variables.struct['tg_val'] = elev_suite_dg 
            #Variable to processed
            vars.append('tg')

        else:
            print "-This type of measurements is not supported yet-"
            sys.exit()

        #Make csv file
        self.Benchmarks = valTable(self.Variables.struct, filename,  vars,
                                   debug=debug, debug_plot=debug_plot)
        #Display csv
        csvName = filename + '_val.csv'
        csv_con = open(csvName, 'r')
        csv_cont = list(csv.reader(csv_con, delimiter=','))
        print "---Validation benchmarks---"
        print(70*'-')
        for row in csv_cont:
           row = [str(e) for e in row[:][1:]]
           print('\t'.join(row))
        print(70*'-')  

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
        if debug: print 'Saving file...'

        #Save as different formats
        if fileformat=='pickle':
            filename = filename + ".p"
            f = open(filename, "wb")
            data = {}
            data['History'] = self.History
            data['Benchmarks'] = self.Benchmarks
            data['Variables'] = self.Variables.__dict__
            #TR: Force caching Variables otherwise error during loading
            #    with 'netcdf4.Variable' type (see above)
            for key in data['Variables']:
                listkeys=['Variable', 'ArrayProxy', 'BaseType'] 
                if any([type(data['Variables'][key]).__name__==x for x in listkeys]):
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
                raise
           
            f.close()
        elif fileformat=='matlab':
            filename = filename + ".mat"
            #TR comment: based on MitchellO'Flaherty-Sproul's code
            dtype = float
            data = {}
            Grd = {}
            Var = {}
            Bch = {}

            data['History'] = self.History
            Bch = self.Benchmarks
            for key in Bch:
                data[key] = Bch[key]   
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
          
            #Save in mat file file
            if debug:
                print 'Dumping in matlab file...'
            savemat(filename, data, oned_as='column')       
        else:
            print "---Wrong file format---"              
