#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import pandas as pd
import csv

#import netCDF4 as nc
#Quick fix
import scipy.io.netcdf as nc

from datetime import datetime, timedelta
import cPickle as pickle
import sys
import os
from utide import ut_solv
import scipy.io as sio

#Local import
from compareData import *
from valTable import valTable
from interpolation_utils import *
from stationClass import Station
from adcpClass import ADCP
from fvcomClass import FVCOM
from tidegaugeClass import Tidegauge

class Validation:
    """
    Validation class/structure.
    Functionality structured as follows:
                 _obs. = measurement/observational variables
    Validation._|_sim. = simulated variables
                |_struct = dictionnary structure for validation purposes
                |_validate. = validation method/function

    Inputs:
    ------
      - observed = any PySeidon measurement object (i.e. ADCP, TideGauge, Drifter,...)
      - simulated = any PySeidon simulation object (i.e. FVCOM or Station)
    """
    def __init__(self, observed, simulated, debug=False):
        self._debug = debug
        self.obs = observed.Variables
        self.sim = simulated.Variables
        self.struct = np.array([])

        #Check if times coincide
        obsMax = self.obs.mtime.max()
        obsMin = self.obs.mtime.min()
        simMax = self.sim.matlabTime.max()
        simMin = self.sim.matlabTime.min()
        absMin = max(obsMin, simMin)
        absMax = min(obsMax, simMax)
        A = set(np.where(self.sim.matlabTime[:] >= absMin)[0].tolist()) 
        B = set(np.where(self.sim.matlabTime[:] <= absMax)[0].tolist())
        test = A.intersection(B) 
        if len(test) == 0:
           print "---Time between simulation and measurement does not match up---"
           sys.exit()

        #Check what kind of simulated data it is
        if simulated.__module__=='pyseidon.stationClass.stationClass':
            #Find closest point to ADCP
            ind = closest_point([self.obs.lon], [self.obs.lat],
                                simulated.Grid.lon[:],
                                simulated.Grid.lat[:])
            nameSite = ''.join(simulated.Grid.name[ind,:][0,:])
            print "Station site: " + nameSite
            el = self.sim.el[:, ind].flatten()
            ua = self.sim.ua[:, ind].flatten()
            va = self.sim.va[:, ind].flatten()
            if self.sim._3D:
                u = np.squeeze(self.sim.u[:, :,ind])
                v = np.squeeze(self.sim.v[:, :,ind])
     
            #Harmonic analysis
            velCoef = ut_solv(self.sim.matlabTime[:],
                              ua[:], va[:],
                              simulated.Grid.lat[ind],
                              cnstit='auto', rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, conf_int=True)

            elCoef = ut_solv(self.sim.matlabTime[:],
                             el, [],
                             simulated.Grid.lat[ind],
                             cnstit='auto', rmin=0.95, notrend=True,
                             method='ols', nodiagn=True, linci=True, conf_int=True)
        #Alternative simulation type
        elif simulated.__module__=='pyseidon.fvcomClass.fvcomClass':
            #Interpolation at measurement location
            el=simulated.Util2D.interpolation_at_point(self.sim.el,
                                                       self.obs.lon, self.obs.lat)
            ua=simulated.Util2D.interpolation_at_point(self.sim.ua,
                                                       self.obs.lon, self.obs.lat)
            va=simulated.Util2D.interpolation_at_point(self.sim.va,
                                                       self.obs.lon, self.obs.lat)
            if self.sim._3D:
               u=simulated.Util3D.interpolation_at_point(self.sim.u,
                                                         self.obs.lon, self.obs.lat)
               v=simulated.Util3D.interpolation_at_point(self.sim.v,
                                                         self.obs.lon, self.obs.lat)
            #Harmonic analysis
            velCoef = ut_solv(self.sim.matlabTime[:],
                              ua[:], va[:], self.obs.lat,
                              cnstit='auto', rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, conf_int=True)

            elCoef = ut_solv(self.sim.matlabTime[:],
                             el[:], [], self.obs.lat,
                             cnstit='auto', rmin=0.95, notrend=True,
                             method='ols', nodiagn=True, linci=True, conf_int=True)

        else:
            print "-This type of simulations is not supported yet-"
            sys.exit()

        #Store in dict structure for compatibility purposes
        if not self.sim._3D:
            sim_mod={'ua':ua[:],
                     'va':va[:],
                     'elev':el[:]}
        else:
            sim_mod={'ua':ua[:],
                     'va':va[:],
                     'elev':el[:],
                     'u':u[:],
                     'v':v[:],
                     'siglay':np.squeeze(simulated.Grid.siglay[:, ind])}
             

        #Check what kind of observed data it is
        if observed.__module__=='pyseidon.adcpClass.adcpClass':
            obstype='ADCP'
            #Harmonic analysis
            self.obs.velCoef = ut_solv(self.obs.mtime, self.obs.ua,
                               self.obs.va, self.obs.lat,
                               cnstit='auto', rmin=0.95, notrend=True,
                               method='ols', nodiagn=True, linci=True, coef_int=True)
            

            self.obs.elCoef = ut_solv(self.obs.mtime, self.obs.surf,
                              [], self.obs.lat,
                              cnstit='auto', rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, coef_int=True)
            #Store in dict structure for compatibility purposes
            obs_mod={'ua':self.obs.ua,
                     'va':self.obs.va,
                     'elev':self.obs.surf,
                     'u':self.obs.east_vel,
                     'v':self.obs.north_vel,
                     'bins':self.obs.bins}

        #Alternative measurement type
        #elif observed.__module__=='pyseidon.tidegaugeClass.tidegaugeClass':
        #    obstype='TideGauge'
        #    obs_mod = {'data':tideData.data, 'elev':tideData.elev}

        else:
            print "-This type of measurements is not supported yet-"
            sys.exit()

        self.struct = {'name': observed.History[0].split(' ')[-1],
                        'type':obstype,
                        'lat':self.obs.lat,
                        'lon':self.obs.lon,
                        'obs_timeseries':obs_mod,
                        'mod_timeseries':sim_mod,
                        'obs_time':self.obs.mtime,
                        'mod_time':self.sim.matlabTime,
                        'vel_obs_harmonics':self.obs.velCoef,
                        'elev_obs_harmonics':self.obs.elCoef,
                        'vel_mod_harmonics':velCoef,
                        'elev_mod_harmonics':elCoef}
 
    def validate(self, filename=[], depth=[], plot=False):
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
        #User input
        if filename==[]:
            filename = input('Enter filename (string) for csv file: ')
            filename = str(filename)
        if (depth==[] and self.sim._3D):
            depth = input('Depth from surface at which the validation will be performed: ')
            depth = float(depth)
            if depth < 0.0: depth = -1.0 * depth
        if depth==[]: depth=5.0

        #initialisation
        vars = []

        if self.struct['type'] == 'ADCP':
    	    (elev_suite, speed_suite, dir_suite, u_suite, v_suite, 
             vel_suite) = compareUV(self.struct, self.sim._3D, plot=plot, depth=depth)
            self.struct['elev_val'] = elev_suite
    	    self.struct['speed_val'] = speed_suite
    	    self.struct['dir_val'] = dir_suite
            self.struct['u_val'] = u_suite
	    self.struct['v_val'] = v_suite
	    self.struct['vel_val'] = vel_suite
            #Variable to processed
            vars.append('elev')
            vars.append('speed')
            vars.append('dir')

        elif self.struct['type'] == 'TideGauge':
     	    elev_suite_dg = compareTG(self.struct)
    	    self.struct['tg_val'] = elev_suite_dg 
            #Variable to processed
            vars.append('tg')

        else:
            print "-This type of measurements is not supported yet-"
            sys.exit()

        #Make csv file
        valTable(self.struct, filename,  vars)
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
