#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import pandas as pd

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
from interpolation_utils import *
from stationClass import Station
from adcpClass import ADCP
from fvcomClass import FVCOM
from tidegaugeClass import Tidegauge

class Validation:
    """ """
    def __init__(self, observed, simulated, debug=False):
        self._debug = debug
        self.obs = observed.Variables
        self.sim = simulated.Variables
        self.struct = np.array([])

        #Check what kind of simulated data it is
        if simulated.__module__=='pyseidon.stationClass.stationClass':
            #Find closest point to ADCP
            ind = closest_point([self.obs.lon], [self.obs.lat],
                                simulated.Grid.lon[:],
                                simulated.Grid.lat[:])
            nameSite = ''.join(simulated.Grid.name[ind,:][0,:])
            print "Station site: " + nameSite
            self.sim.el = self.sim.el[:, ind].flatten()
            self.sim.ua = self.sim.ua[:, ind].flatten()
            self.sim.va=self.sim.va[:, ind].flatten()           
            #Harmonic analysis
            self.sim.velCoef = ut_solv(self.sim.matlabTime[:],
                                       self.sim.ua[:],
                                       self.sim.va[:],
                                       simulated.Grid.lat[ind],
                               cnstit='auto', rmin=0.95, notrend=True,
                               method='ols', nodiagn=True, linci=True, conf_int=True)

            self.sim.elCoef = ut_solv(self.sim.matlabTime[:],
                                      self.sim.el, [],
                                      simulated.Grid.lat[ind],
                              cnstit='auto', rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, conf_int=True)
        #Alternative simulation type
        elif simulated.__module__=='pyseidon.fvcomClass.fvcomClass':
            #Interpolation at measurement location
            self.sim.el=simulated.Util2D.interpolation_at_point(self.sim.el,
                                                       self.obs.lon, self.obs.lat)
            self.sim.ua=simulated.Util2D.interpolation_at_point(self.sim.ua,
                                                       self.obs.lon, self.obs.lat)
            self.sim.va=simulated.Util2D.interpolation_at_point(self.sim.va,
                                                       self.obs.lon, self.obs.lat)
            if self.sim._3D:
               self.sim.u=simulated.Util3D.interpolation_at_point(self.sim.u,
                                                           self.obs.lon, self.obs.lat)
               self.sim.v=simulated.Util3D.interpolation_at_point(self.sim.v,
                                                           self.obs.lon, self.obs.lat)
            #Harmonic analysis
            self.sim.velCoef = ut_solv(self.sim.matlabTime[:],
                                       ua[:],
                                       va[:],
                                       self.obs.lat,
                               cnstit='auto', rmin=0.95, notrend=True,
                               method='ols', nodiagn=True, linci=True, conf_int=True)

            self.sim.elCoef = ut_solv(self.sim.matlabTime[:],
                                      el[:], [],
                                      self.obs.lat,
                              cnstit='auto', rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, conf_int=True)

        else:
            print "-This type of simulations is not supported yet-"
            sys.exit()

        #Store in dict structure for compatibility purposes
        if not self.sim._3D:
            sim_mod={'ua':self.sim.ua[:],
                     'va':self.sim.va[:],
                     'elev':self.sim.el[:]}
        else:
            sim_mod={'ua':self.sim.ua[:],
                     'va':self.sim.va[:],
                     'elev':self.sim.el[:],
                     'u':self.sim.u[:],
                     'v':self.sim.v[:]}
             

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
        #    obstype='tidegauge'
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
                        'vel_mod_harmonics':self.sim.velCoef,
                        'elev_mod_harmonics':self.sim.elCoef}    
    #def compare
