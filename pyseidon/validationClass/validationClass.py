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
from utilities import *
from stationClass import Station
from adcpClass import ADCP
from fvcomClass import FVCOM
from tidegaugeClass import Tidegauge

class Validation:
    """ """
    def __init__(self, observed, simulated):
        self.obs = observed.Variables
        self.sim = simulated.Variables
        self.struct = np.array([])
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
            #Location of the measurement, needed further down
            lonlat = np.array([self.obs.lon, self.obs.lat]).T

        #Alternative measurement type
        #elif observed.__module__=='pyseidon.tidegaugeClass.tidegaugeClass':
        #    obstype='tidegauge'

        else:
            print "-This type of measurements is not supported yet-"
            sys.exit()

        #Check what kind of simulated data it is
        if simulated.__module__=='pyseidon.stationClass.stationClass':
            #Find closest point to ADCP
            ind = closest_point(lonlat, simulated.Grid.lon[:], simulated.Grid.lat[:])
            print "Station site: " + ''.join(simulated.Grid.name[ind,:])
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
                                                       lonlat[0], lonlat[1])
            self.sim.ua=simulated.Util2D.interpolation_at_point(self.sim.ua,
                                                       lonlat[0], lonlat[1])
            self.sim.va=simulated.Util2D.interpolation_at_point(self.sim.va,
                                                       lonlat[0], lonlat[1])
            if self.sim._3D:
               self.sim.u=simulated.Util3D.interpolation_at_point(self.sim.u,
                                                           lonlat[0], lonlat[1])
               self.sim.v=simulated.Util3D.interpolation_at_point(self.sim.v,
                                                           lonlat[0], lonlat[1])
            #Harmonic analysis
            self.sim.velCoef = ut_solv(self.sim.matlabTime[:],
                                       ua[:],
                                       va[:],
                                       lonlat[1],
                               cnstit='auto', rmin=0.95, notrend=True,
                               method='ols', nodiagn=True, linci=True, conf_int=True)

            self.sim.elCoef = ut_solv(self.sim.matlabTime[:],
                                      el[:], [],
                                      lonlat[1],
                              cnstit='auto', rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, conf_int=True)

        else:
            print "-This type of simulations is not supported yet-"
            sys.exit()

       #Store in dict structure for compatibility purposes
       obs_mod={'ua':self.obs.ua,
                'va':self.obs.va,
                'elev':self.obs.surf,
                'u':self.obs.east_vel,
                'v':self.obs.north_vel,
                'bins':self.obs.bins}

        if not self.sim._3D:
            sim_obs={'ua':self.sim.ua[:],
                     'va':self.sim.va[:],
                     'elev':self.sim.el[:]}
        else:
            sim_obs={'ua':self.sim.ua[:],
                     'va':self.sim.va[:],
                     'elev':self.sim.el[:],
                     'u':self.sim.u[:],
                     'v':self.sim.v[:]}
        self.struct = {'name': observed.History[0].split(' ')[-1],
                        'type':obstype,
                        'lat':self.obs.lat[0],
                        'lon':self.obs.lon[0],
                        'obs_timeseries':obs_mod,
                        'mod_timeseries':sim_mod,
                        'obs_time':self.obs.mtime,
                        'mod_time':self.sim.matlabTime,
                        'vel_obs_harmonics':self.obs.velCoef,
                        'elev_obs_harmonics':self.obs.elCoef,
                        'vel_mod_harmonics':self.sim.velCoef,
                        'elev_mod_harmonics':self.sim.elCoef}                 


