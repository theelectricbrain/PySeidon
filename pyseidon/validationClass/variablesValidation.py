#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import pandas as pd
import csv

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
from tidegaugeClass import TideGauge


class _load_validation:
    """
'Variables' subset in Validation class contains the following items:
-----------------------------------------------------------

                           _obs. = measurement/observational variables
    Validation.Variables._|_sim. = simulated variables
                          |_struct. = dictionnary structure for validation purposes 
    """
    def __init__(self, observed, simulated, debug=False):
        if debug: print "..variables.."
        self.obs = observed.Variables
        self.sim = simulated.Variables
        self.struct = np.array([])
        #harmonic constituents to be evaluated
        ut_constits = ['M2','S2','N2','K2','K1','O1','P1','Q1']

        #Check if times coincide
        obsMax = self.obs.matlabTime.max()
        obsMin = self.obs.matlabTime.min()
        simMax = self.sim.matlabTime.max()
        simMin = self.sim.matlabTime.min()
        absMin = max(obsMin, simMin)
        absMax = min(obsMax, simMax)
        A = set(np.where(self.sim.matlabTime[:] >= absMin)[0].tolist()) 
        B = set(np.where(self.sim.matlabTime[:] <= absMax)[0].tolist())
        C = list(A.intersection(B))
        a = set(np.where(self.obs.matlabTime[:] >= absMin)[0].tolist()) 
        b = set(np.where(self.obs.matlabTime[:] <= absMax)[0].tolist())
        c = list(a.intersection(b))        
        if len(C) == 0:
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
                sig = np.squeeze(simulated.Grid.siglay[:, ind])
     
            #Harmonic analysis over matching time
            velCoef = ut_solv(self.sim.matlabTime[C],
                              ua[C], va[C],
                              simulated.Grid.lat[ind],
                              cnstit=ut_constits, rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, conf_int=True)

            elCoef = ut_solv(self.sim.matlabTime[C],
                             el[C], [],
                             simulated.Grid.lat[ind],
                             cnstit=ut_constits, rmin=0.95, notrend=True,
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
               sig=simulated.Util3D.interpolation_at_point(simulated.Grid.siglay,
                                                           self.obs.lon, self.obs.lat)
            #Harmonic analysis
            velCoef = ut_solv(self.sim.matlabTime[C],
                              ua[C], va[C], self.obs.lat,
                              cnstit=ut_constits, rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, conf_int=True)

            elCoef = ut_solv(self.sim.matlabTime[C],
                             el[C], [], self.obs.lat,
                             cnstit=ut_constits, rmin=0.95, notrend=True,
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
                     'siglay':sig[:]}
             

        #Check what kind of observed data it is
        if observed.__module__=='pyseidon.adcpClass.adcpClass':
            obstype='ADCP'
            #Harmonic analysis
            self.obs.velCoef = ut_solv(self.obs.matlabTime[c], self.obs.ua[c],
                               self.obs.va[c], self.obs.lat,
                               cnstit=ut_constits, rmin=0.95, notrend=True,
                               method='ols', nodiagn=True, linci=True, coef_int=True)
            

            self.obs.elCoef = ut_solv(self.obs.matlabTime[c], self.obs.surf[c],
                              [], self.obs.lat,
                              cnstit=ut_constits, rmin=0.95, notrend=True,
                              method='ols', nodiagn=True, linci=True, coef_int=True)

            #Store in dict structure for compatibility purposes
            obs_mod={'ua':self.obs.ua,
                     'va':self.obs.va,
                     'elev':self.obs.surf,
                     'u':self.obs.east_vel,
                     'v':self.obs.north_vel,
                     'bins':self.obs.bins}

        #Alternative measurement type
        elif observed.__module__=='pyseidon.tidegaugeClass.tidegaugeClass':
            obstype='TideGauge'
            self.obs.elCoef = ut_solv(self.obs.matlabTime[c], self.obs.el[c],
                                      [], self.obs.lat,
                                      cnstit=ut_constits, notrend=True,
                                      rmin=0.95, method='ols', nodiagn=True,
                                      #linci=True, ordercnstit='frq')
                                      linci=True, coef_int=True)

            #Store in dict structure for compatibility purposes
            obs_mod = {'data':self.obs.RBR.data, 'elev':self.obs.el}

        else:
            print "-This type of measurements is not supported yet-"
            sys.exit()

        #Store in dict structure for compatibility purposes
        #Common block for 'struct'
        self.struct = {'name': observed.History[0].split(' ')[-1],
                       'type':obstype,
                       'lat':self.obs.lat,
                       'lon':self.obs.lon,
                       'obs_timeseries':obs_mod,
                       'mod_timeseries':sim_mod,
                       'obs_time':self.obs.matlabTime,
                       'mod_time':self.sim.matlabTime,
                       'elev_obs_harmonics':self.obs.elCoef,
                       'elev_mod_harmonics':elCoef}
        #Special blocks for 'struct'
        if self.struct['type']=='ADCP':
            self.struct['vel_obs_harmonics'] = self.obs.velCoef
            self.struct['vel_mod_harmonics'] = velCoef

        if debug: print "..done"
