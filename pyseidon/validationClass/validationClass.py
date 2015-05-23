#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import pandas as pd
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
                |_Variables. = observed and simulated variables and quantities
                |_validate_data = validation method/function against timeseries
    Validation._|_validate_harmonics = validation method/function against 
                |                      harmonic coefficients
                |_Save_as = "save as" function 

    Inputs:
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
 
    def validate_data(self, filename=[], depth=[], plot=False, save_csv=False, 
                      debug=False, debug_plot=False):
        """
        This method computes series of standard validation benchmarks.

        Options:
          - filename: file name of the .csv file to be saved, string.
          - depth: depth at which the validation will be performed, float.
                   Only applicable for 3D simulations.
          - plot: plot series of valiudation graphs, boolean.
          - save_csv: will save both observed and modeled interpolated
                      timeseries into *.csv file  

        References:
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
             vel_suite, csp_suite) = compareUV(self.Variables.struct, self.Variables.sim._3D,
                                    plot=plot, depth=depth, #save_csv=save_csv,
                                    debug=debug, debug_plot=debug_plot)
            self.Variables.struct['elev_val'] = elev_suite
            self.Variables.struct['speed_val'] = speed_suite
            self.Variables.struct['dir_val'] = dir_suite
            self.Variables.struct['u_val'] = u_suite
            self.Variables.struct['v_val'] = v_suite
            self.Variables.struct['vel_val'] = vel_suite
            # custom benchmark
            self.Variables.struct['cubic_speed_val'] = csp_suite
            # Variable to processed
            vars.append('elev')
            vars.append('speed')
            vars.append('dir')
            vars.append('u')
            vars.append('v')
            vars.append('vel')
            # custom var
            vars.append('cubic_speed')

        elif self.Variables.struct['type'] == 'TideGauge':
            elev_suite_dg = compareTG(self.Variables.struct,
                                      plot=plot, #save_csv=save_csv,
                                      debug=debug, debug_plot=debug_plot)
            self.Variables.struct['tg_val'] = elev_suite_dg
            #Variable to processed
            vars.append('tg')

        elif self.Variables.struct['type'] == 'Drifter':
            (elev_suite, speed_suite, dir_suite, u_suite, v_suite,
             vel_suite, csp_suite) = compareUV(self.Variables.struct, self.Variables.sim._3D,
                                    plot=plot, depth=depth, #save_csv=save_csv,
                                    debug=debug, debug_plot=debug_plot)

            self.Variables.struct['speed_val'] = speed_suite
            self.Variables.struct['dir_val'] = dir_suite
            self.Variables.struct['u_val'] = u_suite
            self.Variables.struct['v_val'] = v_suite
            # custom benchmark
            self.Variables.struct['vel_val'] = vel_suite
            self.Variables.struct['cubic_speed_val'] = csp_suite

            # Variable to processed
            vars.append('speed')
            vars.append('dir')
            vars.append('u')
            vars.append('v')
            vars.append('vel')
            # custom var
            vars.append('vel')
            vars.append('cubic_speed')

        else:
            print "-This kind of measurements is not supported yet-"
            sys.exit()

        #Make csv file
        self.Benchmarks = valTable(self.Variables.struct, filename,  vars,
                                   save_csv=save_csv, debug=debug, debug_plot=debug_plot)
        
        #Display csv
        #csvName = filename + '_val.csv'
        #csv_con = open(csvName, 'r')
        #csv_cont = list(csv.reader(csv_con, delimiter=','))
        print "---Validation benchmarks---"
        pd.set_option('display.max_rows', len(self.Benchmarks))
        print(self.Benchmarks)
        pd.reset_option('display.max_rows')
        #print(70*'-')
        #for row in csv_cont:
        #   row = [str(e) for e in row[:][1:]]
        #   print('\t'.join(row))
        #print(70*'-')

    def validate_harmonics(self, filename=[], save_csv=False,
                           debug=False, debug_plot=False):
        """
        This method computes and store in a csv file the error in %
        for each component of the harmonic analysis (i.e. *_error.csv).     

        Options:
          - filename: file name of the .csv file to be saved, string.
          - save_csv: will save both observed and modeled harmonic
                      coefficients into *.csv files (i.e. *_harmo_coef.csv) 
        """
        #User input
        if filename==[]:
            filename = input('Enter filename (string) for csv file: ')
            filename = str(filename)


        #Harmonic analysis over matching time
        if self.Variables._obstype=='adcp':
            time = self.Variables.struct['obs_time']
            lat = self.Variables.struct['lat']
            ua =  self.Variables.struct['obs_timeseries']['ua'][:]
            va =  self.Variables.struct['obs_timeseries']['va'][:]
            el =  self.Variables.struct['obs_timeseries']['elev'] [:]          
            
            self.Variables.obs.velCoef = ut_solv(time, ua, va, lat,
                                         #cnstit=ut_constits, rmin=0.95, notrend=True,
                                         cnstit='auto', rmin=0.95, notrend=True,
                                         method='ols', nodiagn=True, linci=True,
                                         coef_int=True)
            

            self.Variables.obs.elCoef = ut_solv(time, el, [], lat,
                                        #cnstit=ut_constits, rmin=0.95, notrend=True,
                                        cnstit='auto', rmin=0.95, notrend=True,
                                        method='ols', nodiagn=True, linci=True,
                                        coef_int=True)

        elif self.Variables._obstype=='tidegauge':
            time = self.Variables.struct['obs_time']
            lat = self.Variables.struct['lat']
            el =  self.Variables.struct['obs_timeseries']['elev'] [:]
 
            self.Variables.obs.elCoef = ut_solv(time, el, [], lat,
                                        #cnstit=ut_constits, notrend=True,
                                        cnstit='auto', notrend=True,
                                        rmin=0.95, method='ols', nodiagn=True,
                                        #linci=True, ordercnstit='frq')
                                        linci=True, coef_int=True)
        else:
            print "--This kind of observations is not supported---"
            sys.exit()

        if self.Variables._simtype=='fvcom':
            time = self.Variables.struct['mod_time']
            lat = self.Variables.struct['lat']
            el =  self.Variables.struct['mod_timeseries']['elev'][:]           
            
            self.Variables.sim.elCoef = ut_solv(time, el, [], lat,
                             #cnstit=ut_constits, rmin=0.95, notrend=True,
                             cnstit='auto', rmin=0.95, notrend=True,
                             method='ols', nodiagn=True, linci=True, conf_int=True)
            if self.Variables._obstype=='adcp':
                ua =  self.Variables.struct['mod_timeseries']['ua'][:]
                va =  self.Variables.struct['mod_timeseries']['va'][:]
                self.Variables.sim.velCoef = ut_solv(time, ua, va, lat,
                                  #cnstit=ut_constits, rmin=0.95, notrend=True,
                                  cnstit='auto', rmin=0.95, notrend=True,
                                  method='ols', nodiagn=True, linci=True, conf_int=True)

        elif self.Variables._simtype=='station':
            time = self.Variables.struct['mod_time']
            lat = self.Variables.struct['lat']
            el = self.Variables.struct['mod_timeseries']['elev'][:]

            self.Variables.sim.elCoef = ut_solv(time, el, [], lat,
                             #cnstit=ut_constits, rmin=0.95, notrend=True,
                             cnstit='auto', rmin=0.95, notrend=True,
                             method='ols', nodiagn=True, linci=True, conf_int=True)
            if self.Variables._obstype=='adcp':
                ua = self.Variables.struct['mod_timeseries']['ua'][:]
                va = self.Variables.struct['mod_timeseries']['va'][:]
                self.Variables.sim.velCoef = ut_solv(time, ua, va, lat,
                                  #cnstit=ut_constits, rmin=0.95, notrend=True,
                                  cnstit='auto', rmin=0.95, notrend=True,
                                  method='ols', nodiagn=True, linci=True, conf_int=True)

        #find matching and non-matching coef
        matchElCoef = []
        matchElCoefInd = []
        for i1, key1 in enumerate(self.Variables.sim.elCoef['name']):
            for i2, key2 in enumerate(self.Variables.obs.elCoef['name']):
                if key1 == key2:
                   matchElCoefInd.append((i1,i2))
                   matchElCoef.append(key1)
        matchElCoefInd=np.array(matchElCoefInd)
        noMatchElCoef = np.delete(self.Variables.sim.elCoef['name'],
                                  matchElCoefInd[:,0])
        np.hstack((noMatchElCoef,np.delete(self.Variables.obs.elCoef['name'],
                   matchElCoefInd[:,1]) ))

        matchVelCoef = []
        matchVelCoefInd = []
        try:
	    for i1, key1 in enumerate(self.Variables.sim.velCoef['name']):
	        for i2, key2 in enumerate(self.Variables.obs.velCoef['name']):
	            if key1 == key2:
		        matchVelCoefInd.append((i1,i2))
		        matchVelCoef.append(key1)
	    matchVelCoefInd=np.array(matchVelCoefInd)
	    noMatchVelCoef = np.delete(self.Variables.sim.velCoef['name'],
				       matchVelCoefInd[:,0])
	    np.hstack((noMatchVelCoef,np.delete(self.Variables.obs.velCoef['name'],
                       matchVelCoefInd[:,1]) ))
        except AttributeError:
            pass


        #Compare obs. vs. sim. elevation harmo coef
        data = {}
        columns = ['A', 'g', 'A_ci', 'g_ci']

        #Store harmonics in csv files 
        if save_csv:
            #observed elevation coefs
            for key in columns:
                data[key] = self.Variables.obs.elCoef[key]           
            table = pd.DataFrame(data=data, index=self.Variables.obs.elCoef['name'],
                                 columns=columns)
            ##export as .csv file
            out_file = '{}_obs_el_harmo_coef.csv'.format(filename)
            table.to_csv(out_file)            
            data = {}

            #modeled elevation coefs
            for key in columns:
                data[key] = self.Variables.sim.elCoef[key]           
            table = pd.DataFrame(data=data, index=self.Variables.sim.elCoef['name'],
                                 columns=columns)
            ##export as .csv file
            out_file = '{}_sim_el_harmo_coef.csv'.format(filename)
            table.to_csv(out_file)            
            data = {}

        ##error in %
        if not matchElCoef==[]:
            for key in columns:
                b=self.Variables.sim.elCoef[key][matchElCoefInd[:,0]]
                a=self.Variables.obs.elCoef[key][matchElCoefInd[:,1]]
                err = abs((a-b)/a) * 100.0
                data[key] = err

            ##create table
            table = pd.DataFrame(data=data, index=matchElCoef, columns=columns)        
            ##export as .csv file
            out_file = '{}_el_harmo_error.csv'.format(filename)
            table.to_csv(out_file)
            ##print non-matching coefs
            if not noMatchElCoef.shape[0]==0:
                print "Non-matching harmonic coefficients for elevation: ", noMatchElCoef
        else:
            print "-No matching harmonic coefficients for elevation-" 

        #Compare obs. vs. sim. velocity harmo coef
        data = {}
        columns = ['Lsmaj', 'g', 'theta_ci', 'Lsmin_ci',
                   'Lsmaj_ci', 'theta', 'g_ci']
 
        #Store harmonics in csv files 
        if save_csv:
            #observed elevation coefs
            for key in columns:
                data[key] = self.Variables.obs.velCoef[key]          
            table = pd.DataFrame(data=data, index=self.Variables.obs.velCoef['name'],
                                 columns=columns)
            ##export as .csv file
            out_file = '{}_obs_velo_harmo_coef.csv'.format(filename)
            table.to_csv(out_file)            
            data = {}

            #modeled elevation coefs
            for key in columns:
                data[key] = self.Variables.sim.velCoef[key]           
            table = pd.DataFrame(data=data, index=self.Variables.sim.velCoef['name'],
                                 columns=columns)
            ##export as .csv file
            out_file = '{}_sim_velo_harmo_coef.csv'.format(filename)
            table.to_csv(out_file)            
            data = {}

        ##error in %
        if not matchVelCoef==[]:
            for key in columns:
                b=self.Variables.sim.velCoef[key][matchVelCoefInd[:,0]]
                a=self.Variables.obs.velCoef[key][matchVelCoefInd[:,1]]
                err = abs((a-b)/a) * 100.0
                data[key] = err

            ##create table
            table = pd.DataFrame(data=data, index=matchVelCoef, columns=columns)        
            ##export as .csv file
            out_file = '{}_vel0_harmo_error.csv'.format(filename)
            table.to_csv(out_file)
            ##print non-matching coefs
            if not noMatchVelCoef.shape[0]==0:
                print "Non-matching harmonic coefficients for velocity: ", noMatchVelCoef
        else:
            print "-No matching harmonic coefficients for velocity-"      
    

    def Save_as(self, filename, fileformat='pickle', debug=False):
        """
        This method saves the current FVCOM structure as:
           - *.p, i.e. python file
           - *.mat, i.e. Matlab file

        Inputs:
          - filename = path + name of the file to be saved, string

        Keywords:
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
            try:
                data['Benchmarks'] = self.Benchmarks
            except AttributeError:
                pass 
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
            except (SystemError, MemoryError) as e:
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
