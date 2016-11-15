#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

import numpy as np
import pandas as pd
import cPickle as pkl
from os import makedirs
from os.path import exists

# Polar Plots
from math import *
from cmath import *

# Quick fix
from scipy.io import savemat
from utide import solve

# Local import
from compareData import *
from valTable import valTable
from variablesValidation import _load_validation
from pyseidon_dvt.utilities.interpolation_utils import *

# Local import
from plotsValidation import taylorDiagram, benchmarksMap
from valReport import write_report
from pyseidon_dvt.utilities.miscellaneous import mattime_to_datetime

# Custom error
from pyseidon_dvt.utilities.pyseidon_error import PyseidonError


class Validation:
    """
    **Validation class/structure**

    Class structured as follows: ::

                   _History = Quality Control metadata
                  |_Variables. = observed and simulated variables and quantities
                  |_validate_data = validation method/function against timeseries
      Validation._|_validate_harmonics = validation method/function against
                  |                      harmonic coefficients
                  |_Save_as = "save as" function

    Inputs:
      - observed = standalone or list of PySeidon measurement object (i.e. ADCP, TideGauge, Drifter,...)
      - simulated = standalone or list of any PySeidon simulation object (i.e. FVCOM or Station)
    Option:
      - flow = impose flow comparison by surface flow ('sf'), depth-averaged flow ('daf') or at any depth (float)
               , if negative = from sea surface downwards, if positive = from sea bottom upwards
      - nn = if True then use the nearest location in the grid if the location is outside the grid.
      - outpath = specify a path to save validation results (default = './')
    """
    def __init__(self, observed, simulated, flow=[], nn=True, outpath='./', debug=False, debug_plot=False):
        self._debug = debug
        self._flow = flow
        self._coordinates = []
        self._fig = None
        self._ax = None
        if type(observed) in [tuple, list]:
            self._multi_meas = True
        else:
            self._multi_meas = False
        if type(simulated) in [tuple, list]:
            self._multi_sim = True
        else:
            self._multi_sim = False
        self._debug_plot = debug_plot
        if debug: print '-Debug mode on-'   
        self._nn=nn    
        if debug and nn: print '-Using nearest neighbour-'        
        # creates folder to store outputs
        outpath=outpath.replace(" ","_")
        if outpath[-1] is not '/':
            self._outpath = outpath+'/'
        else:
            self._outpath = outpath
        if not self._outpath == './':
            while exists(self._outpath):
                self._outpath = self._outpath[:-1] + '_bis/'

        #Metadata
        if not self._multi_sim:
            sim_origin = simulated._origin_file
        else:
            sim_origin = 'multiple simulation objects'
        if not self._multi_meas:
            self.History = ['Created from ' + observed._origin_file +\
                            ' and ' + sim_origin]
        else:
            self.History = ['Created from multiple measurement sources' +\
                            ' and ' + sim_origin]
        if (not self._multi_meas) and (not self._multi_sim):
            self._observed = observed
            self._simulated = simulated
            self.Variables = _load_validation(self._outpath, self._observed, self._simulated, flow=self._flow, nn=self._nn, debug=self._debug)
            self._coordinates.append([np.mean(self.Variables.obs.lon), np.mean(self.Variables.obs.lat), self.Variables._obstype])

            # -Append message to History field
	    start = mattime_to_datetime(self.Variables.obs.matlabTime[self.Variables._c[0]])
            end = mattime_to_datetime(self.Variables.obs.matlabTime[self.Variables._c[-1]])
            text = 'Temporal domain from ' + str(start) + ' to ' + str(end)
            self.History.append(text)
        elif (self._multi_meas) and (not self._multi_sim):
            self._observed = observed
            self._simulated = [simulated]
        elif (not self._multi_meas) and (self._multi_sim):
            self._observed = [observed]
            self._simulated = simulated
        else:
            self._observed = observed
            self._simulated = simulated

        # Creates folder once compatibility test passed
        if not self._outpath == './':
            makedirs(self._outpath)
            if debug: print '-Saving results to {}-'.format(self._outpath)
        return

    def _validate_data(self, filename=[], depth=[], slack_velo=0.1,
                       plot=False,  save_csv=False, phase_shift=False, harmo_reconstruct=False,
                       debug=False, debug_plot=False):
        """
        This method computes series of standard validation benchmarks.

        Options:
          - filename = file name of the .csv file to be saved, string.
          - depth = depth at which the validation will be performed, float.
                   Only applicable for 3D simulations.
          - slack_velo = slack water's velocity (m/s), float, everything below will be dumped out
          - plot = plot series of validation graphs, boolean.
                       as well as associated plots in specific folder
          - phase_shift = applies phase shift correction to model quantities
          - harmo_reconstruct = uses harmonic reconstruction is there is no matching period between
                                measurements and simulation

        *References*
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
        # User input
        if filename == []:
            filename = raw_input('Enter filename for csv file: ')
            filename = str(filename)
        if type(self._flow) == float:
            depth = self._flow
        if (depth == [] and self.Variables._3D):
            depth = input('Depth from surface at which the validation will be performed: ')
            depth = float(depth)
        if depth == []:
            depth = 5.0

        # Harmonically reconstruct simulation properties at the original observed matlabTime if harmo is on
        if self.Variables.harmo['On'] and harmo_reconstruct:
            observed = self.Variables.harmo['Observed']
            simulated = self.Variables.harmo['Simulated']
            # Harmonic Analysis
            if self.Variables._simtype == 'station':
                harmo_simEl = simulated.Util2D.Harmonic_analysis_at_point(self.Variables.harmo['nameSite'], elevation=True, velocity=False)
                harmo_simVel = simulated.Util2D.Harmonic_analysis_at_point(self.Variables.harmo['nameSite'], elevation=False, velocity=True)
            elif self.Variables._simtype == 'fvcom':
                harmo_simEl = simulated.Util2D.Harmonic_analysis_at_point(observed.Variables.lon, observed.Variables.lat, elevation=True, velocity=False)
                harmo_simVel = simulated.Util2D.Harmonic_analysis_at_point(observed.Variables.lon, observed.Variables.lat, elevation=False, velocity=True)
            else:
                raise PyseidonError("--Harmonic Analysis not possible with this type of measurement--")
            # Reconstruction at recon_time
            simEl = simulated.Util2D.Harmonic_reconstruction(harmo_simEl, recon_time=observed.Variables.matlabTime)
            simVel = simulated.Util2D.Harmonic_reconstruction(harmo_simVel, recon_time=observed.Variables.matlabTime)
            # Pass reconstructed timeseries to variables structure for validation
            self.Variables.struct['mod_timeseries']['ua'] = simVel['u']
            self.Variables.struct['mod_timeseries']['va'] = simVel['v']
            self.Variables.struct['mod_timeseries']['el'] = simEl['h']
            self.Variables.struct['mod_time'] = observed.Variables.matlabTime
        else:
            # raise PyseidonError("-no matching period. Use 'harmo_reconstruct' option-")
            print "-no matching period. Use 'harmo_reconstruct' option-"

        #initialisation
        vars = []
        self.Suites={}
        threeD = self.Variables.sim._3D
        if self._flow == 'daf': threeD = False

        if self.Variables.struct['type'] == 'ADCP':
            suites = compareOBS(self.Variables.struct, self.Variables._save_path, threeD,
                                    plot=plot, depth=depth, slack_velo=slack_velo, save_csv=save_csv,
                                    phase_shift=phase_shift, debug=debug, debug_plot=debug_plot)
                                    
            for key in suites:
                self.Suites[key] = suites[key]
                vars.append(key)

        elif self.Variables.struct['type'] == 'TideGauge':
            suites = compareOBS(self.Variables.struct, self.Variables._save_path,
                                      plot=plot, slack_velo=slack_velo, save_csv=save_csv,
                                      phase_shift=phase_shift, debug=debug, debug_plot=debug_plot)
            for key in suites:
                self.Suites[key] = suites[key]
                vars.append(key)

        elif self.Variables.struct['type'] == 'Drifter':
            suites = compareOBS(self.Variables.struct, self.Variables._save_path, self.Variables._3D,
                                    depth=depth, plot=plot, slack_velo=slack_velo, save_csv=save_csv,
                                    phase_shift=phase_shift, debug=debug, debug_plot=debug_plot)

            for key in suites:
                self.Suites[key] = suites[key]
                vars.append(key)

        else:
            raise PyseidonError("-This kind of measurements is not supported yet-")

        # Make csv file
        self._Benchmarks = valTable(self.Variables.struct, self.Suites, self.Variables._save_path, filename,  vars,
                                    save_csv=save_csv, debug=debug, debug_plot=debug_plot)
	        
    # Display csv
        print "---Validation benchmarks---"
        pd.set_option('display.max_rows', len(self._Benchmarks))
        print(self._Benchmarks)
        pd.reset_option('display.max_rows')

    # Reload _load_validation to return original values
	# Is this even needed? -Jeremy
	#print '-- Reloading Validation Variables --'
	#self.Variables = _load_validation(self._outpath, self._observed, self._simulated, flow=self._flow, nn=self._nn, debug=self._debug)	

    def _validate_harmonics(self, filename='', save_csv=False, debug=False, debug_plot=False):
        """
        This method computes and store in a csv file the error in %
        for each component of the harmonic analysis (i.e. *_error.csv).

        Options:
          filename: file name of the .csv file to be saved, string.
          save_csv: will save both observed and modeled harmonic
                    coefficients into *.csv files (i.e. *_harmo_coef.csv)
        """
        # check if measurement object is a Drifter
        if self.Variables.struct['type'] == 'Drifter':
            print "--- Harmonic analysis does not work with Drifter's data ---"
            return
        # define attributes
        if not hasattr(self, "_HarmonicBenchmarks"):
            self._HarmonicBenchmarks = HarmonicBenchmarks()
        else:
            delattr(self, '_HarmonicBenchmarks')
            self._HarmonicBenchmarks = HarmonicBenchmarks()
        # User input
        if filename==[]:
            filename = raw_input('Enter filename for csv file: ')
            filename = str(filename)
            
        hasEL=False
        hasUV=False
        commonlist_data = self.Variables.struct['_commonlist_data']
        if 'el' in commonlist_data:
            hasEL=True

        ulist=[var for var in ['ua', 'u'] if var in commonlist_data ]
        vlist=[var for var in ['va', 'v'] if var in commonlist_data ]
        if len(ulist)>0 and len(vlist)>0:
            hasUV=True

        # Harmonic analysis over matching & non-matching time
        obs_time = self.Variables.struct['obs_time']
        obs_lat = self.Variables.struct['obs_lat']
        mod_time = self.Variables.struct['mod_time']
        mod_lat = self.Variables.struct['mod_lat']
        if hasEL:     
            obs_el =  self.Variables.struct['obs_timeseries']['el'][:]
            mod_el =  self.Variables.struct['mod_timeseries']['el'][:]
            self.Variables.obs.elCoef = solve(obs_time, obs_el, None, obs_lat,
                                         constit='auto', trend=False, Rayleigh_min=0.95,
                                         method='ols', conf_int='linear')
            self.Variables.sim.elCoef = solve(mod_time, mod_el, None, mod_lat,
                                         constit='auto', trend=False, Rayleigh_min=0.95,
                                         method='ols', conf_int='linear')

        if hasUV:
            obs_ua =  self.Variables.struct['obs_timeseries']['ua'][:]
            obs_va =  self.Variables.struct['obs_timeseries']['va'][:]
            mod_ua =  self.Variables.struct['mod_timeseries']['ua'][:]
            mod_va =  self.Variables.struct['mod_timeseries']['va'][:]
            self.Variables.obs.velCoef = solve(obs_time, obs_ua, obs_va, obs_lat,
                                         constit='auto', trend=False, Rayleigh_min=0.95,
                                         method='ols', conf_int='linear')
            self.Variables.sim.velCoef = solve(mod_time, mod_ua, mod_va, mod_lat,
                                         constit='auto', trend=False, Rayleigh_min=0.95,
                                         method='ols', conf_int='linear')


        # find matching and non-matching coef & hold their magnitude and phase for comparison
        matchElCoef = [] # Coefficient Names
        matchElCoefInd = [] # Coefficient Names' Indexes
	matchEl_Adiff = [] # Coefficents' Maginitude Difference
	matchEl_gdiff = [] # Coefficients' Phase Difference
        for i1, key1 in enumerate(self.Variables.sim.elCoef['name']):
            for i2, key2 in enumerate(self.Variables.obs.elCoef['name']):
                if key1 == key2:
                   matchElCoefInd.append((i1,i2))
                   matchElCoef.append(key1)
		   matchEl_Adiff.append(abs(self.Variables.sim.elCoef['A'][i1]*exp(1j*radians(self.Variables.sim.elCoef['g'][i1])) - 
				       	    self.Variables.obs.elCoef['A'][i2]*exp(1j*radians(self.Variables.obs.elCoef['g'][i2]))))
		   matchEl_gdiff.append(abs(self.Variables.sim.elCoef['g'][i1] - self.Variables.obs.elCoef['g'][i2]))
        matchElCoefInd=np.array(matchElCoefInd)
        noMatchElCoef = np.delete(self.Variables.sim.elCoef['name'], matchElCoefInd[:,0])
        np.hstack((noMatchElCoef, np.delete(self.Variables.obs.elCoef['name'], matchElCoefInd[:,1])))
	
	# Repeat for Velocity Coefficients
        matchVelCoef = []
        matchVelCoefInd = []
	matchVel_LsmajDiff = []
	matchVel_LsminDiff = []
	matchVel_gdiff = []
        try:
            for i1, key1 in enumerate(self.Variables.sim.velCoef['name']):
                for i2, key2 in enumerate(self.Variables.obs.velCoef['name']):
                    if key1 == key2:
                        matchVelCoefInd.append((i1, i2))
                        matchVelCoef.append(key1)
			matchVel_LsmajDiff.append(abs(self.Variables.sim.velCoef['Lsmaj'][i1]*exp(1j*radians(self.Variables.sim.velCoef['g'][i1])) - 
						  self.Variables.obs.velCoef['Lsmaj'][i2]*exp(1j*radians(self.Variables.obs.velCoef['g'][i2]))))
		   	matchVel_LsminDiff.append(abs(self.Variables.sim.velCoef['Lsmin'][i1]*exp(1j*radians(self.Variables.sim.velCoef['g'][i1])) - 
						  self.Variables.obs.velCoef['Lsmin'][i2]*exp(1j*radians(self.Variables.obs.velCoef['g'][i2]))))
			matchVel_gdiff.append(abs(self.Variables.sim.velCoef['g'][i1] - self.Variables.obs.velCoef['g'][i2]))
            matchVelCoefInd = np.array(matchVelCoefInd)
            noMatchVelCoef = np.delete(self.Variables.sim.velCoef['name'], matchVelCoefInd[:, 0])
            np.hstack((noMatchVelCoef, np.delete(self.Variables.obs.velCoef['name'], matchVelCoefInd[:, 1])))
        except AttributeError:
            pass

	# Compare largest ten obs. vs sim. coefficients visually on Polar plots
	# I would prefer to have these plots saved to the directory created by validate_harmonics()
	top10_el = zip(matchEl_Adiff, matchElCoef, matchEl_gdiff)
	top10_el.sort()
	top10_el = top10_el[-10:]
	top10_el = zip(*top10_el)
	plt.axes(polar=True)
	plt.plot(top10_el[2], top10_el[0], 'r.')
	for i, txt in enumerate(top10_el[1]):
	    plt.annotate(top10_el[1][i], (str(top10_el[2][i]), str(top10_el[0][i])), size=8)
	plt.title('Harmonic Analysis Elevation Coefficient comparison')
	plt.savefig('HarmoEl_coeffcompare', format='png')
	plt.clf()

	top10_vel_Lsmaj = zip(matchVel_LsmajDiff, matchVelCoef, matchVel_gdiff)
	top10_vel_Lsmaj.sort()
	top10_vel_Lsmaj = top10_vel_Lsmaj[-10:]
	top10_vel_Lsmaj = zip(*top10_vel_Lsmaj)
	plt.axes(polar=True)
	plt.plot(top10_vel_Lsmaj[2], top10_vel_Lsmaj[0], 'r.')
	for i, txt in enumerate(top10_vel_Lsmaj[1]):
	    plt.annotate(top10_vel_Lsmaj[1][i], (str(top10_vel_Lsmaj[2][i]), str(top10_vel_Lsmaj[0][i])), size=8)
	plt.title('Harmonic Analysis Lsmaj Velocity Coefficient comparison')
	plt.savefig('HarmoVel_Lsmaj_coeffcompare', format='png')
	plt.clf()

	top10_vel_Lsmin = zip(matchVel_LsminDiff, matchVelCoef, matchVel_gdiff)
	top10_vel_Lsmin.sort()
	top10_vel_Lsmin = top10_vel_Lsmin[-10:]
	top10_vel_Lsmin = zip(*top10_vel_Lsmin)	
	plt.axes(polar=True)
	plt.plot(top10_vel_Lsmin[2], top10_vel_Lsmin[0], 'r.')
	for i, txt in enumerate(top10_vel_Lsmin[1]):
	    plt.annotate(top10_vel_Lsmin[1][i], (str(top10_vel_Lsmin[2][i]), str(top10_vel_Lsmin[0][i])), size=8)
	plt.title('Harmonic Analysis Lsmin Velocity Coefficient comparison')
	plt.savefig('HarmoVel_Lsmin_coeffcompare', format='png')
	plt.clf()

        # Compare obs. vs. sim. elevation harmo coef
        data = {}
        columns = ['A', 'g', 'A_ci', 'g_ci']
        measname = self.Variables.struct['name'].split('/')[-1].split('.')[0]

        # Store harmonics in csv files
        if save_csv and hasEL:
            # observed elevation coefs
            for key in columns:
                data[key] = self.Variables.obs.elCoef[key]
            data['constituent'] = self.Variables.obs.elCoef['name']
            measlist = [measname] * len(data['constituent'])
            table = pd.DataFrame(data=data, index=measlist,
                                 columns=columns + ['constituent'])
            # export as .csv file
            out_file = '{}{}_obs_el_harmo_coef.csv'.format(self.Variables._save_path, filename)
            table.to_csv(out_file)
            data = {}

            #modeled elevation coefs
            for key in columns:
                data[key] = self.Variables.sim.elCoef[key]
            data['constituent'] = self.Variables.sim.elCoef['name']
            measlist = [measname] * len(data['constituent'])
            table = pd.DataFrame(data=data, index=measlist,
                                 columns=columns + ['constituent'])
            # export as .csv file
            out_file = '{}{}_sim_el_harmo_coef.csv'.format(self.Variables._save_path, filename)
            table.to_csv(out_file)
            data = {}

        # error in %
        if hasEL:
            if not matchElCoef==[]:
                for key in columns:
                    b=self.Variables.sim.elCoef[key][matchElCoefInd[:,0]]
                    a=self.Variables.obs.elCoef[key][matchElCoefInd[:,1]]
                    err = abs((a-b)/a) * 100.0
                    data[key] = err
                data['constituent'] = matchElCoef
                measlist = [measname] * len(matchElCoef)
                ##create table
                table = pd.DataFrame(data=data, index=measlist, columns=columns + ['constituent'])
                ##export as .csv file
                out_file = '{}{}_el_harmo_error.csv'.format(self.Variables._save_path, filename)
                table.to_csv(out_file)
                ##print non-matching coefs
                if not noMatchElCoef.shape[0]==0:
                    print "Non-matching harmonic coefficients for elevation: ", noMatchElCoef
            else:
                print "-No matching harmonic coefficients for elevation-"

            # save dataframe in attribute
            self._HarmonicBenchmarks.elevation = table

        #Compare obs. vs. sim. velocity harmo coef
        data = {}
        columns = ['Lsmaj', 'g', 'theta_ci', 'Lsmin_ci',
                   'Lsmaj_ci', 'theta', 'g_ci']

        #Store harmonics in csv files
        if save_csv and hasUV:
            #observed elevation coefs
            for key in columns:
                data[key] = self.Variables.obs.velCoef[key]
            data['constituent'] = matchVelCoef
            measlist = [measname] * len(data['constituent'])
            table = pd.DataFrame(data=data, index=measlist,
                                 columns=columns + ['constituent'])
            ##export as .csv file
            out_file = '{}{}_obs_velo_harmo_coef.csv'.format(self.Variables._save_path, filename)
            table.to_csv(out_file)
            data = {}

            #modeled elevation coefs
            for key in columns:
                data[key] = self.Variables.sim.velCoef[key]
            data['constituent'] = matchVelCoef
            measlist = [measname] * len(data['constituent'])
            table = pd.DataFrame(data=data, index=measlist,
                                 columns=columns + ['constituent'])
            ##export as .csv file
            out_file = '{}{}_sim_velo_harmo_coef.csv'.format(self.Variables._save_path, filename)
            table.to_csv(out_file)
            data = {}

        ##error in %
        if hasUV:
            if not matchVelCoef==[]:
                for key in columns:
                    b=self.Variables.sim.velCoef[key][matchVelCoefInd[:,0]]
                    a=self.Variables.obs.velCoef[key][matchVelCoefInd[:,1]]
                    err = abs((a-b)/a) * 100.0
                    data[key] = err
                data['constituent'] = matchVelCoef
                measlist = [measname] * len(matchVelCoef)
                ##create table
                table = pd.DataFrame(data=data, index=measlist, columns=columns + ['constituent'])
                ##export as .csv file
                out_file = '{}{}_vel0_harmo_error.csv'.format(self.Variables._save_path, filename)
                table.to_csv(out_file)
                ##print non-matching coefs
                if not noMatchVelCoef.shape[0]==0:
                    print "Non-matching harmonic coefficients for velocity: ", noMatchVelCoef
            else:
                print "-No matching harmonic coefficients for velocity-"

            # save dataframe in attribute
            self._HarmonicBenchmarks.velocity = table

    def validate_data(self, filename=[], depth=[], slack_velo=0.1, plot=False, save_csv=False, phase_shift=False,
                      debug=False, debug_plot=False):
        """
        This method computes series of standard validation benchmarks.

        Options:
          - filename = file name of the .csv file to be saved, string.
          - depth = depth at which the validation will be performed, float.
                   Only applicable for 3D simulations.
                   If negative = from sea surface downwards, if positive = from sea bottom upwards
          - slack_velo = slack water's velocity (m/s), float, everything below will be dumped out
          - plot = plot series of validation graphs, boolean.
          - save_csv = will save benchmark values into *.csv file
                       as well as associated plots in specific folder
          - phase_shift = applies phase shift correction to model quantities

        *References*
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
        if (not self._multi_meas) and (not self._multi_sim):
            self._validate_data(filename, depth, slack_velo, plot, save_csv, phase_shift, debug, debug_plot)
            self.Benchmarks = self._Benchmarks
        else:
            I=0
            for sim in self._simulated:
                for meas in self._observed:
                    #try:
                        self.Variables = _load_validation(self._outpath, meas, sim, flow=self._flow, debug=self._debug)

                        # -Append message to History field
                        start = mattime_to_datetime(self.Variables.obs.matlabTime[self.Variables._c[0]])
                        end = mattime_to_datetime(self.Variables.obs.matlabTime[self.Variables._c[-1]])
                        text = 'Temporal domain from ' + str(start) + ' to ' + str(end)
                        try:
                            self.History[1] = text
                        except IndexError:
                            self.History.append(text)

                        self._validate_data(filename, depth, slack_velo, plot, save_csv, phase_shift, debug, debug_plot)
                        if I == 0:
                            self.Benchmarks = self._Benchmarks
                            I += 1
                        else:
                            self.Benchmarks = pd.concat([self.Benchmarks, self._Benchmarks])

                        self._coordinates.append([np.mean(self.Variables.obs.lon),
                                                  np.mean(self.Variables.obs.lat),
                                                  self.Variables._obstype])
                    #except PyseidonError:
                    #except:  # making it even more permissive
                    #    print "Error with measurement object "+meas.History[0]
                    #    continue
        if save_csv:
            #if self._multi_meas:
            #    savepath = self.Variables._save_path[:(self.Variables._save_path[:-1].rfind('/')+1)]
            #else:
            #    savepath = self.Variables._save_path
            savepath = self._outpath
            #savepath = self.Variables._save_path
            try:
                out_file = '{}{}_benchmarks.csv'.format(savepath, filename)
                self.Benchmarks.to_csv(out_file)
            except AttributeError:
                raise PyseidonError("-No matching measurement-")

    def validate_harmonics(self, filename=[], save_csv=False, debug=False, debug_plot=False):
        """
        This method computes and store in a csv file the error in %
        for each component of the harmonic analysis (i.e. *_error.csv).

        Options:
          filename: file name of the .csv file to be saved, string.
          save_csv: will save both observed and modeled harmonic
                    coefficients into *.csv files (i.e. *_harmo_coef.csv)
        """
        if (not self._multi_meas) and (not self._multi_sim):
            self.Variables = _load_validation(self._outpath, self._observed, self._simulated, flow=self._flow, debug=self._debug)
            self._validate_harmonics(filename, save_csv, debug, debug_plot)
            self.HarmonicBenchmarks = self._HarmonicBenchmarks
        else:
            I=0
            J=0
            for sim in self._simulated:
                for meas in self._observed:
                    try:
                        self.Variables = _load_validation(self._outpath, meas, sim, flow=self._flow, debug=self._debug)
                        self._validate_harmonics(filename, save_csv, debug, debug_plot)
                        if I == 0 and J == 0:
                            self.HarmonicBenchmarks = HarmonicBenchmarks()
                        if I == 0 and type(self._HarmonicBenchmarks.elevation) != str:
                            self.HarmonicBenchmarks.elevation = self._HarmonicBenchmarks.elevation
                            I += 1
                        elif I != 0 and type(self._HarmonicBenchmarks.elevation) != str:
                            self.HarmonicBenchmarks.elevation = pd.concat([self.HarmonicBenchmarks.elevation,
                                                                           self._HarmonicBenchmarks.elevation])
                        if J == 0 and type(self._HarmonicBenchmarks.velocity) != str:
                            self.HarmonicBenchmarks.velocity = self._HarmonicBenchmarks.velocity
                            J += 1
                        elif J != 0 and type(self._HarmonicBenchmarks.velocity) != str:
                            self.HarmonicBenchmarks.velocity = pd.concat([self.HarmonicBenchmarks.velocity,
                                                                          self._HarmonicBenchmarks.velocity])
                    except PyseidonError:
                        pass

        # Drop duplicated lines
        self.HarmonicBenchmarks.velocity.drop_duplicates(inplace=True)
        self.HarmonicBenchmarks.elevation.drop_duplicates(inplace=True)

        if save_csv:
            #if self._multi_meas:
            #    savepath = self.Variables._save_path[:(self.Variables._save_path[:-1].rfind('/')+1)]
            #else:
            #    savepath = self.Variables._save_path
            savepath = self._outpath
            try:
                try:
                    out_file = '{}{}_elevation_harmonic_benchmarks.csv'.format(savepath, filename)
                    self.HarmonicBenchmarks.elevation.to_csv(out_file)
                except AttributeError:
                    pass
                try:
                    out_file = '{}{}_velocity_harmonic_benchmarks.csv'.format(savepath, filename)
                    self.HarmonicBenchmarks.velocity.to_csv(out_file)
                except AttributeError:
                    pass
            except AttributeError:
                raise PyseidonError("-No matching measurement-")

    def taylor_diagram(self, savepath='', fname="taylor_diagram", labels=True, debug=False):
        """
        Plots Taylor diagram based on the results of 'validate_data'

        Options:
          - savepath = folder path for saving plot, string
          - fname = filename for saving plot, string
          - labels = labels and legend, boolean
        """
        try:
            self._fig, self._ax = taylorDiagram(self.Benchmarks,
                                                savepath=savepath, fname=fname, labels=labels, debug=debug)
        except AttributeError:
            raise PyseidonError("-validate_data needs to be run first-")

    def benchmarks_map(self, savepath='', fname="benchmarks_map", debug=False):
        """
    Plots bathymetric map & model validation benchmarks

    Options:
      - savepath = folder path for saving plot, string
      - fname = filename for saving plot, string

    Note: this function shall work only if ADCP object(s) and FVCOM object
          have been used as inputs
        """
        if not self._simulated.__module__.split('.')[-1] == 'fvcomClass':
            raise PyseidonError("---work only with a combination ADCP object(s) and FVCOM object---")
        try:
            benchmarksMap(self.Benchmarks, self._observed, self._simulated, savepath=savepath, fname=fname, debug=debug)
        except AttributeError:
            raise PyseidonError("---validate_data needs to be run first---")

    def save_as(self, filename, fileformat='pickle', debug=False):
        """
        This method saves the current Validation structure as:
           - *.p, i.e. python file
           - *.mat, i.e. Matlab file

        Inputs:
          - filename = path + name of the file to be saved, string

        Options:
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
            #TR comment: based on Mitchell O'Flaherty-Sproul's code
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
                if any([type(Var[key]).__name__ == x for x in listkeys]):
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

    def write_validation_report(self, report_title="validation_report.pdf", debug=False):
        """
        This method writes a report (*.pdf) based on the validation methods' results

        Kwargs:
          report_title (str): file name
          debug (bool): debug flag
        """
        debug = debug or self._debug
        write_report(self, report_title=report_title, debug=debug)

# utility classes
class HarmonicBenchmarks:
    """
    Storage for hamonic benchmarks
    """
    def __init__(self):
        self.elevation = 'No harmonic benchmarks for the elevation yet. Run validate_harmonics'
        self.velocity = 'No harmonic benchmarks for the velocity yet. Run validate_harmonics'
        return
