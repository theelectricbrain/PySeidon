#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division

import numpy as np
from os import mkdir
from os.path import exists
from tidalStats import TidalStats
from smooth import smooth
from depthInterp import depthFromSurf
from datetime import datetime, timedelta

# Custom error
from pyseidon.utilities.pyseidon_error import PyseidonError

# Local import
from plotsValidation import *

def dn2dt(datenum):
    '''
    Convert matlab datenum to python datetime.
    '''
    return datetime.fromordinal(int(datenum)) + timedelta(days=datenum%1) - \
           timedelta(days=366)

def compareUV(data, threeDim, depth=5, plot=False, save_csv=False,
              debug=False, debug_plot=False):
    """
    Does a comprehensive validation process between modeled and observed.
    Outputs a list of important statistics for each variable, calculated
    using the TidalStats class

    Inputs:
        - data = dictionary containing all necessary observed and model data
        - threeDim = boolean flag, 3D or not
    Outputs:
       - elev_suite = dictionary of useful statistics for sea elevation
       - speed_suite = dictionary of useful statistics for flow speed
       - dir_suite = dictionary of useful statistics for flow direction
       - u_suite = dictionary of useful statistics for u velocity component
       - v_suite = dictionary of useful statistics for v velocity component
       - vel_suite = dictionary of useful statistics for signed flow velocity
       - csp_suite = dictionary of useful statistics for cubic flow speed
    Options:
       - depth = interpolation depth from surface, float
       - plot = boolean flag for plotting results
       - save_csv = boolean flag for saving statistical benchmarks in csv file
    """
    if debug: print "CompareUV..."
    # take data from input dictionary
    mod_time = data['mod_time']
    if not data['type'] == 'Drifter':
        obs_time = data['obs_time']

        mod_el = data['mod_timeseries']['elev']
        obs_el = data['obs_timeseries']['elev']
    else:
        obs_time = data['mod_time']

    # Save path
    name = data['name']
    save_path = name.split('/')[-1].split('.')[0]+'/'

    # Check if 3D simulation
    if threeDim:
        obs_u_all = data['obs_timeseries']['u']
        obs_v_all = data['obs_timeseries']['v']
        mod_u_all = data['mod_timeseries']['u']
        mod_v_all = data['mod_timeseries']['v']
        bins = data['obs_timeseries']['bins']
        siglay = data['mod_timeseries']['siglay']
        # use depth interpolation to get a single timeseries
        mod_depth = mod_el + np.mean(obs_el[~np.isnan(obs_el)])
        (mod_u, obs_u) = depthFromSurf(mod_u_all, mod_depth, siglay,
                                       obs_u_all, obs_el, bins, depth=depth,
                                       debug=debug, debug_plot=debug_plot)
        (mod_v, obs_v) = depthFromSurf(mod_v_all, mod_depth, siglay,
                                       obs_v_all, obs_el, bins, depth=depth,
                                       debug=debug, debug_plot=debug_plot)
    else:
        if not data['type'] == 'Drifter':
            obs_u = data['obs_timeseries']['ua']
            obs_v = data['obs_timeseries']['va']
            mod_u = data['mod_timeseries']['ua']
            mod_v = data['mod_timeseries']['va']
        else:
            obs_u = data['obs_timeseries']['u']
            obs_v = data['obs_timeseries']['v']
            mod_u = data['mod_timeseries']['u']
            mod_v = data['mod_timeseries']['v']


    if debug: print "...convert times to datetime..."
    mod_dt, obs_dt = [], []
    for i in mod_time:
        mod_dt.append(dn2dt(i))
    for j in obs_time:
        obs_dt.append(dn2dt(j))

    if debug: print "...put data into a useful format..."
    mod_spd = np.sqrt(mod_u**2.0 + mod_v**2.0)
    obs_spd = np.sqrt(obs_u**2.0 + obs_v**2.0)
    mod_dir = np.arctan2(mod_v, mod_u) * 180.0 / np.pi
    obs_dir = np.arctan2(obs_v, obs_u) * 180.0 / np.pi
    if not data['type'] == 'Drifter':
        obs_el = obs_el - np.mean(obs_el[~np.isnan(obs_el)])
    # Chose the component with the biggest variance as sign reference
    if np.var(mod_v) > np.var(mod_u):
            mod_signed = np.sign(mod_v)
            obs_signed = np.sign(obs_v)
    else:
            mod_signed = np.sign(mod_u)
            obs_signed = np.sign(obs_u)

    if debug: print "...check if the modeled data lines up with the observed data..."
    if (mod_time[-1] < obs_time[0] or obs_time[-1] < mod_time[0]):
        raise PyseidonError("---time periods do not match up---")

    else:
        if debug: print "...interpolate the data onto a common time step for each data type..."
        if not data['type'] == 'Drifter':
            # elevation
            (mod_el_int, obs_el_int, step_el_int, start_el_int) = smooth(mod_el, mod_dt, obs_el, obs_dt,
                                                                         debug=debug, debug_plot=debug_plot)
            # speed
            (mod_sp_int, obs_sp_int, step_sp_int, start_sp_int) = smooth(mod_spd, mod_dt, obs_spd, obs_dt,
                                                                         debug=debug, debug_plot=debug_plot)
            # direction
            (mod_dr_int, obs_dr_int, step_dr_int, start_dr_int) = smooth(mod_dir, mod_dt, obs_dir, obs_dt,
                                                                         debug=debug, debug_plot=debug_plot)
            # u velocity
            (mod_u_int, obs_u_int, step_u_int, start_u_int) = smooth(mod_u, mod_dt, obs_u, obs_dt,
                                                                     debug=debug, debug_plot=debug_plot)
            # v velocity
            (mod_v_int, obs_v_int, step_v_int, start_v_int) = smooth(mod_v, mod_dt, obs_v, obs_dt,
                                                                     debug=debug, debug_plot=debug_plot)
            # velocity i.e. signed speed
            (mod_ve_int, obs_ve_int, step_ve_int, start_ve_int) = smooth(mod_spd * mod_signed, mod_dt,
                                                                         obs_spd * obs_signed, obs_dt,
                                                                         debug=debug, debug_plot=debug_plot)
            # cubic speed
            #mod_cspd = mod_spd**3.0
            #obs_cspd = obs_spd**3.0
            mod_cspd = mod_signed * mod_spd**3.0
            obs_cspd = mod_signed * obs_spd**3.0
            (mod_cspd_int, obs_cspd_int, step_cspd_int, start_cspd_int) = smooth(mod_cspd, mod_dt, obs_cspd, obs_dt,
                                                                                 debug=debug, debug_plot=debug_plot)
        else:
            # Time steps
            step = mod_time[1] - mod_time[0]
            start = mod_time[0]

            # Already interpolated, so no need to use smooth...
            # speed
            (mod_sp_int, obs_sp_int, step_sp_int, start_sp_int) = (mod_spd, obs_spd, step, start)
            # direction
            (mod_dr_int, obs_dr_int, step_dr_int, start_dr_int) = (mod_dir, obs_dir, step, start)
            # u velocity
            (mod_u_int, obs_u_int, step_u_int, start_u_int) = (mod_u, obs_u, step, start)
            # v velocity
            (mod_v_int, obs_v_int, step_v_int, start_v_int) = (mod_v, obs_v, step, start)
            # velocity i.e. signed speed
            (mod_ve_int, obs_ve_int, step_ve_int, start_ve_int) = (mod_spd, obs_spd, step, start)
            # cubic speed
            mod_cspd = mod_spd**3.0
            obs_cspd = obs_spd**3.0
            (mod_cspd_int, obs_cspd_int, step_cspd_int, start_cspd_int) = (mod_cspd, obs_cspd, step, start)

    if debug: print "...remove directions where velocities are small..."
    MIN_VEL = 0.1
    indexMin = np.where(obs_sp_int < MIN_VEL)
    obs_dr_int[indexMin] = np.nan
    obs_u_int[indexMin] = np.nan
    obs_v_int[indexMin] = np.nan
    obs_ve_int[indexMin] = np.nan
    obs_cspd_int[indexMin] = np.nan

    indexMin = np.where(mod_sp_int < MIN_VEL)
    mod_dr_int[indexMin] = np.nan
    mod_u_int[indexMin] = np.nan
    mod_v_int[indexMin] = np.nan
    mod_ve_int[indexMin] = np.nan
    mod_cspd_int[indexMin] = np.nan

    if debug: print "...get stats for each tidal variable..."
    gear = data['type'] # Type of measurement gear (drifter, adcp,...)
    if not gear == 'Drifter':
        elev_suite = tidalSuite(gear, mod_el_int, obs_el_int, step_el_int, start_el_int,
                                [], [], [], [], [], [],
                                kind='elevation', plot=plot,
                                save_csv=save_csv, save_path=save_path,
                                debug=debug, debug_plot=debug_plot)
    else:
        elev_suite = []
    speed_suite = tidalSuite(gear, mod_sp_int, obs_sp_int, step_sp_int, start_sp_int,
                             [], [], [], [], [], [],
                             kind='speed', plot=plot,
                             save_csv=save_csv, save_path=save_path,
                             debug=debug, debug_plot=debug_plot)
    dir_suite = tidalSuite(gear, mod_dr_int, obs_dr_int, step_dr_int, start_dr_int,
                           [], [], [], [], [], [],
                           kind='direction', plot=plot,
                           save_csv=save_csv, save_path=save_path,
                           debug=debug, debug_plot=debug_plot)
    u_suite = tidalSuite(gear, mod_u_int, obs_u_int, step_u_int, start_u_int,
                         [], [], [], [], [], [],
                         kind='u velocity', plot=plot, save_csv=save_csv, save_path=save_path,
                         debug=debug, debug_plot=debug_plot)
    v_suite = tidalSuite(gear, mod_v_int, obs_v_int, step_v_int, start_v_int,
                         [], [], [], [], [], [],
                         kind='v velocity', plot=plot, save_csv=save_csv, save_path=save_path,
                         debug=debug, debug_plot=debug_plot)

    # TR: requires special treatments from here on
    vel_suite = tidalSuite(gear, mod_ve_int, obs_ve_int, step_ve_int, start_ve_int,
                           mod_u, obs_u, mod_v, obs_v,
                           mod_dt, obs_dt,
                           kind='velocity', plot=plot, save_csv=save_csv, save_path=save_path,
                           debug=debug, debug_plot=debug_plot)
    csp_suite = tidalSuite(gear, mod_cspd_int, obs_cspd_int, step_cspd_int, start_cspd_int,
                           mod_u, obs_u, mod_v, obs_v,
                           mod_dt, obs_dt,
                           kind='cubic speed', plot=plot, save_csv=save_csv, save_path=save_path,
                           debug=debug, debug_plot=debug_plot)

    # output statistics in useful format

    if debug: print "...CompareUV done."

    return (elev_suite, speed_suite, dir_suite, u_suite, v_suite, vel_suite, csp_suite)

def compareTG(data, plot=False, save_csv=False, debug=False, debug_plot=False):
    """
    Does a comprehensive comparison between tide gauge height data and
    modeled data.

    Input:
       - data = dictionary containing all necessary tide gauge and model data.
    Outputs
       - elev_suite = dictionary of useful statistics
    Options:
       - plot = boolean flag for plotting results
       - save_csv = boolean flag for saving statistical benchmarks in csv file
    """
    if debug: print "CompareTG..."
    # load data
    mod_elev = data['mod_timeseries']['elev']
    obs_elev = data['obs_timeseries']['elev']
    obs_datenums = data['obs_time']
    mod_datenums = data['mod_time']
    gear = data['type'] # Type of measurement gear (drifter, adcp,...)
    #TR: comment out
    #mod_harm = data['elev_mod_harmonics']

    # Save path
    name = data['name']
    save_path = name.split('/')[-1].split('.')[0]+'/'

    # convert times and grab values
    obs_time, mod_time = [], []
    for i, v in enumerate(obs_datenums):
        obs_time.append(dn2dt(v))
    for j, w in enumerate(mod_datenums):
        mod_time.append(dn2dt(w))

    if debug: print "...check if they line up in the time domain..."
    if (mod_time[-1] < obs_time[0] or obs_time[-1] < mod_time[0]):
        raise PyseidonError("---time periods do not match up---")

    else:

        if debug: print "...interpolate timeseries onto a common timestep..."
        (mod_elev_int, obs_elev_int, step_int, start_int) = \
            smooth(mod_elev, mod_time, obs_elev, obs_time,
                   debug=debug, debug_plot=debug_plot)

    elev_suite = tidalSuite(gear, mod_elev_int, obs_elev_int, step_int, start_int,
                            [], [], [], [], [], [],
                            kind='elevation', plot=plot, save_csv=save_csv, save_path=save_path,
                            debug=debug, debug_plot=debug_plot)

    if debug: print "...CompareTG done."

    return elev_suite

def tidalSuite(gear, model, observed, step, start,
               model_u, observed_u, model_v, observed_v,
               model_time, observed_time,
               kind='', plot=False, save_csv=False, save_path='./',
               debug=False, debug_plot=False):
    """
    Create stats classes for a given tidal variable.

    Accepts interpolated model and observed data, the timestep, and start
    time. kind is a string representing the kind of data. If plot is set
    to true, a time plot and regression plot will be produced.

    Returns a dictionary containing all the stats.
    """
    if debug: print "tidalSuite..."
    stats = TidalStats(gear, model, observed, step, start,
                       model_u = model_u, observed_u = observed_u, model_v = model_v, observed_v = observed_v,
                       model_time = model_time, observed_time = observed_time,
                       kind=kind, debug=debug, debug_plot=debug_plot)
    stats_suite = stats.getStats()
    stats_suite['r_squared'] = stats.linReg()['r_2']
    try: #Fix for Drifter's data
        stats_suite['phase'] = stats.getPhase(debug=debug)
    except:
        stats_suite['phase'] = 0.0

    if plot or debug_plot:
        plotData(stats)
        plotRegression(stats, stats.linReg())

    if save_csv:
        if not exists(save_path):
            mkdir(save_path)
        else:
            save_path = save_path[:-1] + '_bis/'
        stats.save_data(path=save_path)
        plotData(stats, savepath=save_path, fname=kind+"_"+gear+"_time_series.png")
        plotRegression(stats, stats.linReg(), savepath=save_path, fname=kind+"_"+gear+"_linear_regression.png")

    if debug: print "...tidalSuite done."

    return stats_suite
