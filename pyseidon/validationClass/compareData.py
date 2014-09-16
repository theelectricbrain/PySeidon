import numpy as np
from tidalStats import TidalStats
#from interpolate import interpol
from smooth import smooth
from datetime import datetime, timedelta
from utide import ut_reconstr
import matplotlib.pyplot as plt
from depthInterp import depthFromSurf
from save_FlowFile_BPFormat import sign_speed, get_DirFromN

def dn2dt(datenum):
    '''
    Convert matlab datenum to python datetime.
    '''
    return datetime.fromordinal(int(datenum)) + timedelta(days=datenum%1) - \
           timedelta(days=366)

def compareUV(data):
    '''
    Does a comprehensive validation process between modeled and observed
    data on the following:
        Current speed
        Current direction
        Harmonic constituents (for height and speed)

    Outputs a list of important statistics for each variable, calculated
    using the TidalStats class
    '''
    # take data from input dictionary
    mod_time = data['mod_time']
    obs_time = data['obs_time']
    mod_u_all = data['mod_timeseries']['u']
    mod_v_all = data['mod_timeseries']['v']
    mod_el = data['mod_timeseries']['elev']
    obs_u_all = data['obs_timeseries']['u']
    obs_v_all = data['obs_timeseries']['v']
    obs_el = data['obs_timeseries']['elev']
    v_mod_harm = data['vel_mod_harmonics']
    v_obs_harm = data['vel_obs_harmonics']
    el_mod_harm = data['elev_mod_harmonics']
    el_obs_harm = data['elev_mod_harmonics']
    bins = data['obs_timeseries']['bins']
    siglay = data['mod_timeseries']['siglay']

    # use depth interpolation to get a single timeseries
    mod_depth = mod_el + np.mean(obs_el)
    (mod_u, obs_u) = depthFromSurf(mod_u_all, mod_depth, siglay,
				   obs_u_all, obs_el, bins)
    (mod_v, obs_v) = depthFromSurf(mod_v_all, mod_depth, siglay,
                                   obs_v_all, obs_el, bins)


    # convert times to datetime
    mod_dt, obs_dt = [], []
    for i in mod_time:
	mod_dt.append(dn2dt(i))
    for j in obs_time:
	obs_dt.append(dn2dt(j))

    # put data into a useful format
    mod_spd = np.sqrt(mod_u**2 + mod_v**2)
    obs_spd = np.sqrt(obs_u**2 + obs_v**2)
    mod_dir = np.arctan2(mod_v, mod_u) * 180 / np.pi
    obs_dir = np.arctan2(obs_v, obs_u) * 180 / np.pi
    obs_el = obs_el - np.mean(obs_el)

    # check if the modeled data lines up with the observed data
    if (mod_time[-1] < obs_time[0] or obs_time[-1] < mod_time[0]):

	pred_uv = ut_reconstr(obs_time, v_mod_harm)
	pred_uv = np.asarray(pred_uv)
	pred_h = ut_reconstr(obs_time, el_mod_harm)
	pred_h = np.asarray(pred_h)

	# redo speed and direction and set interpolated variables
	mod_sp_int = np.sqrt(pred_uv[0]**2 + pred_uv[1]**2)
	mod_ve_int = mod_sp_int * np.sign(pred_uv[1])
	mod_dr_int = np.arctan2(pred_uv[1], pred_uv[0]) * 180 / np.pi
	mod_el_int = pred_h[0]
	mod_u_int = pred_uv[0]
	mod_v_int = pred_uv[1]
	obs_sp_int = obs_spd
	obs_ve_int = obs_spd * np.sign(obs_v)
	obs_dr_int = obs_dir
	obs_el_int = obs_el
	obs_u_int = obs_u
	obs_v_int = obs_v
	step_int = obs_dt[1] - obs_dt[0]
	start_int = obs_dt[0]

    else:
        # interpolate the data onto a common time step for each data type
	# elevation
        (mod_el_int, obs_el_int, step_int, start_int) = \
	    smooth(mod_el, mod_dt, obs_el, obs_dt)

	# speed
        (mod_sp_int, obs_sp_int, step_int, start_int) = \
            smooth(mod_spd, mod_dt, obs_spd, obs_dt)

	# direction
        (mod_dr_int, obs_dr_int, step_int, start_int) = \
            smooth(mod_dir, mod_dt, obs_dir, obs_dt)

	# u velocity
	(mod_u_int, obs_u_int, step_int, start_int) = \
	    smooth(mod_u, mod_dt, obs_u, obs_dt)

	# v velocity
	(mod_v_int, obs_v_int, step_int, start_int) = \
	    smooth(mod_v, mod_dt, obs_v, obs_dt)

	# velocity i.e. signed speed
	(mod_ve_int, obs_ve_int, step_int, start_int) = \
	    smooth(mod_spd * np.sign(mod_v), mod_dt, 
		   obs_spd * np.sign(obs_v), obs_dt)
    '''
    # separate into ebb and flow
    mod_dir_n = get_DirFromN(mod_u_int, mod_v_int)
    obs_dir_n = get_DirFromN(obs_u_int, mod_v_int)
    mod_signed_s, mod_PA = sign_speed(mod_u_int, mod_v_int, mod_sp_int,
				      mod_dr_int, 0)
    obs_signed_s, obs_PA = sign_speed(obs_u_int, obs_v_int, obs_sp_int,
				      obs_dr_int, 0)
    print mod_signed_s[:20], mod_PA[:20]
    print obs_signed_s[:20], obs_PA[:20]
    '''

    # remove directions where velocities are small
    MIN_VEL = 0.5
    for i in np.arange(obs_sp_int.size):
 	if (obs_sp_int[i] < MIN_VEL):
	    obs_dr_int[i] = np.nan
	if (mod_sp_int[i] < MIN_VEL):
	    mod_dr_int[i] = np.nan

    # get stats for each tidal variable
    elev_suite = tidalSuite(mod_el_int, obs_el_int, step_int, start_int,
			    type='elevation', plot=False)
    speed_suite = tidalSuite(mod_sp_int, obs_sp_int, step_int, start_int,
			    type='speed', plot=False)
    dir_suite = tidalSuite(mod_dr_int, obs_dr_int, step_int, start_int,
			    type='direction', plot=False)
    u_suite = tidalSuite(mod_u_int, obs_u_int, step_int, start_int,
			    type='u velocity', plot=False)
    v_suite = tidalSuite(mod_v_int, obs_v_int, step_int, start_int,
			    type='v velocity', plot=False)
    vel_suite = tidalSuite(mod_ve_int, obs_ve_int, step_int, start_int,
			    type='velocity', plot=False)
    #ebb_suite = tidalSuite(mod_ebb, obs_ebb, step_int, start_int,
	#		    type='ebb', plot=True)
    #flo_suite = tidalSuite(mod_flo, obs_flo, step_int, start_int,
	#		    type='flow', plot=True)
    # output statistics in useful format
    return (elev_suite, speed_suite, dir_suite, u_suite, v_suite, vel_suite)

def tidalSuite(model, observed, step, start, type, 
		  plot=False):
    '''
    Create stats classes for a given tidal variable.

    Accepts interpolated model and observed data, the timestep, and start
    time. Type is a string representing the type of data. If plot is set
    to true, a time plot and regression plot will be produced.
    
    Returns a dictionary containing all the stats.
    '''
    stats = TidalStats(model, observed, step, start, type=type)
    stats_suite = stats.getStats()
    stats_suite['r_squared'] = stats.linReg()['r_2']
    stats_suite['phase'] = stats.getPhase()

    if plot:
	stats.plotData()
	stats.plotRegression(stats.linReg())

    return stats_suite

def compareTG(data):
    '''
    Does a comprehensive comparison between tide gauge height data and
    modeled data, much like the above function.

    Input is a dictionary containing all necessary tide gauge and model data.
    Outputs a dictionary of useful statistics.
    '''
    # load data
    mod_elev = data['mod_timeseries']['elev']
    obs_elev = data['obs_timeseries']['elev']
    obs_datenums = data['obs_time']
    mod_datenums = data['mod_time']
    mod_harm = data['elev_mod_harmonics']

    # convert times and grab values
    obs_time, mod_time = [], []
    for i, v in enumerate(obs_datenums):
	obs_time.append(dn2dt(v))
    for j, w in enumerate(mod_datenums):
	mod_time.append(dn2dt(w))

    # check if they line up in the time domain
    if (mod_time[-1] < obs_time[0] or obs_time[-1] < mod_time[0]):

	# use ut_reconstr to create a new timeseries
	mod_elev_int = ut_reconstr(obs_datenums, mod_harm)[0]
	obs_elev_int = obs_elev
	step_int = obs_time[1] - obs_time[0]
	start_int = obs_time[0]

    else:

        # interpolate timeseries onto a common timestep
        (obs_elev_int, mod_elev_int, step_int, start_int) = \
            smooth(mod_elev, mod_time, obs_elev, obs_time)

    # get validation statistics
    stats = TidalStats(mod_elev_int, obs_elev_int, step_int, start_int,
		       debug=True, type='height')
    elev_suite = stats.getStats()
    elev_suite['r_squared'] = stats.linReg()['r_2']
    elev_suite['phase'] = stats.getPhase(debug=False)

    return elev_suite
