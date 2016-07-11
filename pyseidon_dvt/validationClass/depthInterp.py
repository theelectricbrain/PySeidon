#!/usr/bin/python2.7
# encoding: utf-8
import numpy as np
from scipy.interpolate import interp1d

'''
ASSUMPTIONS:
first dimension of matrices identify the timestep, second the depth
column vectors are organized from bottom to top
i.e. data[3] is the column at the third timestep
     data[3][10] is the tenth layer from the bottom at the third timestep
ADCP bins are the depths of each ADCP layer
ADCP bin width is constant
The top ADCP value of any column is no greater than 95% of the total depth
'''

ADCP_TOP_SURF = 0.95

def depthToSigma(obs_data, obs_depth, siglay, bins, debug=False, debug_plot=False):
    '''
    Performs linear interpolation on 3D ADCP data to change it into a sigma
    layer format, similar to an FVCOM run.

    Outputs a 2D numpy array representing the ADCP data in sigma layer
    format.
    '''
    if debug: print "depthToSigma..."
    sig_obs = np.zeros(obs_data.shape[0], siglay.size)

    if debug: print "...map old depths to between 0 and 1, then interpol..."
    # loop through columns/steps
    for i, column in enumerate(obs_data):
        # map old depths to between 0 and 1, make interpolation function
        col_nonan = column[np.where(~np.isnan(column))[0]]
        old_depths = bins[np.where(~np.isnan(column))[0]]
        mapped_depths = old_depths / obs_depth[i]
        f_obs = interp1d(col_nonan, mapped_depths)
        # perform interpolation
        sig_obs[i] = f_obs(siglay)

    if debug: print "...depthToSigma done."

    return sig_obs

def sigmaToDepth(mod_data, mod_depth, siglay, bins, debug=False, debug_plot=False):
    '''
    Performs linear interpolation on 3D FVCOM output to change it into the
    same format as ADCP output (i.e. constant depths, NaNs above surface)

    Outputs a 2D numpy array representing the FVCOM matrix in ADCP format.
    '''
    if debug: print "sigmaToDepth..."
    bin_mod = np.zeros(mod_data.shape[0], bins.size)
    bin_width = bins[1] - bins[0]

    if debug: print "...loop through columns/steps and interpol..."
    # loop through columns/steps
    for i, column in enumerate(mod_data):
        # create interpolation function
        f_mod = interp1d(column, siglay)
        depth = mod_depth[i]
        # loop through bins
        for j in np.arange(bins.size):
            # check if location is above ADCP_TOP_SURF
            loc = float(bins_width * j + bins[0]) / float(depth)
            if (loc <= ADCP_TOP_SURF):
                bin_mod[i][j] = f_mod(loc)
            else:
                bin_mod[i][j] = np.nan

    if debug: print "...sigmaToDepth done."

    return bin_mod

def depthFromSurf(mod_data, mod_depth, siglay, obs_data, obs_depth, bins, depth=5,
                  debug=False, debug_plot=False):
    '''
    Performs linear interpolation on 3D ocean data to obtain data at a
    specific distance from the surface.

    Inputs:
        - mod_data = 2D numpy array of FVCOM model data
        - mod_depth = 1D numpy array of model depths at each timestep
        - siglay = array containing values between 0 and 1 representing the
          respective percentage of depths for each sigma layer
        - obs_data = 2D numpy array of observed ADCP data
        - obs_depth = 1D numpy array of observed depths at each timestep
        - depth = number of metres from surface of output timeseries.
          Defaults to 5m

    Outputs:
        - (new_mod, new_obs) = timeseries representing model and observed data
                               at 'depth' metres from the surface.
    '''
    if debug: print "depthFromSurf..."
    new_mod = np.zeros(mod_data.shape[0])
    new_obs = np.zeros(obs_data.shape[0])
    depth = np.abs(depth)

    if debug: print "...loop through simulation columns and interpol at specified depth..."
    # loop over mod_data columns
    for i, step in enumerate(mod_data):      
        # create interpolation function
        #TR: quick fix
        try:
            f_mod = interp1d(np.abs(siglay), step) #, bounds_error=False)
            # find location of specified depth and perform interpolation
            location = mod_depth[i] - depth
            sig_loc = float(location) / float(mod_depth[i])
            new_mod[i] = f_mod(sig_loc)
        except ValueError:
            f_mod = interp1d(np.abs(siglay)[::-1], step[::-1], bounds_error=False)
            # find location of specified depth and perform interpolation
            location = mod_depth[i] - depth
            sig_loc = float(location) / float(mod_depth[i])
            new_mod[i] = f_mod(sig_loc)

    if debug: print "...loop through measurement columns and interpol at specified depth..."
    # loop over obs_data columns
    for ii, column in enumerate(obs_data):
        # create interpolation function
        col_nonan = column[np.where(~np.isnan(column))[0]]
        bin_nonan = bins[np.where(~np.isnan(column))[0]]

        if not col_nonan.shape[0]==0:
            # find location of specified depth and perform interpolation
            try:
                f_obs = interp1d(bin_nonan, col_nonan)
                location = obs_depth[ii] - depth
                new_obs[ii] = f_obs(location)
            except ValueError:
                f_obs = interp1d(bin_nonan[::-1], col_nonan[::-1], bounds_error=False)
                location = obs_depth[ii] - depth
                new_obs[ii] = f_obs(location)
        else:
            new_obs[ii] = np.nan

    if debug: print "...depthFromSurf done."

    return (new_mod, new_obs)

def depthFromBott(mod_data, mod_depth, siglay, obs_data, obs_depth, bins, depth=5,
                  debug=False, debug_plot=False):
    '''
    Performs linear interpolation on 3D ocean data to obtain data at a
    specific distance from the surface.

    Inputs:
        - mod_data = 2D numpy array of FVCOM model data
        - mod_depth = 1D numpy array of model depths at each timestep
        - siglay = array containing values between 0 and 1 representing the
          respective percentage of depths for each sigma layer
        - obs_data = 2D numpy array of observed ADCP data
        - obs_depth = 1D numpy array of observed depths at each timestep
        - depth = number of metres from surface of output timeseries.
          Defaults to 5m

    Outputs:
        - (new_mod, new_obs) = timeseries representing model and observed data
                               at 'depth' metres from the surface.
    '''
    if debug: print "depthFromBott..."
    new_mod = np.zeros(mod_data.shape[0])
    new_obs = np.zeros(obs_data.shape[0])
    depth = np.abs(depth)

    if debug: print "...loop through simulation columns and interpol at specified depth..."
    # loop over mod_data columns
    for i, step in enumerate(mod_data):
        # create interpolation function
        #TR: quick fix
        try:
            f_mod = interp1d(np.abs(siglay), step) #, bounds_error=False)
            # find location of specified depth and perform interpolation
            sig_loc = float(depth) / float(mod_depth[i])
            new_mod[i] = f_mod(sig_loc)
        except ValueError:
            f_mod = interp1d(np.abs(siglay)[::-1], step[::-1], bounds_error=False)
            # find location of specified depth and perform interpolation
            sig_loc = float(depth) / float(mod_depth[i])
            new_mod[i] = f_mod(sig_loc)

    if debug: print "...loop through measurement columns and interpol at specified depth..."
    # loop over obs_data columns
    for ii, column in enumerate(obs_data):
        # create interpolation function
        col_nonan = column[np.where(~np.isnan(column))[0]]
        bin_nonan = bins[np.where(~np.isnan(column))[0]]

        if not col_nonan.shape[0]==0:
            # find location of specified depth and perform interpolation
            try:
                f_obs = interp1d(bin_nonan, col_nonan) #, bounds_error=False)
                new_obs[ii] = f_obs(depth)
            except ValueError:
                f_obs = interp1d(bin_nonan[::-1], col_nonan[::-1], bounds_error=False)
                new_obs[ii] = f_obs(depth)
        else:
            new_obs[ii] = np.nan

    if debug: print "...depthFromBott done."

    return (new_mod, new_obs)