#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

import numpy as np
from scipy.stats import t, pearsonr
from datetime import datetime, timedelta
from scipy.interpolate import interp1d
from scipy.signal import correlate
import time
import pandas as pd
from sys import exit

# local imports
from pyseidon_dvt.utilities.BP_tools import principal_axis

# Custom error
from pyseidon_dvt.utilities.pyseidon_error import PyseidonError

class TidalStats:
    """
    An object representing a set of statistics on tidal heights used
    to determine the skill of a model in comparison to observed data.
    Standards are from NOAA's Standard Suite of Statistics.

    Instantiated with two arrays containing predicted and observed
    data which have already been interpolated so they line up, the
    time step between points, and the start time of the data.

    To remove NaNs in observed data, linear interpolation is performed to
    fill gaps. Additionally, NaNs are trimmed from the start and end.

    Functions are used to calculate statistics and to output
    visualizations and tables.
    """
    def __init__(self, gear, model_data, observed_data, time_step, start_time,
                 model_u = [], observed_u = [], model_v= [], observed_v = [],
                 model_time = [], observed_time = [], phase_shift=False,
                 kind='', debug=False, debug_plot=False):
        if debug: print "TidalStats initialisation..."
        self._debug = debug
        self._debug_plot = debug_plot
        self.gear = gear # Type of measurement gear (drifter, adcp,...), str
        self.model = np.asarray(model_data)
        self.model = self.model.astype(np.float64)
        self.observed = np.asarray(observed_data)
        self.observed = self.observed.astype(np.float64)
        self.model_u = np.asarray(model_u)
        self.model_u = self.model_u.astype(np.float64)
        self.observed_u = np.asarray(observed_u)
        self.observed_u = self.observed_u.astype(np.float64)
        self.model_v = np.asarray(model_v)
        self.model_v = self.model_v.astype(np.float64)
        self.observed_v = np.asarray(observed_v)
        self.observed_v = self.observed_v.astype(np.float64)
        self.kind = kind
        self.step = time_step
        self.length = model_data.size

        # TR: pass this step if dealing with Drifter's data
        if not self.gear == 'Drifter':
            try:
                # TR: fix for interpolation pb when 0 index or -1 index array values = nan
                if debug: print "...trim nans at start and end of data.."
                start_index, end_index = 0, -1
                while np.isnan(self.observed[start_index]) or np.isnan(self.model[start_index]):
                    start_index += 1
                while np.isnan(self.observed[end_index]) or np.isnan(self.model[end_index]):
                    end_index -= 1
            except IndexError:  # due too nans everywhere, itself due to no obs data at requested depth
                raise  PyseidonError("-No matching measurement-")

            # Correction for bound index call
            if end_index == -1:
                end_index = None
            else:
                end_index += 1
            if debug: print "Start index: ", start_index
            if debug: print "End index: ", end_index

            m = self.model[start_index:end_index]
            o = self.observed[start_index:end_index]

            setattr(self, 'model', m)
            setattr(self, 'observed', o)
        
            # set up array of datetimes corresponding to the data (and timestamps)
            self.times = start_time + np.arange(self.model.size) * time_step
            timestamps = np.zeros(len(self.times))
            for j, jj in enumerate(self.times):
                timestamps[j] = time.mktime(jj.timetuple())

            if debug: print "...uses linear interpolation to eliminate any NaNs in the data..."

            if True in np.isnan(self.observed):
                time_nonan = timestamps[np.where(~np.isnan(self.observed))[0]]
                obs_nonan = self.observed[np.where(~np.isnan(self.observed))[0]]                
                func = interp1d(time_nonan, obs_nonan)
                self.observed = func(timestamps)
            if True in np.isnan(self.model):
                time_nonan = timestamps[np.where(~np.isnan(self.model))[0]]
                mod_nonan = self.model[np.where(~np.isnan(self.model))[0]]                
                func = interp1d(time_nonan, mod_nonan)
                self.model = func(timestamps)

        # TR: pass this step if dealing with Drifter's data
        else:
            self.times = start_time + np.arange(self.model.size) * time_step # needed for plots, en seconds
            #TR: those are not the real times though

        # Applying phase shift correction if needed
        if phase_shift:
            if debug: print "...applying phase shift..."
            try:  # Fix for Drifter's data
                step_sec = time_step.seconds
            except AttributeError:
                step_sec = time_step * 24.0 * 60.0 * 60.0  # converts matlabtime to seconds
            phase = self.getPhase()
            phaseIndex = int(phase * 60.0 / step_sec)
            if debug: print "Phase index = "+str(phaseIndex)
            # create shifted data)
            # if (phaseIndex < 0):
            #     # left shift
            #     iM = np.s_[-phaseIndex:]
            #     iO = np.s_[:self.length + phaseIndex]
            # elif (phaseIndex > 0):
            #     # right shift
            #     iM = np.s_[:self.length - phaseIndex]
            #     iO = np.s_[phaseIndex:]
            # else:  # if (phaseIndex == 0):
            #     # no shift
            #     iM = np.s_[:]
            #     iO = np.s_[:]
            # self.model = self.model[iM]
            # self.observed = self.observed[iO]
            # if not model_u == []:
            #     self.model_u = self.model_u[iM]
            # if not model_v == []:
            #     self.model_v = self.model_v[iM]
            # if not observed_u == []:
            #     self.observed_u = self.observed_v[iO]
            # if not observed_v == []:
            #     self.observed_v = self.observed_v[iO]
            self.model = np.roll(self.model, phaseIndex)
            if not model_u == []:
                self.model_u = np.roll(self.model_u, phaseIndex)
            if not model_v == []:
                self.model_v = np.roll(self.model_v, phaseIndex)
        if debug: print "...TidalStats initialisation done."

        # Error attributes
        if self.kind in ['cubic speed', 'velocity', 'direction']:
            # TR: pass this step if dealing with Drifter's data
            if not self.gear == 'Drifter':
                # interpolate cubic speed, u and v on same time steps
                model_timestamps = np.zeros(len(model_time))
                for j, jj in enumerate(model_time):
                    model_timestamps[j] = time.mktime(jj.timetuple())
                obs_timestamps = np.zeros(len(observed_time))
                for j, jj in enumerate(observed_time):
                    obs_timestamps[j] = time.mktime(jj.timetuple())
                func_u = interp1d(model_timestamps, model_u)
                self.model_u = func_u(timestamps)
                func_v = interp1d(model_timestamps, model_v)
                self.model_v = func_v(timestamps)
                func_u = interp1d(obs_timestamps, observed_u)
                self.observed_u = func_u(timestamps)
                func_v = interp1d(obs_timestamps, observed_v)
                self.observed_v = func_v(timestamps)

            # if self.kind == 'cubic speed':
            #     # R.Karsten formula
            #     self.error = ((self.model_u**2.0 + self.model_v**2.0)**(3.0/2.0)) - \
            #                  ((self.observed_u**2.0 + self.observed_v**2.0)**(3.0/2.0))
            # else:
            #     self.error = self.observed - self.model
            self.error = self.observed - self.model
        elif self.kind in ['speed', 'elevation', 'u velocity', 'v velocity', 'Phase']:
            self.error = self.observed - self.model
        else:
            print "---Data kind not supported---"
            exit()

        # if debug: print "...establish limits as defined by NOAA standard..."
        # if self.kind == 'velocity':
        #     self.ERROR_BOUND = 0.26
        # elif (self.kind == 'speed' or self.kind == 'velocity'):
        #     self.ERROR_BOUND = 0.26
        # elif (self.kind == 'elevation'):
        #     self.ERROR_BOUND = 0.15
        # elif (self.kind == 'direction'):
        #     self.ERROR_BOUND = 22.5
        # elif (self.kind == 'u velocity' or self.kind == 'v velocity'):
        #     self.ERROR_BOUND = 0.35
        # elif self.kind == 'cubic speed':
        #     self.ERROR_BOUND = 0.26**3.0
        # else:
        #     self.ERROR_BOUND = 0.5

        # instead of using the NOAA errors, use 10% of the data range
        obs_range = 0.1 * (np.nanmax(self.observed) - np.nanmin(self.observed))
        mod_range = 0.1 * (np.nanmax(self.model) - np.nanmin(self.model))
        self.ERROR_BOUND = (obs_range + mod_range) / 2.

        return

    def getRMSE(self, debug=False):
        '''
        Returns the root mean squared error of the data.
        '''
        if debug or self._debug: print "...getRMSE..."
        if self.kind == 'velocity':
            # Special definition of rmse - R.Karsten
            rmse = np.sqrt(np.nanmean((self.model_u - self.observed_u)**2.0 + (self.model_v - self.observed_v)**2.0))
        else:
            rmse = np.sqrt(np.nanmean(self.error**2))
        return rmse


    def getNRMSE(self, debug=False):
        """
        Returns the normalized root mean squared error between the model and
        observed data in %.
        """
        if debug or self._debug: print "...getNRMSE..."
        if self.kind == 'velocity':
            # Special definition of rmse - R.Karsten
            rmse0 = np.sqrt(np.nanmean((self.observed_u)**2.0 + (self.observed_v)**2.0))
        else:
            rmse0 = np.sqrt(np.mean(self.observed**2.0))
        # return 100. * self.getRMSE() / (max(self.observed) - min(self.observed))
        return 100. * self.getRMSE() / rmse0

    def getSD(self, debug=False):
        '''
        Returns the standard deviation of the error.
        '''
        if debug or self._debug: print "...getSD..."
        return np.sqrt(np.nanmean(abs(self.error - np.nanmean(self.error)**2)))

    def getBias(self, debug=False):
        """
        Returns the bias of the model, a measure of over/under-estimation.
        """
        if debug or self._debug: print "...getBias..."
        return np.nanmean(self.error)

    def getSI(self, debug=False):
        """
        Returns the scatter index of the model, a weighted measure of data
        scattering.
        """
        if debug or self._debug: print "...getSI..."
        return self.getRMSE() / np.nanmean(self.observed)

    def getPBIAS(self, debug=False):
        """
        Returns the percent bias between the model and the observed data.

        References:
          Yapo P. O., Gupta H. V., Sorooshian S., 1996.
          Automatic calibration of conceptual rainfall-runoff models: sensitivity to calibration data.
          Journal of Hydrology. v181 i1-4. 23-48

          Sorooshian, S., Q. Duan, and V. K. Gupta. 1993.
          Calibration of rainfall-runoff models: Application of global optimization
          to the Sacramento Soil Moisture Accounting Model.
          Water Resources Research, 29 (4), 1185-1194, doi:10.1029/92WR02617.
        """
        if debug or self._debug: print "...getPBIAS..."

        # if self.kind in ['elevation', 'direction', 'u velocity', 'v velocity', 'velocity']:
        #     norm_error = self.error / self.observed
        #     pbias = 100. * np.sum(norm_error) / norm_error.size
        # else:
        #     norm_error = self.model - self.observed
        #     pbias = 100. * (np.sum(norm_error) / np.sum(self.observed))
        #     # TR: this formula may lead to overly large values when used with sinusoidal signals

        pbias = 100. * (np.nansum(self.error) / np.nansum(self.observed))

        return pbias


    def getNSE(self, debug=False):
        """
        Returns the Nash-Sutcliffe Efficiency coefficient of the model vs.
        the observed data. Identifies if the model is better for
        approximation than the mean of the observed data.
        """
        SSE_mod = np.nansum((self.observed - self.model)**2)
        SSE_mean = np.nansum((self.observed - np.nanmean(self.observed))**2)
        return 1 - SSE_mod / SSE_mean

    def getCORR(self, debug=False):
        """
        Returns the Pearson correlation coefficient for the model vs.
        the observed data, a number between -1 and 1. -1 implies perfect
        negative correlation, 1 implies perfect correlation.
        """
        return pearsonr(self.observed, self.model)[0]

    def statsForDirection(self, debug=False):
        """
        Special stats for direction

        Outputs:
          - err = absolute error
          - nerr = absolute error divided by standard deviation in %
        """
        if debug: print "Computing special stats for direction..."
        pr_axis_mod, pr_ax_var_mod = principal_axis(self.model_u, self.model_v)
        pr_axis_obs, pr_ax_var_obs = principal_axis(self.observed_u, self.observed_v)

        # Defines intervals
        dir_all_mod = self.model[:]
        dir_all_obs = self.observed[:]
        ind_mod = np.where(dir_all_mod<0)
        ind_obs = np.where(dir_all_obs<0)
        dir_all_mod[ind_mod] = dir_all_mod[ind_mod] + 360
        dir_all_obs[ind_obs] = dir_all_obs[ind_obs] + 360

        # sign speed - eliminating wrap-around
        dir_PA_mod = dir_all_mod - pr_axis_mod
        dir_PA_mod[dir_PA_mod < -90] += 360
        dir_PA_mod[dir_PA_mod > 270] -= 360
        dir_PA_obs = dir_all_obs - pr_axis_obs
        dir_PA_obs[dir_PA_obs < -90] += 360
        dir_PA_obs[dir_PA_obs > 270] -= 360

        #general direction of flood passed as input argument
        floodIndex_mod = np.where((dir_PA_mod >= -90) & (dir_PA_mod<90))[0]
        ebbIndex_mod = np.arange(dir_PA_mod.shape[0])
        ebbIndex_mod = np.delete(ebbIndex_mod, floodIndex_mod[:])
        floodIndex_obs = np.where((dir_PA_obs >= -90) & (dir_PA_obs<90))[0]
        ebbIndex_obs = np.arange(dir_PA_obs.shape[0])
        ebbIndex_obs = np.delete(ebbIndex_obs, floodIndex_obs[:])

        # principal axis for ebb and flood
        np.delete(floodIndex_mod, np.where(floodIndex_mod >= self.model_u.shape[0]))
        np.delete(ebbIndex_mod, np.where(ebbIndex_mod >= self.model_u.shape[0]))
        np.delete(floodIndex_obs, np.where(floodIndex_obs >= self.observed_u.shape[0]))
        np.delete(ebbIndex_obs, np.where(ebbIndex_obs >= self.observed_u.shape[0]))

        pr_axis_mod_flood, pr_ax_var_mod_flood = principal_axis(self.model_u[floodIndex_mod],
                                                                self.model_v[floodIndex_mod])
        pr_axis_mod_ebb, pr_ax_var_mod_ebb = principal_axis(self.model_u[ebbIndex_mod],
                                                            self.model_v[ebbIndex_mod])
        pr_axis_obs_flood, pr_ax_var_obs_flood = principal_axis(self.observed_u[floodIndex_obs],
                                                                self.observed_v[floodIndex_obs])
        pr_axis_obs_ebb, pr_ax_var_obs_ebb = principal_axis(self.observed_u[ebbIndex_obs],
                                                            self.observed_v[ebbIndex_obs])
        err_flood = np.abs(pr_axis_mod_flood - pr_axis_obs_flood)
        err_ebb = np.abs(pr_axis_mod_ebb - pr_axis_obs_ebb)

        if debug: print "...ebb direction error: " + str(err_ebb) + " degrees..."
        if debug: print "...flood direction error: " + str(err_flood) + " degrees..."

        err = 0.5 * (err_flood + err_ebb)

        #std_flood = np.asarray(dir_all_obs[floodIndex_obs]).std()
        #std_ebb = np.asarray(dir_all_obs[ebbIndex_obs]).std()
        #std_flood = dir_all_obs[floodIndex_obs].max() - dir_all_obs[floodIndex_obs].min()
        #std_ebb = dir_all_obs[ebbIndex_obs].max() - dir_all_obs[ebbIndex_obs].min()
        #nerr = 0.5 * (np.abs(err_flood / std_flood) + np.abs(err_ebb / std_ebb))
        nerr = err / 90.0

        return err, nerr * 100.0


    def getCF(self, debug=False):
        '''
        Returns the central frequency of the data, i.e. the fraction of
        errors that lie within the defined limit.
        '''
        central_err = [i for i in self.error if abs(i) < self.ERROR_BOUND]
        central_num = len(central_err)
        if debug or self._debug: print "...getCF..."
        return (float(central_num) / float(self.length)) * 100

    def getPOF(self, debug=False):
        '''
        Returns the positive outlier frequency of the data, i.e. the
        fraction of errors that lie above the defined limit.
        '''
        upper_err = [i for i in self.error if i > 2 * self.ERROR_BOUND]
        upper_num = len(upper_err)
        if debug or self._debug: print "...getPOF..."
        return (float(upper_num) / float(self.length)) * 100

    def getNOF(self, debug=False):
        '''
        Returns the negative outlier frequency of the data, i.e. the
        fraction of errors that lie below the defined limit.
        '''
        lower_err = [i for i in self.error if i < -2 * self.ERROR_BOUND]
        lower_num = len(lower_err)
        if debug or self._debug: print "...getNOF..."
        return (float(lower_num) / float(self.length)) * 100

    def getMDPO(self, debug=False):
        '''
        Returns the maximum duration of positive outliers, i.e. the
        longest amount of time across the data where the model data
        exceeds the observed data by a specified limit.

        Takes one parameter: the number of minutes between consecutive
        data points.
        '''
        try: #Fix for Drifter's data
            timestep = self.step.seconds / 60
        except AttributeError:
            timestep = self.step * 24.0 * 60.0 # converts matlabtime (in days) to minutes

        max_duration = 0
        current_duration = 0
        for i in np.arange(self.error.size):
            if (self.error[i] > self.ERROR_BOUND):
                current_duration += timestep
            else:
                if (current_duration > max_duration):
                    max_duration = current_duration
                current_duration = 0
        if debug or self._debug: print "...getMDPO..."
        return max(max_duration, current_duration)

    def getMDNO(self, debug=False):
        '''
        Returns the maximum duration of negative outliers, i.e. the
        longest amount of time across the data where the observed
        data exceeds the model data by a specified limit.

        Takes one parameter: the number of minutes between consecutive
        data points.
        '''
        try: #Fix for Drifter's data
            timestep = self.step.seconds / 60
        except AttributeError:
            timestep = self.step * 24.0 * 60.0 # converts matlabtime (in days) to minutes

        max_duration = 0
        current_duration = 0
        for i in np.arange(self.error.size):
            if (self.error[i] < -self.ERROR_BOUND):
                current_duration += timestep
            else:
                if (current_duration > max_duration):
                    max_duration = current_duration
                current_duration = 0
        if debug or self._debug: print "...getMDNO..."
        return max(max_duration, current_duration)

    def getMSE(self, debug=False):
        """
        Returns the mean square error (float)
        """
        mse = np.nanmean(self.error**2.0)
        if debug or self._debug: print "...getMSE..."
        return mse

    def getNMSE(self, debug=False):
        """
        Returns the normalized mean square error in % (float)
        """
        mse0 = np.nanmean(self.observed**2.0)
        nmse = 100. * self.getMSE() / mse0
        if debug or self._debug: print "...getNMSE..."
        return nmse


    def getWillmott(self, debug=False):
        '''
        Returns the Willmott skill statistic.
        '''

        # start by calculating MSE
        MSE = self.getMSE()

        # now calculate the rest of it
        obs_mean = np.nanmean(self.observed)
        skill = 1 - MSE / np.nanmean((abs(self.model - obs_mean) +
                                   abs(self.observed - obs_mean))**2)
        if debug or self._debug: print "...getWillmott..."
        return skill

    def getPhase(self, max_phase=timedelta(hours=3), debug=False):
        '''
        Attempts to find the phase shift between the model data and the
        observed data.

        Iteratively tests different phase shifts, and calculates the RMSE
        for each one. The shift with the smallest RMSE is returned.

        Argument max_phase is the span of time across which the phase shifts
        will be tested. If debug is set to True, a plot of the RMSE for each
        phase shift will be shown.
        '''
        if debug or self._debug: print "getPhase..."
        # grab the length of the timesteps in seconds
        max_phase_sec = max_phase.seconds
        try:  # Fix for Drifter's data
            step_sec = self.step.seconds
        except AttributeError:
            step_sec = self.step * 24.0 * 60.0 * 60.0 # converts matlabtime to seconds

        num_steps = max_phase_sec / step_sec

        if debug or self._debug: print "...iterate through the phase shifts and check RMSE..."
        errors = []
        phases = np.arange(-num_steps, num_steps).astype(int)
        for i in phases:
            # create shifted data
            shift_mod = np.roll(self.model, i)
            # if (i < 0):
            #     # left shift
            #     #shift_mod = self.model[-i:]
            #     shift_obs = self.observed[:self.length + i]
            # if (i > 0):
            #     # right shift
            #     shift_mod = self.model[:self.length - i]
            #     shift_obs = self.observed[i:]
            # if (i == 0):
            #     # no shift
            #     shift_mod = self.model
            #     shift_obs = self.observed

            start = self.times[abs(i)]
            step = self.times[1] - self.times[0]

            # create TidalStats class for shifted data and get the RMSE
            #stats = TidalStats(self.gear, shift_mod, shift_obs, step, start, kind='Phase')
            stats = TidalStats(self.gear, shift_mod, self.observed, step, start, kind='Phase')
            nrms_error = stats.getNRMSE()
            errors.append(nrms_error)

        if debug or self._debug: print "...find the minimum rmse, and thus the minimum phase..."
        min_index = errors.index(min(errors))
        best_phase = phases[min_index]
        phase_minutes = best_phase * step_sec / 60

        return phase_minutes

    def altPhase(self, debug=False):
        """
        Alternate version of lag detection using scipy's cross correlation function.
        """
        if debug or self._debug: print "altPhase..."
        # normalize arrays
        mod = self.model
        mod -= self.model.mean()
        mod /= mod.std()
        obs = self.observed
        obs -= self.observed.mean()
        obs /= obs.std()

        if debug or self._debug: print "...get cross correlation and find number of timesteps of shift..."
        xcorr = correlate(mod, obs)
        samples = np.arange(1 - self.length, self.length)
        time_shift = samples[xcorr.argmax()]

        # find number of minutes in time shift
        try: #Fix for Drifter's data
            step_sec = self.step.seconds
        except AttributeError:
            step_sec = self.step * 24.0 * 60.0 * 60.0 # converts matlabtime (in days) to seconds
        lag = time_shift * step_sec / 60

        if debug or self._debug: print "...altPhase done."

        return lag

    def getStats(self, phase_shift=False, debug=False):
        """
        Returns each of the statistics in a dictionary.
        """

        stats = {}
        stats['gear'] = self.gear
        stats['RMSE'] = self.getRMSE()
        stats['CF'] = self.getCF()
        stats['SD'] = self.getSD()
        stats['POF'] = self.getPOF()
        stats['NOF'] = self.getNOF()
        stats['MDPO'] = self.getMDPO()
        stats['MDNO'] = self.getMDNO()
        stats['skill'] = self.getWillmott()
        if not phase_shift:
            stats['phase'] = self.getPhase(debug=debug)
        else:
            stats['phase'] = 0.0
        #stats['phase'] = self.getPhase()
        stats['CORR'] = self.getCORR()
        stats['NRMSE'] = self.getNRMSE()
        stats['NSE'] = self.getNSE()
        stats['bias'] = self.getBias()
        stats['SI'] = self.getSI()
        stats['pbias'] = self.getPBIAS()
        stats['MSE'] = self.getMSE()
        stats['NMSE'] = self.getNMSE()

        if debug or self._debug: print "...getStats..."

        return stats

    def linReg(self, alpha=0.05, debug=False):
        '''
        Does linear regression on the model data vs. recorded data.

        Gives a 100(1-alpha)% confidence interval for the slope
        '''
        if debug or self._debug: print "linReg..."
        # set stuff up to make the code cleaner
        obs = self.observed
        mod = self.model
        obs_mean = np.mean(obs)
        mod_mean = np.mean(mod)
        n = mod.size
        df = n - 2

        # calculate square sums
        SSxx = np.sum(mod**2) - np.sum(mod)**2 / n
        SSyy = np.sum(obs**2) - np.sum(obs)**2 / n
        SSxy = np.sum(mod * obs) - np.sum(mod) * np.sum(obs) / n
        SSE = SSyy - SSxy**2 / SSxx
        MSE = SSE / df

        # estimate parameters
        slope = SSxy / SSxx
        intercept = obs_mean - slope * mod_mean
        sd_slope = np.sqrt(MSE / SSxx)
        r_squared = 1 - SSE / SSyy

        # calculate 100(1 - alpha)% CI for slope
        width = t.isf(0.5 * alpha, df) * sd_slope
        lower_bound = slope - width
        upper_bound = slope + width
        slope_CI = (lower_bound, upper_bound)

        # calculate 100(1 - alpha)% CI for intercept
        lower_intercept = obs_mean - lower_bound * mod_mean
        upper_intercept = obs_mean - upper_bound * mod_mean
        intercept_CI = (lower_intercept, upper_intercept)

        # estimate 100(1 - alpha)% CI for predictands
        predictands = slope * mod + intercept
        sd_resid = np.std(obs - predictands)
        y_CI_width = t.isf(0.5 * alpha, df) * sd_resid * \
            np.sqrt(1 - 1 / n)

        # return data in a dictionary
        data = {}
        data['slope'] = slope
        data['intercept'] = intercept
        data['r_2'] = r_squared
        data['slope_CI'] = slope_CI
        data['intercept_CI'] = intercept_CI
        data['pred_CI_width'] = y_CI_width
        data['conf_level'] = 100 * (1 - alpha)

        if debug or self._debug: print "...linReg done."

        return data

    def crossVal(self, alpha=0.05, debug=False):
        '''
        Performs leave-one-out cross validation on the linear regression.

        i.e. removes one datum from the set, redoes linreg on the training
        set, and uses the results to attempt to predict the missing datum.
        '''
        if debug or self._debug: print "crossVal..."
        cross_error = np.zeros(self.model.size)
        cross_pred = np.zeros(self.model.size)
        model_orig = self.model
        obs_orig = self.observed
        time_orig = self.time

        if debug or self._debug: print "...loop through each element, remove it..."
        for i in np.arange(self.model.size):
            train_mod = np.delete(model_orig, i)
            train_obs = np.delete(obs_orig, i)
            train_time = np.delete(time_orig, i)
            train_stats = TidalStats(train_mod, train_obs, train_time)

            # redo the linear regression and get parameters
            param = train_stats.linReg(alpha)
            slope = param['slope']
            intercept = param['intercept']

            # predict the missing observed value and calculate error
            pred_obs = slope * model_orig[i] + intercept
            cross_pred[i] = pred_obs
            cross_error[i] = abs(pred_obs - obs_orig[i])

        # calculate PRESS and PRRMSE statistics for predicted data
        if debug or self._debug: print "...predicted residual sum of squares and predicted RMSE..."
        PRESS = np.sum(cross_error**2)
        PRRMSE = np.sqrt(PRESS) / self.model.size

        # return data in a dictionary
        data = {}
        data['PRESS'] = PRESS
        data['PRRMSE'] = PRRMSE
        data['cross_pred'] = cross_pred

        if debug or self._debug: print "...crossVal done."

        return data

    def save_data(self, path=''):
            df = pd.DataFrame(data={'time': self.times.ravel(),
                                    'observed':self.observed.ravel(),
                                    'modeled':self.model.ravel() })
            df.to_csv(path+str(self.kind)+'.csv')
