#!/usr/bin/python2.7
# encoding: utf-8
import numpy as np
from scipy.stats import t, pearsonr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy.interpolate import interp1d
from scipy.signal import correlate
import time
import seaborn
import pandas as pd
from sys import exit

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
                 model_time = [], observed_time = [],
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

        # TR: pass this step if dealing with Drifter's data
        if not self.gear == 'Drifter':
            # TR: fix for interpolation pb when 0 index or -1 index array values = nan
            if debug: print "...trim nans at start and end of data.."
            start_index, end_index = 0, -1
            while (np.isnan(self.observed[start_index]) or np.isnan(self.model[start_index])):
                start_index += 1
            while (np.isnan(self.observed[end_index]) or np.isnan(self.model[end_index])):
                end_index -= 1

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
            self.step = time_step
            timestamps = np.zeros(len(self.times))
            for j, jj in enumerate(self.times):
                timestamps[j] = time.mktime(jj.timetuple())

            if debug: print "...uses linear interpolation to eliminate any NaNs in the data..."

            if (True in np.isnan(self.observed)):
                obs_nonan = self.observed[np.where(~np.isnan(self.observed))[0]]
                time_nonan = timestamps[np.where(~np.isnan(self.observed))[0]]
                func = interp1d(time_nonan, obs_nonan)
                self.observed = func(timestamps)
            if (True in np.isnan(self.model)):
                mod_nonan = self.model[np.where(~np.isnan(self.model))[0]]
                time_nonan = timestamps[np.where(~np.isnan(self.model))[0]]
                func = interp1d(time_nonan, mod_nonan)
                self.model = func(timestamps)

        # TR: pass this step if dealing with Drifter's data
        else:
            self.step = time_step # needed for getMDPO, getMDNO, getPhase & altPhase
            self.times = start_time + np.arange(self.model.size) * time_step # needed for plots, en seconds
            #TR: those are not the real times though

        # Error attributes
        if self.kind in ['power density', 'velocity']:
            # TR: pass this step if dealing with Drifter's data
            if not self.gear == 'Drifter':
            # interpolate power density, u and v on same time steps
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

            if self.kind == 'power density':
                # R.Karsten formula
                self.error = ((self.model_u**2.0 + self.model_v**2.0)**(3.0/2.0)) - \
                             ((self.observed_u**2.0 + self.observed_v**2.0)**(3.0/2.0))
                self.error = 0.5 * 1025.0 * self.error
            else:
                self.error = self.observed - self.model
        elif self.kind in ['speed', 'elevation', 'direction', 'u velocity', 'v velocity', 'Phase']:
            self.error = self.observed - self.model
        else:
            print "---Data kind not supported---"
            exit()
        self.length = self.error.size

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
        # elif self.kind == 'power density':
        #     self.ERROR_BOUND = 0.5*1025.0*0.26**3.0
        # else:
        #     self.ERROR_BOUND = 0.5

        # instead of using the NOAA errors, use 10% of the data range
        obs_range = 0.1 * (np.max(self.observed) - np.min(self.observed))
        mod_range = 0.1 * (np.max(self.model) - np.min(self.model))
        self.ERROR_BOUND = (obs_range + mod_range) / 2.

        if debug: print "...TidalStats initialisation done."

    def getOverUnder(self, debug=False):
        """
        Determines if model over or under estimate the reference

        Returns:
          - ovORun = character, '+' if overestimation and '-' if underestimation
        """
        if debug or self._debug: print "...getOverUnder..."
        if np.var(self.model) > np.var(self.observed[:]):
            ovORun = '+'
        else:
            ovORun = '-'
        return ovORun

    def getRMSE(self, debug=False):
        '''
        Returns the root mean squared error of the data.
        '''
        if debug or self._debug: print "...getRMSE..."
        if self.kind == 'velocity':
            # Special definition of rmse - R.Karsten
            rmse = np.sqrt(np.mean((self.model_u - self.observed_u)**2.0 + (self.model_v - self.observed_v)**2.0))
        elif self.kind == 'power density':
            rmse = np.sqrt(np.mean(self.error))
        else:
            rmse = np.sqrt(np.mean(self.error**2))
        return rmse

    def getSD(self, debug=False):
        '''
        Returns the standard deviation of the error.
        '''
        if debug or self._debug: print "...getSD..."
        return np.sqrt(np.mean(abs(self.error - np.mean(self.error)**2)))

    def getBias(self, debug=False):
        """
        Returns the bias of the model, a measure of over/under-estimation.
        """
        if debug or self._debug: print "...getBias..."
        return np.mean(self.error)

    def getSI(self, debug=False):
        """
        Returns the scatter index of the model, a weighted measure of data
        scattering.
        """
        if debug or self._debug: print "...getSI..."
        return self.getRMSE() / np.mean(self.observed)


    def getNRMSE(self, debug=False):
        """
        Returns the normalized root mean squared error between the model and
        observed data in %.
        """
        if debug or self._debug: print "...getNRMSE..."
        return 100. * self.getRMSE() / (max(self.observed) - min(self.observed))

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

        if self.kind in ['elevation', 'direction', 'u velocity', 'v velocity', 'velocity']:
            norm_error = self.error / self.observed
            pbias = 100. * np.sum(norm_error) / norm_error.size
        else:
            norm_error = self.model - self.observed
            pbias = 100. * (np.sum(norm_error) / np.sum(self.observed))
            # TR: this formula leads to overly large values when used with sinusoidal signals
        return pbias


    def getNSE(self, debug=False):
        """
        Returns the Nash-Sutcliffe Efficiency coefficient of the model vs.
        the observed data. Identifies if the model is better for
        approximation than the mean of the observed data.
        """
        SSE_mod = np.sum((self.observed - self.model)**2)
        SSE_mean = np.sum((self.observed - np.mean(self.observed))**2)
        return 1 - SSE_mod / SSE_mean

    def getCORR(self, debug=False):
        """
        Returns the Pearson correlation coefficient for the model vs.
        the observed data, a number between -1 and 1. -1 implies perfect
        negative correlation, 1 implies perfect correlation.
        """
        return pearsonr(self.observed, self.model)[0]

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
            timestep = self.step / 60            

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
            timestep = self.step / 60 

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

    def getWillmott(self, debug=False):
        '''
        Returns the Willmott skill statistic.
        '''

        # start by calculating MSE
        MSE = np.mean(self.error**2)

        # now calculate the rest of it
        obs_mean = np.mean(self.observed)
        skill = 1 - MSE / np.mean((abs(self.model - obs_mean) +
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
        try: #Fix for Drifter's data
            step_sec = self.step.seconds
        except AttributeError:
            step_sec = self.step / 60

        num_steps = max_phase_sec / step_sec

        if debug or self._debug: print "...iterate through the phase shifts and check RMSE..."
        errors = []
        phases = np.arange(-num_steps, num_steps)
        for i in phases:
            # create shifted data
            if (i < 0):
                # left shift
                shift_mod = self.model[-i:]
                shift_obs = self.observed[:self.length + i]
            if (i > 0):
                # right shift
                shift_mod = self.model[:self.length - i]
                shift_obs = self.observed[i:]
            if (i == 0):
                # no shift
                shift_mod = self.model
                shift_obs = self.observed

            start = self.times[abs(i)]
            step = self.times[1] - self.times[0]

        # create TidalStats class for shifted data and get the RMSE
        stats = TidalStats(self.gear, shift_mod, shift_obs, step, start, kind='Phase')

        rms_error = stats.getRMSE()
        errors.append(rms_error)

        if debug or self._debug: print "...find the minimum rmse, and thus the minimum phase..."
        min_index = errors.index(min(errors))
        best_phase = phases[min_index]
        phase_minutes = best_phase * step_sec / 60

        return phase_minutes

    def altPhase(self, debug=False):
        """
        Alternate version of lag detection using scipy's cross correlation
        function.
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
            step_sec = self.step
        lag = time_shift * step_sec / 60

        if debug or self._debug: print "...altPhase done."

        return lag

    def getStats(self, debug=False):
        """
        Returns each of the statistics in a dictionary.
        """

        stats = {}
        stats['ovORun'] = self.getOverUnder()
        stats['RMSE'] = self.getRMSE()
        stats['CF'] = self.getCF()
        stats['SD'] = self.getSD()
        stats['POF'] = self.getPOF()
        stats['NOF'] = self.getNOF()
        stats['MDPO'] = self.getMDPO()
        stats['MDNO'] = self.getMDNO()
        stats['skill'] = self.getWillmott()
        try: #Fix for Drifter's data
            stats['phase'] = self.getPhase(debug=debug)
        except:
            stats['phase'] = 0.0
        stats['CORR'] = self.getCORR()
        stats['NRMSE'] = self.getNRMSE()
        stats['NSE'] = self.getNSE()
        stats['bias'] = self.getBias()
        stats['SI'] = self.getSI()
        stats['pbias'] = self.getPBIAS()
        stats['phase'] = self.getPhase(debug=debug)

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

    def plotRegression(self, lr, save=False, out_f='', debug=False):
        """
        Plots a visualization of the output from linear regression,
        including confidence intervals for predictands and slope.

        If save is set to True, exports the plot as an image file to out_f.
        """
        df = pd.DataFrame(data={'model': self.model.flatten(),
                                'observed':self.observed.flatten()})
        #define figure frame
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        ax = fig.add_subplot(111)

        ax.scatter(self.model, self.observed, c='b', marker='+', alpha=0.5)

        ## plot regression line
        mod_max = np.amax(self.model)
        mod_min = np.amin(self.model)
        upper_intercept = lr['intercept'] + lr['pred_CI_width']
        lower_intercept = lr['intercept'] - lr['pred_CI_width']
        ax.plot([mod_min, mod_max], [mod_min * lr['slope'] + lr['intercept'],
                mod_max * lr['slope'] + lr['intercept']],
                color='k', linestyle='-', linewidth=2, label='Linear fit')

        ## plot CI's for slope
        ax.plot([mod_min, mod_max], [mod_min * lr['slope_CI'][0] + lr['intercept_CI'][0],
                                     mod_max * lr['slope_CI'][0] + lr['intercept_CI'][0]],
                 color='r', linestyle='--', linewidth=2)
        ax.plot([mod_min, mod_max], [mod_min * lr['slope_CI'][1] + lr['intercept_CI'][1],
                                     mod_max * lr['slope_CI'][1] + lr['intercept_CI'][1]],
                 color='r', linestyle='--', linewidth=2, label='Slope CI')

        ## plot CI's for predictands
        ax.plot([mod_min, mod_max], [mod_min * lr['slope'] + upper_intercept,
                                     mod_max * lr['slope'] + upper_intercept],
                 color='g', linestyle='--', linewidth=2)
        ax.plot([mod_min, mod_max], [mod_min * lr['slope'] + lower_intercept,
                                     mod_max * lr['slope'] + lower_intercept],
                 color='g', linestyle='--', linewidth=2, label='Predictand CI')

        ax.set_xlabel('Modeled Data')
        ax.set_ylabel('Observed Data')
        fig.suptitle('Modeled vs. Observed {}: Linear Fit'.format(self.kind))
        plt.legend(loc='lower right', shadow=True)

        r_string = 'R Squared: {}'.format(np.around(lr['r_2'], decimals=3))
        plt.title(r_string)

        # Pretty plot
        seaborn.set(style="darkgrid")
        color = seaborn.color_palette()[2]
        g = seaborn.jointplot("model", "observed", data=df, kind="reg",
                              xlim=(df.model.min(), df.model.max()),
                              ylim=(df.observed.min(), df.observed.max()),
                              color=color, size=7)
        plt.suptitle('Modeled vs. Observed {}: Linear Fit'.format(self.kind))

        if save:
            fig.savefig(out_f)
        else:
            fig.show()

    def plotData(self, graph='time', save=False, out_f='', debug=False):
        """
        Provides a visualization of the data.

        Takes an option which determines the kind of graph to be made.
        time: plots the model data against the observed data over time
        scatter : plots the model data vs. observed data

        If save is set to True, saves the image file in out_f.
        """
        #define figure frame
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        ax = fig.add_subplot(111)

        if (graph == 'time'):
            ax.plot(self.times, self.model, label='Model Predictions')
            ax.plot(self.times, self.observed, color='r',
                     label='Observed Data')
            ax.set_xlabel('Time')
            if self.kind == 'elevation':
                ax.set_ylabel('Elevation (m)')
            if self.kind == 'speed':
                ax.set_ylabel('Flow speed (m/s)')
            if self.kind == 'direction':
                ax.set_ylabel('Flow direction (deg.)')
            if self.kind == 'u velocity':
                ax.set_ylabel('U velocity (m/s)')
            if self.kind == 'v velocity':
                ax.set_ylabel('V velocity (m/s)')
            if self.kind == 'velocity':
                ax.set_ylabel('Signed flow speed (m/s)')
            if self.kind == 'power density':
                ax.set_ylabel('Power density (W/m2)')

            fig.suptitle('Predicted and Observed {}'.format(self.kind))
            ax.legend(shadow=True)

        if (graph == 'scatter'):
            ax.scatter(self.model, self.observed, c='b', alpha=0.5)
            ax.set_xlabel('Predicted Height')
            ax.set_ylabel('Observed Height')
            fig.suptitle('Predicted vs. Observed {}'.format(self.kind))

        if save:
            fig.savefig(out_f)
        else:
            fig.show()

    def save_data(self):
            df = pd.DataFrame(data={'time': self.times.flatten(),
                                    'observed':self.observed.flatten(),
                                    'modeled':self.model.flatten() })
            df.to_csv(str(self.kind)+'.csv')

