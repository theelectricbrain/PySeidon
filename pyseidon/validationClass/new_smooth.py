#!/usr/bin/python2.7
# encoding: utf-8
from datetime import timedelta
import numpy as np
import time
from scipy.ndimage.measurements import find_objects


def smooth(data_1, dt_1, data_2, dt_2, debug=False, debug_plot=False):
    '''
    Smooths a dataset by taking the average of all datapoints within
    a certain timestep to reduce noise. Lines up two datasets in the
    time domain, as well.
    Accepts four variables representing the data. data_1 and data_2 are the
    data points, dt_1 and dt_2 are the datetimes corresponding to the points.
    '''
    if debug: print "smooth..."

    time_step = timedelta(minutes=10)

    # create POSIX timestamp array corresponding to each dataset
    '''
    times_1, times_2 = np.zeros(len(dt_1)), np.zeros(len(dt_2))
    for i in np.arange(times_1.size):
        times_1[i] = time.mktime(dt_1[i].timetuple())
    for i in np.arange(times_2.size):
        times_2[i] = time.mktime(dt_2[i].timetuple())
    '''
    make_posix = lambda x: time.mktime(x.timetuple())
    times_1 = map(make_posix, dt_1)
    times_2 = map(make_posix, dt_2)
    times_1, times_2 = np.array(times_1), np.array(times_2)

    # choose smoothing interval
    start = max(times_1[0], times_2[0])
    end = min(times_1[-1], times_2[-1])
    length = end - start
    dt_start = max(dt_1[0], dt_2[0])

    # grab number of steps and timestamp for start time
    step_sec = time_step.total_seconds()
    steps = int(length / step_sec) + 1

    # sort times into bins then identify bin indices
    series_1, series_2 = np.zeros(steps), np.zeroes(steps)
    time_bins = np.arange(steps) * step_sec + start
    inds_1 = np.digitize(times_1, time_bins)
    inds_2 = np.digitize(times_2, time_bins)
    slices_1 = find_objects(inds_1)
    slices_2 = find_objects(inds_2)

    # should be a faster for loop than before, iterate through the slices
    for i, slc in enumerate(slices_1):
        series_1[i] = np.nanmean(data_1[slc])
    for i, slc in enumerate(slices_2):
        series_2[i] = np.nanmean(data_2[slc])

    if debug: print "...smooth done."
    return (series_1, series_2, time_step, dt_start)
