from scipy.interpolate import interp1d
import numpy as np
from datetime import timedelta
import time


def interpol(data_1, data_2, time_step=timedelta(minutes=5)):
    '''
    Interpolates between two datasets so their points line up in the time
    domain.

    Accepts two sets of data, each of which are dictionaries containing two
    values:
    time-   array containing the datetimes corresponding to the points
    pts-    1D numpy array containing the data

    Third optional argument sets the time between data points in the output
    data. Is a timedelta object, defaults to 10 minutes.
    '''
    dt_1 = data_1['time']
    dt_2 = data_2['time']

    # create POSIX timestamp array corresponding to each dataset
    times_1, times_2 = np.zeros(len(dt_1)), np.zeros(len(dt_2))
    for i in np.arange(times_1.size):
        times_1[i] = time.mktime(dt_1[i].timetuple())
    for v in np.arange(times_2.size):
        times_2[v] = time.mktime(dt_2[v].timetuple())

    # generate interpolation functions using linear interpolation
    f1 = interp1d(times_1, data_1['pts'])
    f2 = interp1d(times_2, data_2['pts'])

    # choose interval on which to interpolate
    start = max(times_1[0], times_2[0])
    end = min(times_1[-1], times_2[-1])
    length = end - start

    # determine number of steps in the interpolation interval
    step_sec = time_step.total_seconds()
    steps = int(length / step_sec)

    # create POSIX timestamp array for new data and perform interpolation
    output_times = start + np.arange(steps) * step_sec

    series_1 = f1(output_times)
    series_2 = f2(output_times)

    dt_start = max(dt_1[0], dt_2[0])

    return (series_1, series_2, time_step, dt_start)
