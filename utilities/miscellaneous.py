#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from datetime import datetime
from datetime import timedelta

def date2py(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + \
        timedelta(days=matlab_datenum%1) - timedelta(days = 366)

    return python_datetime

def signal_extremum(signal):
    """
    This function spots the extremum of a random signal(x).
    Inputs:
      - signal: 1D array of n elements
    Outputs:
      - extremum: 1D array of n elements
      - indices: list containing extremum indices
    """

    extremum = np.zeros(signal.shape[0])
    indices = []

    S = np.sign(signal[0] - signal[1])
    N = (np.arange(signal.shape[0]-2)) + 1

    for i in N:
        E = np.sign(signal[i] - signal[i+1])
        if (E != S):
            extremum[i] = 1.0
            indices.append(i)
            S = np.sign(S*(-1.0))
    
    return extremum, indices

def exceedance_SSE(x, signal, extremum, indices):
    """
    This function calculate the excedence curve of a Sea Surface Elevation signal(x).
    Inputs:
      - x: 1D array of n elements
      - signal: 1D array of n elements
      - extremum: 1D array of n elements
      - indices: list containing extremum indices
    Outputs:
      - Exceedance: list of % of occurences
      - Ranges: list of signal amplitude bins
    """

    N = len(indices)
    period = []
    amp = []
    time_stamp = []

    for i in range(N-1):
        p = x[indices[i+1]]-x[indices[i]]
        if p > (60*60*3): #exceeds 3 hour ramping
            period.append(p)
            amp.append(abs(signal[indices[i]]-signal[indices[i+1]]))
            time_stamp.append(indices[i])

    Max = round(max(amp),1)	
    dy = round((Max/10.0),1)
    Ranges = np.arange(0,(Max + dy), dy)
    Exceedance = np.zeros(Ranges.shape[0])
    Period = np.sum(period)

    N = len(amp)
    M = len(Ranges)

    for i in range(M):
        r = Ranges[i]
        for j in range(N):
            if amp[j] > r:
                Exceedance[i] = Exceedance[i] + period[j]

    Exceedance = (Exceedance * 100) / Period

    return Exceedance, Ranges, amp, time_stamp
