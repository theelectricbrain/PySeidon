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


