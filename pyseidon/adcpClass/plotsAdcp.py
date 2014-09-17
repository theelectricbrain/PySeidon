#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn

class PlotsAdcp:
    """'Plots' subset of FVCOM class gathers plotting functions"""
    def __init__(self, variable, debug=False):
        self._var = variable

    def plot_xy(self, x, y, title=' ', xLabel=' ', yLabel=' '):
        """
        Simple X vs Y plot

        Inputs:
        ------
          - x = 1D array
          - y = 1D array
        """
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        self._fig = plt.plot(x, y, label=title)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        plt.ylabel(yLabel)
        plt.xlabel(xLabel)
        #plt.legend()
        plt.show()      

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
