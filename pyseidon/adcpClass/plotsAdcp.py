#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn
from windrose import WindroseAxes
from interpolation_utils import *

class PlotsAdcp:
    """'Plots' subset of FVCOM class gathers plotting functions"""
    def __init__(self, variable, debug=False):
        self._var = variable

    def plot_xy(self, x, y, xerror=[], yerror=[], title=' ', xLabel=' ', yLabel=' '):
        """
        Simple X vs Y plot

        Inputs:
        ------
          - x = 1D array
          - y = 1D array

        Keywords:
        --------
          - xerror = error on 'x', 1D array
          - yerror = error on 'y', 1D array
          - title = plot title, string
          - xLabel = title of the x-axis, string
          - yLabel = title of the y-axis, string
        """
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        self._fig = plt.plot(x, y, label=title)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        plt.ylabel(yLabel)
        plt.xlabel(xLabel)
        if not yerror==[]:
            #plt.errorbar(x, y, yerr=yerror, fmt='o', ecolor='k')
            plt.fill_between(x, y-yerror, y+yerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)
        if not xerror==[]:
            #plt.errorbar(x, y, xerr=xerror, fmt='o', ecolor='k')
            plt.fill_betweenx(y, x-xerror, x+xerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)

        plt.show() 

    def Histogram(self, y, title=' ', xLabel=' ', yLabel=' '):
        """
        Histogram plot

        Inputs:
        ------
          - bins = list of bin edges
          - y = 1D array

        Keywords:
        --------
          - title = plot title, string
          - xLabel = title of the x-axis, string
          - yLabel = title of the y-axis, string
        """
        fig = plt.figure(figsize=(18,10))
        density, bins = np.histogram(y, bins=50, normed=True, density=True)
        unity_density = density / density.sum()
        widths = bins[:-1] - bins[1:]
        # To plot correct percentages in the y axis 
        plt.bar(bins[1:], unity_density, width=widths)
        formatter = ticker.FuncFormatter(lambda v, pos: str(v * 100))
        plt.gca().yaxis.set_major_formatter(formatter)

        plt.ylabel(yLabel)
        plt.xlabel(xLabel)

        plt.show()   
  
    def rose_diagram(self, direction, norm):

        """
        Plot rose diagram

        Inputs:
        ------
          - direction = 1D array
          - norm = 1D array
        """
        #Convertion
        #TR: not quite sure here, seems to change from location to location
        #    express principal axis in compass
        direction = np.mod(90.0 - direction, 360.0)

        #Create new figure
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        rect = [0.1, 0.1, 0.8, 0.8]
        ax = WindroseAxes(fig, rect)#, axisbg='w')
        fig.add_axes(ax)
        #Rose
        ax.bar(direction, norm , normed=True, opening=0.8, edgecolor='white')
        #adjust legend
        l = ax.legend(shadow=True, bbox_to_anchor=[-0.1, 0], loc='lower left')
        plt.setp(l.get_texts(), fontsize=10)
        plt.xlabel('Rose diagram in % of occurrences - Colormap of norms')
        plt.show() 

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
