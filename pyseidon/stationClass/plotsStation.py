#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn
from windrose import WindroseAxes
from interpolation_utils import *
from miscellaneous import depth_at_FVCOM_element as depth_at_ind

class PlotsStation:
    """
    Description:
    -----------
    'Plots' subset of Station class gathers plotting functions
    """
    def __init__(self, variable, grid, debug):
        self._debug = debug
        self._var = variable
        self._grid = grid
        #Back pointer
        grid = self._grid
        #self._grid._ax = grid._ax

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
        plt.legend()
        plt.show()      

    def add_points(self, x, y, label=' ', color='black'):
        """
        Add scattered points (x,y) on current figure,
        where x and y are 1D arrays of the same lengths.

        Inputs:
        ------
          - x = float number
          - y = float numbe

        Keywords:
        --------
          - Label = a string
          - Color = a string, 'red', 'green', etc. or gray shades like '0.5' 
        """
        plt.scatter(x, y, s=100, color=color)
        #TR : annotate does not work on my machine !?
        plt.annotate(label, xy=(x, y), xycoords='data', xytext=(-20, 20),
                     textcoords='offset points', ha='right',
                     arrowprops=dict(arrowstyle="->", shrinkA=0))


#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
