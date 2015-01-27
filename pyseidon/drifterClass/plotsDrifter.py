#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn

class PlotsDrifter:
    """
    Description:
    -----------
    'Plots' subset of Drifter class gathers plotting functions
    """
    def __init__(self, variable, debug):
        self._debug = debug
        self._var = variable

    def trajectories(self, title='Drifter trajectories & speed (m/s)',
                     cmin=[], cmax=[], debug=False):
        '''
        2D xy colormap plot of all the trajectories.
        Colors represent the drifter velocity

        Keywords:
        --------
          - title = plot title, string
          - cmin = minimum limit colorbar
          - cmax = maximum limit colorbar
        '''   
        debug = debug or self._debug
        if debug: print "Plotting drifter's trajectories..."
        lon = self._var.lon
        lat = self._var.lat

        #Compute drifter's speeds
        u=self._var.u
        v=self._var.v
        norm=np.sqrt(u**2.0 + v**2.0) 

        #setting limits and levels of colormap
        if cmin==[]:
            if debug:
                print "Computing cmin..."
            cmin=norm[:].min()
        if cmax==[]:
            if debug:
                print "Computing cmax..."
            cmax=norm[:].max()

        #Figure window params
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        self._fig = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))

        #Scatter plot     
        #sc = plt.scatter(lon, lat, c=norm, lw=0, cmap=plt.cm.gist_earth)
        #sc.set_clim([cmin,cmax])

        #Quiver plot     
        sc = plt.quiver(lon, lat, u, v, norm, lw=0.0, scale=100.0,
                        cmap=plt.cm.gist_earth)
        sc.set_clim([cmin,cmax])

        #Label and axis parameters
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.gca().patch.set_facecolor('0.5')
        cbar=fig.colorbar(sc,ax=self._fig)
        cbar.set_label(title, rotation=-90,labelpad=30)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        self._fig.xaxis.set_major_formatter(ticks)
        self._fig.yaxis.set_major_formatter(ticks)
        plt.show()
        if debug or self._debug:
            print '...Passed'

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
