#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn
import pandas as pd

class PlotsDrifter:
    """
    **'Plots' subset of Drifter class gathers plotting functions**
    """
    def __init__(self, variable, debug):
        self._debug = debug
        # Pointer
        setattr(self, '_var', variable)

        return

    def _def_fig(self):
        """Defines figure window"""
        self._fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')

    def trajectories(self, title='Drifter trajectories & speed (m/s)',
                     cmin=[], cmax=[], debug=False, dump=False, **kwargs):
        """
        2D xy colormap plot of all the trajectories.
        Colors represent the drifter velocity

        Options:
          - title = plot title, string
          - cmin = minimum limit colorbar
          - cmax = maximum limit colorbar
          - dump = boolean, dump profile data in csv file
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
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
        self._def_fig()
        self._ax = self._fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))

        #Scatter plot     
        #sc = plt.scatter(lon, lat, c=norm, lw=0, cmap=plt.cm.gist_earth)
        #sc.set_clim([cmin,cmax])

        #Quiver plot     
        sc = plt.quiver(lon, lat, u, v, norm, lw=0.0, scale=100.0,
                        cmap=plt.cm.gist_earth)
        sc.set_clim([cmin,cmax])

        #Label and axis parameters
        self._ax.set_ylabel('Latitude')
        self._ax.set_xlabel('Longitude')
        plt.gca().patch.set_facecolor('0.5')
        cbar=self._fig.colorbar(sc,ax=self._ax)
        cbar.set_label(title, rotation=-90,labelpad=30)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        self._ax.xaxis.set_major_formatter(ticks)
        self._ax.yaxis.set_major_formatter(ticks)
        self._ax.grid()
        self._fig.show()
        if dump:
            self._dump_data_as_csv(norm, u, v, lon, lat, title='drifter_velocity', **kwargs)
        if debug or self._debug:
            print '...Passed'

    def _dump_data_as_csv(self, var1, var2, var3, x, y, title=' ', **kwargs):
        """
        Dumps map data in csv file

        Inputs:
          - var1 = 1 D numpy array oh n elements
          - var2 = 1 D numpy array oh n elements
          - var3 = 1 D numpy array oh n elements
          - x = coordinates, 1 D numpy array (nele or nnode)
          - y = coordinates, 1 D numpy array (nele or nnode)

        Options:
          - title = file name, string
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        if title == ' ': title = 'dump_map_data'
        filename=title + '.csv'
        df = pd.DataFrame({'norm': var1, 'u':var2, 'v':var3,
                           'lon':x[:], 'lat':y[:]})

        df.to_csv(filename, encoding='utf-8', **kwargs)

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
