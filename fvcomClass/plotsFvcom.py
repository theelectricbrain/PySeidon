#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn
from windrose import WindroseAxes
from interpolation_utils import *
from miscellaneous import depth_at_FVCOM_element as depth_at_ind

class PlotsFvcom:
    """'Plots' subset of FVCOM class gathers plotting functions"""
    def __init__(self, variable, grid, debug):
        self._debug = debug
        self._var = variable
        self._grid = grid
    def colormap_var(self, var, title='Title', cmin=[], cmax=[], mesh=True, debug=False):
        '''
        2D xy colormap plot of any given variable and mesh.

        Input:
        -----
          - var = gridded variable, 1 dimensional numpy array

        Keywords:
        --------
          - title = plot title, string
          - cmin = minimum limit colorbar
          - cmax = maximum limit colorbar
          - mesh = True, with mesh; False, without mesh
        '''
        if debug or self._debug:
            print 'Plotting grid...'
        # Figure if var had nele or node dimensions
        if var.shape[0] == self._grid.nele:
            dim = self._grid.nele
        elif var.shape[0] == self._grid.node:
            dim = self._grid.node
        else:
            print "Var has the wrong dimension, var.shape[0]= Grid.nele or node"
            return

        # Bounding box nodes, elements and variable
        bb = self._grid._ax     
        
        # Mesh triangle
        trinodes = self._grid.trinodes 
        lon = self._grid.lon
        lat = self._grid.lat
        tri = Tri.Triangulation(lon, lat, triangles=trinodes)

        #setting limits and levels of colormap
        if cmin==[]:
            cmin=var[:].min()
        if cmax==[]:
            cmax=var[:].max()
        step = (cmax-cmin) / 50.0
        levels=np.arange(cmin, (cmax+step), step)   # depth contours to plot

        #Figure window params
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        self._fig = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
        #CS = self._fig.tricontour(tri, -1.0*self._grid.h, levels=levels,
        #                          shading='faceted',cmap=plt.cm.gist_earth)

        #Plotting functions
        f = self._fig.tripcolor(tri, var,vmax=cmax,vmin=cmin,cmap=plt.cm.gist_earth)
        if mesh:
            plt.triplot(tri)

        #Label and axis parameters
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.gca().patch.set_facecolor('0.5')
        cbar=fig.colorbar(f, ax=self._fig)
        cbar.set_label(title, rotation=-90,labelpad=30)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        self._fig.xaxis.set_major_formatter(ticks)
        self._fig.yaxis.set_major_formatter(ticks)
        self._fig.set_xlim([bb[0],bb[1]])
        self._fig.set_ylim([bb[2],bb[3]])
        plt.grid()
        plt.show()
        if debug or self._debug:
            print '...Passed'

    def rose_diagram(self, direction, norm):

        """Plot rose diagram. direction and norm = 1D arrays"""

        #Create new figure
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        rect = [0.1, 0.1, 0.8, 0.8]
        ax = WindroseAxes(fig, rect)#, axisbg='w')
        fig.add_axes(ax)
        #Rose
        ax.bar(direction, norm , normed=True, opening=0.8, edgecolor='white')
        #adjust legend
        l = ax.legend(shadow=True)
        plt.setp(l.get_texts(), fontsize=10)
        plt.xlabel('Rose diagram in % of occurrences - Colormap of norms')
        plt.show() 

    def plot_xy(self, x, y, title=' ', xLabel=' ', yLabel=' '):
        """Simple Y vs X plot"""
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        self._fig = plt.plot(x, y)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        plt.ylabel(yLabel)
        plt.xlabel(xLabel)
        plt.show()      

    def add_points(self, x, y, label=' ', color='black'):
        """
        Add scattered points (x,y) on current figure,
        where x and y are 1D arrays of the same lengths.

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

                                
