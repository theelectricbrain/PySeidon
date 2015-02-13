#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import seaborn
from windrose import WindroseAxes
from interpolation_utils import *
from miscellaneous import depth_at_FVCOM_element as depth_at_ind

class PlotsFvcom:
    """
    Description:
    -----------
    'Plots' subset of FVCOM class gathers plotting functions
    """
    def __init__(self, variable, grid, debug):
        self._debug = debug
        self._var = variable
        self._grid = grid
        #Back pointer
        grid = self._grid
        #self._grid._ax = grid._ax

    def _def_fig(self):
        """Defines figure window"""
        self._fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')


    def colormap_var(self, var, title='Title', cmin=[], cmax=[], cmap=[],
                     mesh=True, debug=False):
        '''
        2D xy colormap plot of any given variable and mesh.

        Input:
        -----
          - var = gridded variable, 1 D numpy array (nele or nnode)

        Keywords:
        --------
          - title = plot title, string
          - cmin = minimum limit colorbar
          - cmax = maximum limit colorbar
          - cmap = matplolib colormap
          - mesh = True, with mesh; False, without mesh
        '''
        debug = debug or self._debug
        if debug:
            print 'Plotting grid...'
        # Figure if var had nele or nnode dimensions
        if var.shape[0] == self._grid.nele:
            dim = self._grid.nele
        elif var.shape[0] == self._grid.nnode:
            dim = self._grid.nnode
        else:
            print "Var has the wrong dimension, var.shape[0]= Grid.nele or nnode"
            return

        # Bounding box nodes, elements and variable 
        lon = self._grid.lon[:]
        lat = self._grid.lat[:]
        if debug:
            print "Computing bounding box..."
        if self._grid._ax == []:
            self._grid._ax = [lon.min(), lon.max(),
                             lat.min(), lat.max()]
        bb = self._grid._ax  

        if not hasattr(self._grid, 'triangle'):        
            # Mesh triangle
            if debug:
                print "Computing triangulation..."
            trinodes = self._grid.trinodes[:] 
            tri = Tri.Triangulation(lon, lat, triangles=trinodes)
            self._grid.triangle = tri
        else:
            tri = self._grid.triangle

        #setting limits and levels of colormap
        if cmin==[]:
            if debug:
                print "Computing cmin..."
            cmin=var[:].min()
        if cmax==[]:
            if debug:
                print "Computing cmax..."
            cmax=var[:].max()
        step = (cmax-cmin) / 50.0
        levels=np.arange(cmin, (cmax+step), step)   # depth contours to plot

        #Figure window params
        self._def_fig()
        self._ax = self._fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))

        #Plotting functions
        if debug:
            print "Computing colormap..."
        if cmap==[]:
            f = self._ax.tripcolor(tri, var[:],vmax=cmax,vmin=cmin,cmap=plt.cm.gist_earth)
        else:
            f = self._ax.tripcolor(tri, var[:],vmax=cmax,vmin=cmin,cmap=cmap) 
        if mesh:
            plt.triplot(tri, color='white', linewidth=0.5)

        #Label and axis parameters
        self._ax.set_ylabel('Latitude')
        self._ax.set_xlabel('Longitude')
        self._ax.patch.set_facecolor('0.5')
        cbar=self._fig.colorbar(f, ax=self._ax)
        cbar.set_label(title, rotation=-90,labelpad=30)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        self._ax.xaxis.set_major_formatter(ticks)
        self._ax.yaxis.set_major_formatter(ticks)
        self._ax.set_xlim([bb[0],bb[1]])
        self._ax.set_ylim([bb[2],bb[3]])
        self._ax.grid()
        self._fig.show()
        if debug or self._debug:
            print '...Passed'

    def rose_diagram(self, direction, norm):

        """
        Plots rose diagram

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
        #fig = plt.figure(figsize=(18,10))
        #plt.rc('font',size='22')
        self._def_fig()      
        rect = [0.1, 0.1, 0.8, 0.8]
        ax = WindroseAxes(self._fig, rect)#, axisbg='w')
        self._fig.add_axes(ax)
        #Rose
        ax.bar(direction, norm , normed=True, opening=0.8, edgecolor='white')
        #adjust legend
        l = ax.legend(shadow=True, bbox_to_anchor=[-0.1, 0], loc='lower left')
        plt.setp(l.get_texts(), fontsize=10)
        plt.xlabel('Rose diagram in % of occurrences - Colormap of norms')
        self._fig.show() 

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
        #fig = plt.figure(figsize=(18,10))
        #plt.rc('font',size='22')
        self._def_fig()
        self._ax = self._fig.add_subplot(111)         
        self._ax.plot(x, y, label=title)
        scale = 1
        self._ax.set_ylabel(yLabel)
        self._ax.set_xlabel(xLabel)
        self._ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
	self._ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
	self._ax.grid(b=True, which='major', color='w', linewidth=1.5)
	self._ax.grid(b=True, which='minor', color='w', linewidth=0.5)
        if not yerror==[]:
            self._ax.fill_between(x, y-yerror, y+yerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)        
        if not xerror==[]:
            self._ax.fill_betweenx(y, x-xerror, x+xerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)
        if (not xerror==[]) or (not yerror==[]): 
            blue_patch = mpatches.Patch(color='#089FFF',
                         label='Standard deviation',alpha=0.2)
            plt.legend(handles=[blue_patch],loc=1, fontsize=12)

        self._fig.show()      

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
        ## the histogram of the data
        #fig = plt.figure(figsize=(18,10))
        self._def_fig()
        self._ax = self._fig.add_subplot(111)        
        density, bins = np.histogram(y, bins=50, normed=True, density=True)
        unity_density = density / density.sum()
        widths = bins[:-1] - bins[1:]
        # To plot correct percentages in the y axis 
        self._ax.bar(bins[1:], unity_density, width=widths)
        formatter = ticker.FuncFormatter(lambda v, pos: str(v * 100))
        self._ax.yaxis.set_major_formatter(formatter)

        plt.ylabel(yLabel)
        plt.xlabel(xLabel)

        self._fig.show()  

    def add_points(self, x, y, label=' ', color='black'):
        """
        Adds scattered points (x,y) on current figure,
        where x and y are 1D arrays of the same lengths.

        Inputs:
        ------
          - x = float number or list of float numbers
          - y = float number or list of float numbers

        Keywords:
        --------
          - Label = a string
          - Color = a string, 'red', 'green', etc. or gray shades like '0.5' 
        """
        plt.scatter(x, y, s=50, color=color)
        #TR : annotate does not work on my machine !?
        plt.annotate(label, xy=(x, y), xycoords='data', xytext=(-20, 20),
                     textcoords='offset points', ha='right',
                     arrowprops=dict(arrowstyle="->", shrinkA=0),
                     fontsize=12)

    def _save_as_pickle(ax, filename='saved_plot.p'):
        """
        Saves figure as pickle file which can be then re-loaded 
        with: ax = pickle.load(file('saved_plot.p'))
        """ 
        pickle.dump(ax, file(filename, 'w'))
                                
