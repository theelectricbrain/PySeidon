#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn
from windrose import WindroseAxes
from interpolation_utils import *
from shortest_element_path import *

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
          var = gridded variable, 1 dimensional numpy array

        Keywords:
        --------
          cmin = min limit colorbar
          cmax = max limit colorbar
          mesh = True, with mesh; False, without mesh
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
          Label = a string
          Color = a string, 'red', 'green', etc. or gray shades like '0.5' 
        """
        plt.scatter(x, y, s=100, color=color)
        #TR : annotate does not work on my machine !?
        plt.annotate(label, xy=(x, y), xycoords='data', xytext=(-20, 20),
                     textcoords='offset points', ha='right',
                     arrowprops=dict(arrowstyle="->", shrinkA=0))

    def vertical_slice(self, var, start_pt, end_pt, cmax=[], cmin=[], debug=False):
        """
        Draw vertical slice in var along the shortest path between
        start_point, end_pt.
 
        Inputs:
        ------
          -var = 2D dimensional (sigma level, element) variable, array
          -start_pt = starting point, [longitude, latitude]
          -end_pt = ending point, [longitude, latitude]
        """
        debug = debug or self._debug
        if not self._var._3D:
            print "Error: Only available for 3D runs."
            raise
        else: 
            lons = [start_pt[0], end_pt[0]]
            lats = [start_pt[1], end_pt[1]]
            #Finding the closest elements to start and end points
            ind = closest_point(lons, lats, self._grid.lonc, self._grid.latc, debug)
            #Finding the shortest path between start and end points
            print "Computing shortest path..."
            short_path = shortest_element_path(self._grid.lonc,
                                               self._grid.latc,
                                               self._grid.lon,
                                               self._grid.lat,
                                               self._grid.trinodes,
                                               self._grid.h)
            el, _ = short_path.getTargets([ind])
            #Extract along path
            varP = var[:,el]
            x = self._grid.xc[el]
            y = self._grid.yc[el]
            siglay = self._grid.siglay[:, el]
            line = np.zeros(siglay.shape)
            #Pythagore
            line[:,1:] = np.sqrt(np.square(x[1:]-x[:-1]) + np.square(y[1:]-y[:-1]))
            return el, varP, siglay, line
            #Computing depth
            #print "Computing depth..."
            #TR comment: the interpolation here is not perfect but doesn't
            #            matter too much for plotting purposes
            #size = self._grid.trinodes[el].shape[0]
            #size1 = self._var.el.shape[0]
            #elc = np.zeros((size1, size))
            #hc = np.zeros((size))
            #for ind,value in enumerate(self._grid.trinodes.T[el[0]]):
            #    elc[:, ind] = np.mean(self._var.el[:, value-1], axis=1)
            #    hc[ind] = np.mean(self._grid.h[value-1])
            #Plot features
            #if not cmax:
            #    cmax = np.max(var)
            #if not cmin:
            #    cmin = np.min(var)
            #plt.clf()
            #fig,ax = plt.subplots()
            #plt.rc('font',size='22')
            #levels = np.linspace(0,3.3,34)
            #cs = ax.contourf(line,siglay,var,levels=levels, cmap=plt.cm.jet)
            #ax.contour(line,siglay,mean_vel,cs.levels, colors='k')
            #cbar = fig.colorbar(cs,ax=ax)
            #cbar.set_label(r'Velocity $(m/s)$', rotation=-90,labelpad=30)
            #if eastwest:
            #    ax.set_xlabel('Longitude')
            #else:
            #    ax.set_xlabel('Latitude')
            #ax.set_title('vel_mean')
            #scale = 1
            #ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
            #ax.xaxis.set_major_formatter(ticks)
            #ax.yaxis.set_major_formatter(ticks)                                   
