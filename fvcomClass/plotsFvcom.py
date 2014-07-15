#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn

class PlotsFvcom:
    """'Plots' subset of FVCOM class gathers plotting functions"""
    def __init__(self, cls):
        self._var = cls.Variables
        self._grid = cls.Grid
        self._debug = cls._debug

    def colormap_var(self, var, title='Title', cmin=[], cmax=[], mesh=True, debug=False):
        '''
        2D xy colormap plot of any given variable and mesh.
        Input:
        -----
          var = variable located at nodes, 1 dimensional numpy array

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
        lim = self._grid.region_n
        if dim==self._grid.nele:
            lim = self._grid.region_e
        bb = self._grid.ax
        
        
        # Mesh triangle
        nv = self._grid.nv.T -1
        if var.shape[0] == self._grid.nele:
            lon = self._grid.lonc
            lat = self._grid.latc
        else var.shape[0] == self._grid.node:     
            lon = self._grid.lon
            lat = self._grid.lat
        tri = Tri.Triangulation(lon, lat, triangles=nv)

        #setting limits and levels of colormap
        if cmin==[]:
            cmin=round(var[lim].min())
            if cmin==0.0:
                cmin = -1.0
        if cmax==[]:
            cmax=round(var[lim].max())
            if cmax==0.0:
                cmax==1.0
        step = round(round(cmax-cmin) / 50.0)
        levels=np.arange(cmin, (cmax+step), step)   # depth contours to plot

        #Figure window params
        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))

        #Plotting functions
        plt.tricontourf(tri, var,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
        if mesh:
            plt.triplot(tri)

        #Label and axis parameters
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.gca().patch.set_facecolor('0.5')
        cbar=plt.colorbar()
        cbar.set_label(title, rotation=-90,labelpad=30)
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        ax.xaxis.set_major_formatter(ticks)
        ax.yaxis.set_major_formatter(ticks)
        ax.set_xlim([bb[0],bb[1]])
        ax.set_ylim([bb[2],bb[3]])
        plt.grid()
        plt.show()
        if debug or self._debug:
            print '...Passed'

    #def spatial

