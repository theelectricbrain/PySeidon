#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np
import numexpr as ne

class FunctionsFvcom:
    ''''Utils' subset of FVCOM class gathers useful functions""" '''
    def __init__(self, cls):
        self._debug = cls._debug
        self._var = cls.Variables
        self._grid = cls.Grid
        #Create pointer to FVCOM class
        cls.Variables = self._var
        cls.Grid = self._grid
  
    def centers(self, debug=False):
        '''Use/Inputs/Outputs of this method has to be clarified !!!'''
        if debug or self._debug:
            print 'Computing elc, hc...'

        size = self._grid.trinodes.T.shape[0]
        size1 = self._var.el.shape[0]
        elc = np.zeros((size1, size))
        hc = np.zeros((size))
        for i,v in enumerate(self._grid.trinodes.T):
            elc[:, i] = np.mean(self._var.el[:, v], axis=1)
            hc[i] = np.mean(self._grid.h[v], axis=1)

        if debug or self._debug:
            print '...Passed'

        return elc, hc

    def hc(self, debug=False):
        '''Use/Inputs/Outputs of this method has to be clarified !!!'''
        if debug or self._debug:
            print 'Computing hc...'

        size = self._grid.trinodes.T.shape[0]
        size1 = self._var.el.shape[0]
        elc = np.zeros((size1, size))
        for i,v in enumerate(self._grid.trinodes.T):
            elc[:, i] = np.mean(self._var.el[:, v], axis=1)

        if debug or self._debug:
            print '...Passed'
        #TR_comments: No return ???

    def closest_point(self, pt_lon, pt_lat,debug=False):
        '''
        Finds the closest exact lon, lat centre indexes of an FVCOM class
        to given lon, lat coordinates.

        Inputs:
          - pt_lon = list of longitudes in degrees
          - pt_lat = list of latitudes in degrees
        Outputs:
          - closest_point_indexes = numpy array of grid indexes
        '''
        if debug or self._debug:
            print 'Computing closest_point_indexes...'

        points = np.array([pt_lon, pt_lat]).T

        #point_list = np.array([lon,lat]).T
        point_list = np.array([self._var.lonc, self._var.latc]).T

        closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                        (point_list[:, 1] - points[:, 1, None])**2)

        closest_point_indexes = np.argmin(closest_dist, axis=1)

        if debug or self._debug:
            print '...Passed'

        return closest_point_indexes

    def ele_region(self, debug=False):
        '''Return element indexes included in bounding box, aka ax'''
        
        if debug or self._debug:
            print 'Computing region_e...'

        region_e = np.argwhere((self._grid.lonc >= self._grid.ax[0]) &
                                     (self._grid.lonc <= self._grid.ax[1]) &
                                     (self._grid.latc >= self._grid.ax[2]) &
                                     (self._grid.latc <= self._grid.ax[3]))

        # Create new grid variable through pointer
        self._grid.region_e = region_e[:,0]

        if debug or self._debug:
            print '...Passed'

        return region_e

    def node_region(self, debug=False):
        '''Return node indexes included in bounding box, aka ax'''
        if debug or self._debug:
            print 'Computing region_n...'

        region_n = np.argwhere((self._grid.lon >= self._grid.ax[0]) &
                                     (self._grid.lon <= self._grid.ax[1]) &
                                     (self._grid.lat >= self._grid.ax[2]) &
                                     (self._grid.lat <= self._grid.ax[3]))

        # Create new grid var through pointer
        self._grid.region_n = region_n[:,0]

        if debug or self._debug:
            print '...Passed'

        return region_n

    def bounding_box(self, ax=[]):
        """
        Define bounding box and reset the box by default.
        Input ex:
        --------
          .bounding_box(ax=[min lon, max lon, min lat, max lat])
        """
        # reset through pointer
        if ax:
            self._grid.ax = ax
        else:
            self._grid.ax = [min(self._grid.lon), max(self._grid.lon),
                             min(self._grid.lat), max(self._grid.lat)]
        self.node_region()
        self.ele_region()
          


    def velo_norm(self, debug=False):
        """Compute velocity norm -> FVCOM.Variables.velo_norm"""
        if debug or self._debug:
            print 'Computing velocity norm...'
        if self._var._3D:
            u = self._var.u[:, :, self._grid.region_e[:]]
            v = self._var.v[:, :, self._grid.region_e[:]]
            ww = self._var.ww[:, :, self._grid.region_e[:]]
            vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')
        else:
            u = self._var.ua[:, self._grid.region_e[:]]
            v = self._var.va[:, self._grid.region_e[:]]
            vel = ne.evaluate('sqrt(u**2 + v**2)')

        self._var.velo_norm = vel           

        if debug or self._debug:
            print '...Passed'
        

