#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np

class Functions:
    ''' '''
    def __init__(self,cls):
        self._var = cls.Variables
        self._grid = cls.Grid
        self._ax = cls._ax
        self._debug = cls._debug
    
    def el_region(self, debug=False):
        '''Use/Inputs/Outputs of this method has to be clarified !!!'''
        
        if debug or self._debug:
            print 'Computing region_e...'

        region_e = np.argwhere((self._var.lonc >= self._ax[0]) &
                                     (self._var.lonc <= self._ax[1]) &
                                     (self._var.latc >= self._ax[2]) &
                                     (self._var.latc <= self._ax[3]))
        if debug or self._debug:
            print '...Passed'

        return region_e

    def node_region(self, debug=False):
        '''Use/Inputs/Outputs of this method has to be clarified !!!'''
        if debug or self._debug:
            print 'Computing region_n...'

        region_n = np.argwhere((self._var.lon >= self._ax[0]) &
                                     (self._var.lon <= self._ax[1]) &
                                     (self._var.lat >= self._ax[2]) &
                                     (self._var.lat <= self._ax[3]))
        if debug or self._debug:
            print '...Passed'

        return region_n

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
            hc[i] = np.mean(self._var.h[v], axis=1)

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
        #No return ???

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

