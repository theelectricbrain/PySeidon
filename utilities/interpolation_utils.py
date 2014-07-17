#!/usr/bin/python2.7
# encoding: utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
from matplotlib.path import Path

def closest_point( pt_lon, pt_lat, lon, lat, debug=False):
    '''
    Finds the closest exact lon, lat centre indexes of an FVCOM class
    to given lon, lat coordinates.

    Inputs:
      - pt_lon = list of longitudes in degrees to find
      - pt_lat = list of latitudes in degrees to find
      - lon = list of longitudes in degrees to search in
      - lat = list of latitudes in degrees to search in
    Outputs:
      - closest_point_indexes = numpy array of grid indexes
    '''
    if debug:
        print 'Computing closest_point_indexes...'

    points = np.array([pt_lon, pt_lat]).T
    point_list = np.array([lon, lat]).T

    closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                    (point_list[:, 1] - points[:, 1, None])**2)

    closest_point_indexes = np.argmin(closest_dist, axis=1)
    if debug:
        print '...Passed'

    return closest_point_indexes

def interp_at_point(var, pt_lon, pt_lat, lon, lat, index, trinodes , debug=False):
    """
    Interpol any given variables at any give location.
    Inputs:
      - var = variable, numpy array, dim=(time, nele or node)
      - pt_lon = longitude in degrees to find
      - pt_lat = latitude in degrees to find
      - lon = list of longitudes of var, numpy array, dim=(nele or node)
      - lat = list of latitudes of var, numpy array, dim=(nele or node)
      - trinodes = FVCOM trinodes, numpy array, dim=(3,nele)

    Outputs:
      - varInterp = var interpolate at (pt_lon, pt_lat)
    """
    if debug:
        print 'Computing var_interp...'

    #Triangulation
    triIndex = trinodes[index]
    if debug:
        print triIndex, lon[triIndex], lat[triIndex]
    newtri = Tri.Triangulation(lon[triIndex], lat[triIndex], np.array([[0,1,2]]))
    trif = newtri.get_trifinder()
    trif.__call__(pt_lon, pt_lat)
    if debug:
        if len(var)==1:
            averEl = var[triIndex]
            print 'Var', averEl
            inter = Tri.LinearTriInterpolator(newtri, averEl)
            zi = inter(pt_lon, pt_lat)
            print zi
            plt.tricontourf(newtri, averEl, cmap=plt.cm.jet)
            cbar = plt.colorbar()
        fig = plt.figure(figsize=(18,10))
        ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat[triIndex])*np.pi/180.0)))

        plt.triplot(newtri, 'ko-')
        plt.plot(pt_lon, pt_lat, 'ko-')

        plt.show()

    #Choose the right interpolation depending on the variable
    if len(var.shape)==1:
        #triVar = np.zeros(triIndex.shape)
        triVar = var[triIndex]
        inter = Tri.LinearTriInterpolator(newtri, triVar[:])
        varInterp = inter(pt_lon, pt_lat)      
    elif len(var.shape)==2:
        #triVar = np.zeros((var.shape[0], triIndex.shape[0]))
        triVar = var[:, triIndex]
        varInterp = np.ones(triVar.shape[0])
        for i in range(triVar.shape[0]):
            inter = Tri.LinearTriInterpolator(newtri, triVar[i,:])
            varInterp[i] = inter(pt_lon, pt_lat)       
    else:
        #triVar = np.zeros((var.shape[0], var.shape[1], triIndex.shape[0]))
        triVar = var[:, :, triIndex]
        varInterp = np.ones(triVar.shape[:-1])
        for i in range(triVar.shape[0]):
           for j in range(triVar.shape[1]):
               inter = Tri.LinearTriInterpolator(newtri, triVar[i,j,:])
               varInterp[i,j] = inter(pt_lon, pt_lat)

    if debug:
        print '...Passed'

    return varInterp


