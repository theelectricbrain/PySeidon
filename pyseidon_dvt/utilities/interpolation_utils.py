#!/usr/bin/python2.7
# encoding: utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
from matplotlib.path import Path
from scipy.spatial import KDTree
import scipy.interpolate as interpolate

def closest_point(pt_lon, pt_lat, lon, lat, lonc, latc, tri,
                  debug=False):
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
    points = np.array([[pt_lon], [pt_lat]]).T
    point_list = np.array([lonc[:], latc[:]]).T

    #closest_dist = (np.square((point_list[:, 0] - points[:, 0, np.newaxis])) +
    #                np.square((point_list[:, 1] - points[:, 1, np.newaxis])))

    #closest_point_indexes = np.argmin(closest_dist, axis=1)
    
    #Wesley's optimized version of this bottleneck
    point_list0 = point_list[:, 0]
    points0 = points[:, 0, np.newaxis]
    point_list1 = point_list[:, 1]
    points1 = points[:, 1, np.newaxis]
    
    closest_dist = ((point_list0 - points0) *
                    (point_list0 - points0) +
                    (point_list1 - points1) *
                    (point_list1 - points1)
                    )

    #Check if it is the right index
    if tri[:].max() == (lonc[:].shape[0]-1):
        lo = lonc[:]
        la = latc[:]
    else:
        lo = lon[:]
        la = lat[:]

    index = np.argmin(closest_dist, axis=1)[0]
    triIndex = tri[index]
    triIndex.sort()#due to new version of netCDF4
    trig = Tri.Triangulation(lo[triIndex], la[triIndex], np.array([[0,1,2]]))
    trif = trig.get_trifinder()
    test = -1 * trif.__call__(pt_lon, pt_lat)
    if test: index = np.nan #freak point

    #Thomas' optimized version of this bottleneck    
    #closest_point_indexes = np.zeros(points.shape[0])
    #for i in range(points.shape[0]):
    #    dist=((point_list-points[i,:])**2).sum(axis=1)    
    #    ndx = d.argsort()
    #    closest_point_indexes[i] = ndx[0]
    #if debug: print 'Closest dist: ', closest_dist
    if debug:
        print 'closest_point_indexes', index
        print '...Passed'

    return index

def closest_points( pt_lon, pt_lat, lon, lat, debug=False):
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

    lonc=lon[:]
    latc=lat[:]
    if not type(pt_lon)==list:
        points = np.array([[pt_lon], [pt_lat]]).T
    else:
        points = np.array([pt_lon, pt_lat]).T
    point_list = np.array([lonc, latc]).T

    #closest_dist = (np.square((point_list[:, 0] - points[:, 0, np.newaxis])) +
    #                np.square((point_list[:, 1] - points[:, 1, np.newaxis])))

    #closest_point_indexes = np.argmin(closest_dist, axis=1)
    
    #Wesley's optimized version of this bottleneck
    point_list0 = point_list[:, 0]
    points0 = points[:, 0, np.newaxis]
    point_list1 = point_list[:, 1]
    points1 = points[:, 1, np.newaxis]
    
    closest_dist = ((point_list0 - points0) *
                    (point_list0 - points0) +
                    (point_list1 - points1) *
                    (point_list1 - points1)
                    )
    closest_point_indexes = np.argmin(closest_dist, axis=1)

    #Thomas' optimized version of this bottleneck    
    #closest_point_indexes = np.zeros(points.shape[0])
    #for i in range(points.shape[0]):
    #    dist=((point_list-points[i,:])**2).sum(axis=1)    
    #    ndx = d.argsort()
    #    closest_point_indexes[i] = ndx[0]
    if debug:
        print 'Closest dist: ', closest_dist


    if debug:
        print 'closest_point_indexes', closest_point_indexes
        print '...Passed'

    return closest_point_indexes

def interpN_at_pt(var, pt_x, pt_y, index, trinodes,
                  aw0, awx, awy, debug=False):
    """
    Interpol node variable any given variables at any give location.
    Inputs:
      - var = variable, numpy array, dim=(node) or (time, node) or (time, level, node)
      - pt_x = x coordinate in m to find
      - pt_y = y coordinate in m to find
      - xc = list of x coordinates of var, numpy array, dim= ele
      - yc = list of y coordinates of var, numpy array, dim= ele
      - trinodes = FVCOM trinodes, numpy array, dim=(3,nele)
      - index = index of the nearest element
      - aw0, awx, awy = grid parameters
    Outputs:
      - varInterp = var interpolate at (pt_lon, pt_lat)
    """
    if debug:
        print 'Interpolating at node...'

    n1 = int(trinodes[index,0])
    n2 = int(trinodes[index,1])
    n3 = int(trinodes[index,2])
    # x0 = pt_x - xc[index]
    # y0 = pt_y - yc[index]
    #due to Mitchell's alternative, conversion made in functionsFvcom.py
    x0 = pt_x
    y0 = pt_y

    if len(var.shape)==1:
        var0 = (aw0[0,index] * var[n1]) \
             + (aw0[1,index] * var[n2]) \
             + (aw0[2,index] * var[n3])
        varX = (awx[0,index] * var[n1]) \
             + (awx[1,index] * var[n2]) \
             + (awx[2,index] * var[n3])
        varY = (awy[0,index] * var[n1]) \
             + (awy[1,index] * var[n2]) \
             + (awy[2,index] * var[n3])
    elif len(var.shape)==2:
        var0 = (aw0[0,index] * var[:,n1]) \
             + (aw0[1,index] * var[:,n2]) \
             + (aw0[2,index] * var[:,n3])
        varX = (awx[0,index] * var[:,n1]) \
             + (awx[1,index] * var[:,n2]) \
             + (awx[2,index] * var[:,n3])
        varY = (awy[0,index] * var[:,n1]) \
             + (awy[1,index] * var[:,n2]) \
             + (awy[2,index] * var[:,n3])
    else:
        var0 = (aw0[0,index] * var[:,:,n1]) \
             + (aw0[1,index] * var[:,:,n2]) \
             + (aw0[2,index] * var[:,:,n3])
        varX = (awx[0,index] * var[:,:,n1]) \
             + (awx[1,index] * var[:,:,n2]) \
             + (awx[2,index] * var[:,:,n3])
        varY = (awy[0,index] * var[:,:,n1]) \
             + (awy[1,index] * var[:,:,n2]) \
             + (awy[2,index] * var[:,:,n3])
    varPt = var0 + (varX * x0) + (varY * y0)

    if debug:
        if len(var.shape)==1:
            zi = varPt
            print 'Varpt: ', zi
        print '...Passed'

    #TR comment: squeeze seems to resolve my problem with pydap
    return varPt.squeeze()

def interpN(var,trinodes,aw0,debug=False):
    """
    Interpol node variable at elements.
    Inputs:
      - var = variable, numpy array, dim=(node) or (time, node) or (time, level, node)
      - trinodes = FVCOM trinodes, numpy array, dim=(3,nele)
      - aw0, awx, awy = grid parameters
    Outputs:
      - varInterp = var interpolated
    """
    if debug:
        print 'Interpolating at nodes...'

    n1 = [int(number) for number in trinodes[:,0]]
    n2 = [int(number) for number in trinodes[:,1]]
    n3 = [int(number) for number in trinodes[:,2]]

    if len(var.shape)==1:
        var0 = (aw0[0,:] * var[n1]) \
             + (aw0[1,:] * var[n2]) \
             + (aw0[2,:] * var[n3])
    elif len(var.shape)==2:
        var0 = (aw0[0,:] * var[:,n1]) \
             + (aw0[1,:] * var[:,n2]) \
             + (aw0[2,:] * var[:,n3])
    else:
        var0 = (aw0[0,:] * var[:,:,n1]) \
             + (aw0[1,:] * var[:,:,n2]) \
             + (aw0[2,:] * var[:,:,n3])
    varPt = var0

    if debug: print '...Passed'

    #TR comment: squeeze seems to resolve my problem with pydap
    return varPt.squeeze()
    
def interpE_at_pt(var, pt_x, pt_y, index, triele,
                  a1u, a2u, debug=False):
    """
    Interpol element variable any given variables at any give location.
    Inputs:
      - var = variable, numpy array, dim=(nele) or (time, nele) or (time, level, nele)
      - pt_x = x coordinate in m to find
      - pt_y = y coordinate in m to find
      - xc = list of x coordinates of var, numpy array, dim= nele
      - yc = list of y coordinates of var, numpy array, dim= nele
      - triele = FVCOM triele, numpy array, dim=(3,nele)
      - index = index of the nearest element
      - a1u, a2u = grid parameters
    Outputs:
      - varInterp = var interpolate at (pt_lon, pt_lat)
    """
    if debug:
        print 'Interpolating at element...'

    n1 = int(triele[index,0])
    n2 = int(triele[index,1])
    n3 = int(triele[index,2])

    #Test for ghost points
    test = [-1, var.shape[-1]]

    #TR quick fix: due to error with pydap.proxy.ArrayProxy
    #              not able to cop with numpy.int
    n1 = int(n1)
    n2 = int(n2)
    n3 = int(n3)

    # x0 = pt_x - xc[index]
    # y0 = pt_y - yc[index]
    #due to Mitchell's alternative, conversion made in functionsFvcom.py
    x0 = pt_x
    y0 = pt_y

    if len(var.shape)==1:
        # Treatment of ghost points
        if n1 in test:
            V1 = 0.0
        else:
            V1 = var[n1]
        if n2 in test:
            V2 = 0.0
        else:
            V2 = var[n2]
        if n3 in test:
            V3 = 0.0
        else:
            V3 = var[n3]

        dvardx = (a1u[0,index] * var[index]) \
               + (a1u[1,index] * V1) \
               + (a1u[2,index] * V2) \
               + (a1u[3,index] * V3)
        dvardy = (a2u[0,index] * var[index]) \
               + (a2u[1,index] * V1) \
               + (a2u[2,index] * V2) \
               + (a2u[3,index] * V3)
        varPt = var[index] + (dvardx * x0) + (dvardy * y0)
    elif len(var.shape)==2:
        # Treatment of ghost points
        if n1 in test:
            V1 = np.zeros(var.shape[0])
        else:
            V1 = var[:,n1]
        if n2 in test:
            V2 = np.zeros(var.shape[0])
        else:
            V2 = var[:,n2]
        if n3 in test:
            V3 = np.zeros(var.shape[0])
        else:
            V3 = var[:,n3]

        dvardx = (a1u[0,index] * var[:,index]) \
               + (a1u[1,index] * V1) \
               + (a1u[2,index] * V2) \
               + (a1u[3,index] * V3)
        dvardy = (a2u[0,index] * var[:,index]) \
               + (a2u[1,index] * V1) \
               + (a2u[2,index] * V2) \
               + (a2u[3,index] * V3)
        varPt = var[:,index] + (dvardx * x0) + (dvardy * y0)
    else:
                # Treatment of ghost points
        if n1 in test:
            V1 = np.zeros((var.shape[0], var.shape[1]))
        else:
            V1 = var[:,:,n1]
        if n2 in test:
            V2 = np.zeros((var.shape[0], var.shape[1]))
        else:
            V2 = var[:,:,n2]
        if n3 in test:
            V3 = np.zeros((var.shape[0], var.shape[1]))
        else:
            V3 = var[:,:,n3]

        dvardx = (a1u[0,index] * var[:,:,index]) \
               + (a1u[1,index] * V1) \
               + (a1u[2,index] * V2) \
               + (a1u[3,index] * V3)
        dvardy = (a2u[0,index] * var[:,:,index]) \
               + (a2u[1,index] * V1) \
               + (a2u[2,index] * V2) \
               + (a2u[3,index] * V3)
        varPt = var[:,:,index] + (dvardx * x0) + (dvardy * y0)

    if debug:
        if len(var.shape)==1:
            zi = varPt
            print 'Varpt: ', zi
        print '...Passed'
    #TR comment: squeeze seems to resolve my problem with pydap
    return varPt.squeeze()

def interpE(var, xc, yc, triele,
            a1u, a2u, debug=False):
    """
    Interpol element variable at node locations.
    Inputs:
      - var = variable, numpy array, dim=(nele) or (time, nele) or (time, level, nele)
      - xc = list of x coordinates of var, numpy array, dim= nele
      - yc = list of y coordinates of var, numpy array, dim= nele
      - triele = FVCOM triele, numpy array, dim=(3,nele)
      - a1u, a2u = grid parameters
    Outputs:
      - varInterp = var interpolate at (pt_lon, pt_lat)
    """
    if debug:
        print 'Interpolating at element...'

    n1 = [int(number) for number in triele[:,0]]
    n2 = [int(number) for number in triele[:,1]]
    n3 = [int(number) for number in triele[:,2]]

    #TR quick fix: due to error with pydap.proxy.ArrayProxy
    #              not able to cop with numpy.int

    x0 = xc[:]
    y0 = yc[:]

    if len(var.shape)==1:
        # Treatment of ghost points
        Var = np.hstack((var,0))
        V1 = Var[n1]
        V2 = Var[n2]
        V3 = Var[n3]

        dvardx = (a1u[0,:] * var[:]) \
               + (a1u[1,:] * V1) \
               + (a1u[2,:] * V2) \
               + (a1u[3,:] * V3)
        dvardy = (a2u[0,:] * var[:]) \
               + (a2u[1,:] * V1) \
               + (a2u[2,:] * V2) \
               + (a2u[3,:] * V3)
        varPt = var[:] + (dvardx * x0) + (dvardy * y0)
    elif len(var.shape)==2:
        # Treatment of ghost points
        Var = np.zeros((var.shape[0], var.shape[1]+1))
        Var[:,:-1] = var[:]
        V1 = Var[:,n1]
        V2 = Var[:,n2]
        V3 = Var[:,n3]

        dvardx = (a1u[0,:] * var[:,:]) \
               + (a1u[1,:] * V1) \
               + (a1u[2,:] * V2) \
               + (a1u[3,:] * V3)
        dvardy = (a2u[0,:] * var[:,:]) \
               + (a2u[1,:] * V1) \
               + (a2u[2,:] * V2) \
               + (a2u[3,:] * V3)
        varPt = var[:,:] + (dvardx * x0) + (dvardy * y0)
    else:
        # Treatment of ghost points
        Var = np.zeros((var.shape[0], var.shape[1], var.shape[2]+1))
        Var[:,:,:-1] = var[:]
        V1 = Var[:,:,n1]
        V2 = Var[:,:,n2]
        V3 = Var[:,:,n3]

        dvardx = (a1u[0,:] * var[:,:,:]) \
               + (a1u[1,:] * V1) \
               + (a1u[2,:] * V2) \
               + (a1u[3,:] * V3)
        dvardy = (a2u[0,:] * var[:,:,:]) \
               + (a2u[1,:] * V1) \
               + (a2u[2,:] * V2) \
               + (a2u[3,:] * V3)
        varPt = var[:,:,:] + (dvardx * x0) + (dvardy * y0)

    if debug:
        if len(var.shape)==1:
            zi = varPt
            print 'Varpt: ', zi
        print '...Passed'
    #TR comment: squeeze seems to resolve my problem with pydap
    return varPt.squeeze()

def interp_at_point(var, pt_lon, pt_lat, lon, lat,
                    index, trinodes, tri=[], debug=False):
    """
    Interpol any given variables at any give location.

    Inputs:
      - var = variable, numpy array, dim=(time, nele or node)
      - pt_lon = longitude in degrees to find
      - pt_lat = latitude in degrees to find
      - lon = list of longitudes of var, numpy array, dim=(nele or node)
      - lat = list of latitudes of var, numpy array, dim=(nele or node)
      - trinodes = FVCOM trinodes, numpy array, dim=(3,nele)
    Keywords:
      - tri = triangulation object
    Outputs:
      - varInterp = var interpolate at (pt_lon, pt_lat)
    """
    if debug:
        print 'Interpolating at point...'
    #Finding the right indexes

    #Triangulation
    #if debug:
    #    print triIndex, lon[triIndex], lat[triIndex]
    triIndex = trinodes[index]
    triIndex.sort()#due to new version of netCDF4
    if tri==[]:
        tri = Tri.Triangulation(lon[triIndex], lat[triIndex], np.array([[0,1,2]]))

    trif = tri.get_trifinder()
    trif.__call__(pt_lon, pt_lat)
    if debug:
        if len(var.shape)==1:
            averEl = var[triIndex]
            print 'Var', averEl
            inter = Tri.LinearTriInterpolator(tri, averEl)
            zi = inter(pt_lon, pt_lat)
            print 'zi', zi

    #Choose the right interpolation depending on the variable
    if len(var.shape)==1:
        triVar = np.zeros(triIndex.shape)
        triVar[:] = var[triIndex]
        inter = Tri.LinearTriInterpolator(tri, triVar[:])
        varInterp = inter(pt_lon, pt_lat)      
    elif len(var.shape)==2:
        triVar = np.zeros((var.shape[0], triIndex.shape[0]))
        triVar[:] = var[:, triIndex]
        varInterp = np.ones(triVar.shape[0])
        for i in range(triVar.shape[0]):
            inter = Tri.LinearTriInterpolator(tri, triVar[i,:])
            varInterp[i] = inter(pt_lon, pt_lat)       
    else:
        triVar = np.zeros((var.shape[0], var.shape[1], triIndex.shape[0]))
        triVar[:] = var[:, :, triIndex]
        varInterp = np.ones(triVar.shape[:-1])
        for i in range(triVar.shape[0]):
           for j in range(triVar.shape[1]):
               inter = Tri.LinearTriInterpolator(tri, triVar[i,j,:])
               varInterp[i,j] = inter(pt_lon, pt_lat)

    if debug:
        print '...Passed'

    #TR comment: squeeze seems to resolve my problem with pydap
    return varInterp.squeeze()

def interpE_at_point_bis(var, pt_lon, pt_lat, lonc, latc, debug=False):
    """
    Interpol at awkward locations.

    Inputs:
      - var = variable, numpy array, dim=(time, nele)
      - pt_lon = longitude to find
      - pt_lat = latitude to find
      - lonc = list of longitudes of var, numpy array, dim=(nele)
      - lonc = list of latitudes of var, numpy array, dim=(nele)

    Outputs:
      - varInterp = var interpolate at (pt_lon, pt_lat)
    """
    if debug:
        print 'Interpolating at awkward point...'
    #Finding the right indexes
    points = np.array([[pt_lon], [pt_lat]]).T
    point_list = np.array([lonc[:], latc[:]]).T
    point_list0 = point_list[:, 0]
    points0 = points[:, 0, np.newaxis]
    point_list1 = point_list[:, 1]
    points1 = points[:, 1, np.newaxis]
    
    closest_dist = ((point_list0 - points0) *
                    (point_list0 - points0) +
                    (point_list1 - points1) *
                    (point_list1 - points1)
                    )
    #Finding closest elements
    triIndex = [0,0,0]
    triIndex[0] = np.argmin(closest_dist, axis=1)[0]
    closest_dist[:,triIndex[0]]=np.inf
    triIndex[1] = np.argmin(closest_dist, axis=1)[0]
    closest_dist[:,triIndex[1]]=np.inf
    #test if inside triangle
    test=1
    while test:
        trig = Tri.Triangulation(lonc[triIndex], latc[triIndex], np.array([[0,1,2]]))
        triIndex[2] = np.argmin(closest_dist, axis=1)[0]
        trif = trig.get_trifinder()
        test = -1 * trif.__call__(pt_lon, pt_lat)
        if test: closest_dist[:,triIndex[2]]=np.inf
    #new scheme
    #linear equation based on plane equation
    x1 = lonc[triIndex[0]]; x2 = lonc[triIndex[1]]; x3 = lonc[triIndex[2]]
    y1 = latc[triIndex[0]]; y2 = latc[triIndex[1]]; y3 = latc[triIndex[2]]
    a = np.array([[x2-x1,x3-x1],[y2-y1,y3-y1]])
    b = np.array([pt_lon-x1,pt_lat-y1])
    A, B = np.linalg.solve(a, b)
    #applying weights
    if len(var.shape)==1:
        z1 = var[triIndex[0]]
        z2 = var[triIndex[1]]
        z3 = var[triIndex[2]]
    elif len(var.shape)==2:
        z1 = var[:,triIndex[0]]
        z2 = var[:,triIndex[1]]
        z3 = var[:,triIndex[2]]
    else:
        z1 = var[:,:,triIndex[0]]
        z2 = var[:,:,triIndex[1]]
        z3 = var[:,:,triIndex[2]]
    varInterp = z1 + A * (z2 - z1) + B * (z3 - z1) 
    #end new scheme
    if debug:
        print '...Passed'

    #TR comment: squeeze seems to resolve my problem with pydap
    return varInterp.squeeze()

def interp_linear_to_nodes(var, xc, yc, x, y):
    """Linear interpolation from elements to nodes"""
    L = xc.shape[0]
    M = x.shape[0]
    orig = np.zeros((L, 2))
    ask = np.zeros((M, 2))
    orig[:, 0] = xc
    orig[:, 1] = yc
    ask[:, 0] = x
    ask[:, 1] = y
    interpol = interpolate.LinearNDInterpolator(orig, var)
    varinterp = interpol(ask)
    varinterp[np.where(varinterp==np.nan)]=0.0

    return  varinterp