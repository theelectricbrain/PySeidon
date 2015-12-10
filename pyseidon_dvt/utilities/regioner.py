from __future__ import division
import numpy as np
from bisect import bisect_left, bisect_right
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
#quick fix
#import netCDF4 as nc
import scipy.io.netcdf as nc

def node_region(ax, lon, lat):

    region_n = np.argwhere((lon >= ax[0]) &
                            (lon <= ax[1]) &
                            (lat >= ax[2]) &
                            (lat <= ax[3]))

    region_n = region_n.ravel()

    return region_n

def element_region(ax, lonc, latc):

    region_e = np.argwhere((lonc >= ax[0]) &
                            (lonc <= ax[1]) &
                            (latc >= ax[2]) &
                            (latc <= ax[3]))

    region_e = region_e.ravel()

    return region_e


def regioner(gridVar, ax, debug=False):
    """
    Takes as input a region (given by a four elemenTakes as input a region
    (given by a four element NumPy array),
    and some standard data output by ncdatasort and loadnc2d_python
    and returns only the data that lies within the region specified
    in the region arrayt NumPy array),
    and some standard data output by ncdatasort and loadnc2d_python
    and returns only the data that lies within the region specified
    in the region array

    Inputs:
      - region = four element array containing the four corners of the
        region box. Entires should be in the following form:
        [long1, long2, lat1, lat2] with the following property:
        abs(long1) < abs(long2), etc.
      - data = standard python data dictionary for these files
      - name = what should the region be called in the output file
      - savedir = where should the resultant data be saved? Default is
        none, i.e. the data will not be saved, only returned.

    **dim = {'2D', '3D'}** the dimension of the data to use regioner
    on. Default is 2D.
    """
    if debug:
        print 'Reindexing...'
    lon = gridVar.lon[:]
    lat = gridVar.lat[:]
    nbe = gridVar.triele[:]
    nv = gridVar.trinodes[:]
    a1u = gridVar.a1u[:]
    a2u = gridVar.a2u[:]
    aw0 = gridVar.aw0[:]
    awx = gridVar.awx[:]
    awy = gridVar.awy[:]
    x = gridVar.x[:]
    xc = gridVar.xc[:]
    y = gridVar.y[:]
    yc = gridVar.yc[:]
    lonc = gridVar.lonc[:]
    latc = gridVar.latc[:]

    l = nv.shape[0]

    idx = node_region(ax, lon, lat)

    #first, reindex elements in the region
    element_index_tmp = np.zeros(l, int)
    nv_rs = nv.reshape(l*3, order='F')
    #find indices that sort nv_rs
    nv_sortedind = nv_rs.argsort()
    #sort the array
    nv_sortd = nv_rs[nv_sortedind]
    #pick out the values in the region
    if debug:
        print 'Extracting values from box...'
    #TR comment: very slow...gonna need optimisation down the line
    for i in xrange(len(idx)):
        i1 = bisect_left(nv_sortd, idx[i])
        i2 = bisect_right(nv_sortd, idx[i])
        inds = nv_sortedind[i1:i2]
        element_index_tmp[inds % l] = 1
        element_index = np.where(element_index_tmp == 1)[0]
    element_index = element_index.astype(int)

    #TR needs to be inside the loop?
    node_index = np.unique(nv[element_index,:]).astype(int)
    #create new linkage arrays
    nv_tmp = nv[element_index,:]
    L = len(nv_tmp[:,0])
    nv_tmp2 = np.empty((1, L*3))

    #make a new array of the node labellings for the tri's in the region
    if debug:
        print 'Re-labelling nodes...'
    nv2 = nv_tmp.reshape(L * 3, order='F')
    nv2_sortedind = nv2.argsort()
    nv2_sortd = nv2[nv2_sortedind]

    for i in xrange(len(node_index)):
        i1 = bisect_left(nv2_sortd, node_index[i])
        i2 = bisect_right(nv2_sortd, node_index[i])
        inds = nv2_sortedind[i1:i2]
        nv_tmp2[0, inds] = i

    nv_new = np.reshape(nv_tmp2, (L, 3), 'F')

    #now do the same for nbe...sort of
    nbe_index = np.unique(nbe[element_index, :])
    #ghost points
    ghost=np.asarray(list(set(nbe_index) - set(element_index)))
    
    nbe_tmp = nbe[element_index,:]
    lnbe = len(nbe_tmp[:,0])
    #nbe_tmp2 = np.empty((1, lnbe*3))
    #TR: np.empty sometimes generates freak values
    nbe_tmp2 = np.zeros((1, lnbe*3))

    if debug:
        print 'Re-labelling elements...'
    nbe2 = nbe_tmp.reshape(lnbe*3, order='F')
    nbe_sortedind = nbe2.argsort()
    nbe_sortd = nbe2[nbe_sortedind]
    
    #TR: iterator
    I = 0
    for i in xrange(len(nbe_index)):
        #TR: check if ghost point
        if not nbe_index[i] in ghost:
            i1 = bisect_left(nbe_sortd, nbe_index[i])
            i2 = bisect_right(nbe_sortd, nbe_index[i])
            inds = nbe_sortedind[i1:i2]
            nbe_tmp2[0, inds] = I
            I += 1 

    nbe_new = np.reshape(nbe_tmp2, (lnbe,3), 'F')

    #create new variables for the region

    data = {}
    data['node_index'] = node_index
    data['element_index'] = element_index
    data['nbe'] = nbe_new.astype(int)
    data['nv'] = nv_new.astype(int)

    data['a1u'] = a1u[:, element_index]
    data['a2u'] = a2u[:, element_index]
    data['aw0'] = aw0[:, element_index]
    data['awx'] = awx[:, element_index]
    data['awy'] = awy[:, element_index]

    data['x'] = x[node_index]
    data['y'] = y[node_index]
    data['xc'] = xc[element_index]
    data['yc'] = yc[element_index]

    data['lon'] = lon[node_index]
    data['lat'] = lat[node_index]
    data['lonc'] = lonc[element_index]
    data['latc'] = latc[element_index]

    data['triangle'] = Tri.Triangulation(data['lon'], data['lat'], \
                                        data['nv'])

    return data

