from __future__ import division
import numpy as np
from bisect import bisect_left, bisect_right
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import netCDF4 as nc
from createNC import createNC



#def regioner(region, data, name=None, savedir=None, dim='2D'):
#def regioner(filename, data, name=None, savedir=None, dim='2D'):


def loadGrid(filename):
    data = nc.Dataset(filename, 'r')
    x = data.variables['x'][:]
    y = data.variables['y'][:]
    xc = data.variables['xc'][:]
    yc = data.variables['yc'][:]
    lon = data.variables['lon'][:]
    lat = data.variables['lat'][:]
    lonc = data.variables['lonc'][:]
    latc = data.variables['latc'][:]
    a1u = data.variables['a1u'][:]
    a2u = data.variables['a2u'][:]
    aw0 = data.variables['aw0'][:]
    awx = data.variables['awx'][:]
    awy = data.variables['awy'][:]


    nbe = data.variables['nbe'][:].T
    nv = data.variables['nv'][:].T -1

    # Make Trinode available for Python indexing
    #trinodes = nv.T - 1

    return lon, lat, nbe, nv, a1u, a2u, aw0, awx, awy, x, xc, y, yc, lonc, latc


def node_region(ax, lon, lat):

    region_n = np.argwhere((lon >= ax[0]) &
                            (lon <= ax[1]) &
                            (lat >= ax[2]) &
                            (lat <= ax[3]))

    region_n = region_n.flatten()

    return region_n


def regioner(filename, ax):
    """
    Takes as input a region (given by a four elemenTakes as input a region
    (given by a four element NumPy array),
    and some standard data output by ncdatasort and loadnc2d_python
    and returns only the data that lies within the region specified
    in the region arrayt NumPy array),
    and some standard data output by ncdatasort and loadnc2d_python
    and returns only the data that lies within the region specified
    in the region array

    :Parameters:
        **region** -- four element array containing the four corners of the
        region box.  Entires should be in the following form:
            [long1, long2, lat1, lat2] with the following property:
                abs(long1) < abs(long2), etc.

        **data** -- standard python data dictionary for these files

        **name** -- what should  the region be called in the output file?

        **savedir** -- where should the resultant data be saved? Default is
        none, i.e. the data will not be saved, only returned.

        **dim = {'2D', '3D'}** the dimension of the data to use regioner
        on.  Default is 2D.
        """

    lon, lat, nbe, nv, a1u, a2u, aw0, awx, awy, x, xc, y, yc, lonc, latc = loadGrid(filename)

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
    for i in xrange(len(idx)):
        i1 = bisect_left(nv_sortd, idx[i])
        i2 = bisect_right(nv_sortd, idx[i])
        inds = nv_sortedind[i1:i2]
        element_index_tmp[inds % l] = 1
        element_index = np.where(element_index_tmp == 1)[0]
        node_index = np.unique(nv[element_index,:])
        #create new linkage arrays
        nv_tmp = nv[element_index,:]
        L = len(nv_tmp[:,0])
        #nv_tmp2 = np.empty((1, L*3.0))
        nv_tmp2 = np.empty((1, L*3))

    #make a new array of the node labellings for the tri's in the region

    nv2 = nv_tmp.reshape(L * 3, order='F')
    nv2_sortedind = nv2.argsort()
    nv2_sortd = nv2[nv2_sortedind]

    for i in xrange(len(node_index)):
        i1 = bisect_left(nv2_sortd, node_index[i])
        i2 = bisect_right(nv2_sortd, node_index[i])
        inds = nv2_sortedind[i1:i2]
        nv_tmp2[0, inds] = i

    nv_new = np.reshape(nv_tmp2, (L, 3), 'F')
    #now do the same for nbe
    nbe_index = np.unique(nbe[element_index, :])
    nbe_tmp = nbe[element_index,:]
    lnbe = len(nbe_tmp[:,0])
    nbe_tmp2 = np.empty((1, lnbe*3))

    nbe2 = nbe_tmp.reshape(lnbe*3, order='F')
    nbe_sortedind = nbe2.argsort()
    nbe_sortd = nbe2[nbe_sortedind]

    for i in xrange(len(nbe_index)):
        i1 = bisect_left(nbe_sortd, nbe_index[i])
        i2 = bisect_right(nbe_sortd, nbe_index[i])
        inds = nbe_sortedind[i1:i2]
        nbe_tmp2[0, inds] = i

    nbe_new = np.reshape(nbe_tmp2, (lnbe,3), 'F')
    #nbe_new[nbe_new > len(nv_new[:,0]), :] = 0
    nbe_new[nbe_new > len(nv_new[:,0])] = 0

    #create new variables for the region

    data = {}
    data['node_index'] = node_index
    data['element_index'] = element_index
    data['nbe'] = nbe_new
    data['nv'] = nv_new

    data['a1u'] = a1u[:, element_index]
    data['a2u'] = a2u[:, element_index]
    data['aw0'] = aw0[:, element_index]
    data['awx'] = awx[:, element_index]
    data['awy'] = awy[:, element_index]

    data['x'] = x[node_index]
    data['y'] = y[node_index]
    data['xc'] = xc[node_index]
    data['yc'] = yc[node_index]

    data['lon'] = lon[node_index]
    data['lat'] = lat[node_index]
    data['lonc'] = lonc[node_index]
    data['latc'] = latc[node_index]

    data['trigrid'] = Tri.Triangulation(data['lon'], data['lat'], \
                                            data['nv'])


#        data['node_index'] = node_index
#        data['element_index'] = element_index
#        data['nbe'] = nbe_new
#        data['nv'] = nv_new
#        data['a1u'] = a1u[element_index, :]
#        data['a2u'] = a2u[element_index, :]
#        data['h'] = data['h'][node_index]
#        data['uvnodell'] = data['uvnodell'][element_index,:]
#        data['x'] = data['x'][node_index]
#        data['y'] = data['y'][node_index]
#        data['zeta'] = data['zeta'][:,node_index]
#        data['ua'] = data['ua'][:,element_index]
#        data['va'] = data['va'][:,element_index]
#        data['lon'] = data['lon'][node_index]
#        data['lat'] = data['lat'][node_index]
#        data['trigrid'] = mplt.Triangulation(data['lon'], data['lat'], \
#                                             data['nv'])
#        #take care of extra variables if data file is 3D
#        if dim=='3D':
#
#            data['u'] = data['u'][:,:,element_index]
#            data['v'] = data['v'][:,:,element_index]
#            data['ww'] = data['ww'][:,:,element_index]
#            #save the data if that was requested.
#            if savedir != None and name != None:
#                mat_save(data, regionData)
#                return data
    return data

if __name__ == '__main__':

    filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
    ind = [-66.3419, -66.3324, 44.2755, 44.2815]

    print 'Regioning'
    data = regioner(filename, ind)
    print 'Regioned'

    print 'Data loading....'
    ncFile = nc.Dataset(filename, 'r')
    data['h'] = ncFile.variables['h'][data['node_index']]
    data['elev'] = ncFile.variables['zeta'][:, data['node_index']]
    data['ua'] = ncFile.variables['ua'][:, data['element_index']]
    data['va'] = ncFile.variables['va'][:, data['element_index']]
    data['time'] = ncFile.variables['time'][:]
    data['siglay'] = ncFile.variables['siglay'][:, data['node_index']]
    data['siglev'] = ncFile.variables['siglev'][:, data['node_index']]

    print 'Data loaded'

    print 'Creating new NC file for chunked data'
    createNC(data)
