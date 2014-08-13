from __future__ import division
import numpy as np
import netCDF4 as nc
import os
import fnmatch
# Need to add closest point

class station:
    def __init__(self, filename, elements=slice(None)):

        self.isMulti(filename, elements)

    def isMulti(self, filename, elements):
        split = filename.split('/')
        if split[-1]:
            self.multi = False
            self.load(filename, elements)
        else:
            self.multi = True
            self.loadMulti(filename, elements)

    def load(self, filename, elements):
        self.data = nc.Dataset(filename)
        self.x = self.data.variables['x'][:]
        self.y = self.data.variables['y'][:]
        self.lon = self.data.variables['lon'][:]
        self.lat = self.data.variables['lat'][:]
        self.siglay = self.data.variables['siglay'][:]
        self.siglev = self.data.variables['siglev'][:]
        self.h = self.data.variables['h'][:]
        self.time_JD = self.data.variables['time_JD'][:]
        self.time_second = self.data.variables['time_second'][:]
        self.time = self.time_JD + 678942 + self.time_second / (24*3600)
        self.u = self.data.variables['u'][:, :, elements]
        self.v = self.data.variables['v'][:, :, elements]
        self.ww = self.data.variables['ww'][:, :, elements]
        self.ua = self.data.variables['ua'][:, elements]
        self.va = self.data.variables['va'][:, elements]
        self.elev = self.data.variables['zeta'][:, elements]

    def loadMulti(self, filename, elements):
        self.matches = self.findFiles(filename)

        self.x = np.array([])
        self.y = np.array([])
        self.lon = np.array([])
        self.lat = np.array([])
        self.siglay = np.array([])
        self.siglev = np.array([])
        self.h = np.array([])
        self.time_JD = np.array([])
        self.time_second = np.array([])
        self.time = np.array([])
        self.u = np.array([])
        self.v = np.array([])
        self.ww = np.array([])
        self.ua = np.array([])
        self.va = np.array([])
        self.elev = np.array([])

        for i, v in enumerate(self.matches):
            data = nc.Dataset(v, 'r')
            x = data.variables['x'][:]
            y = data.variables['y'][:]
            lon = data.variables['lon'][:]
            lat = data.variables['lat'][:]
            siglay = data.variables['siglay'][:]
            siglev = data.variables['siglev'][:]
            h = data.variables['h'][:]
            t_JD = data.variables['time_JD'][:]
            t_second = data.variables['time_second'][:]
            u = data.variables['u'][:, :, elements]
            v = data.variables['v'][:, :, elements]
            ww = data.variables['ww'][:, :, elements]
            ua = data.variables['ua'][:, elements]
            va = data.variables['va'][:, elements]
            elev = data.variables['zeta'][:, elements]

            self.x = np.hstack((self.x, x))
            self.y = np.hstack((self.y, y))
            self.lon = np.hstack((self.lon, lon))
            self.lat = np.hstack((self.lat, lat))
            self.h = np.hstack((self.h, h))
            self.time_JD = np.hstack((self.time_JD, t_JD))
            self.time_second = np.hstack((self.time_second, t_second))

            if i == 0:
                self.siglay = siglay
                self.siglev = siglev
                self.u = u
                self.v = v
                self.ww = ww
                self.ua = ua
                self.va = va
                self.elev = elev

            else:
                self.siglay = np.vstack((self.siglay, siglay))
                self.siglev = np.vstack((self.siglev, siglev))
                self.u = np.vstack((self.u, u))
                self.v = np.vstack((self.v, v))
                self.ww = np.vstack((self.ww, ww))
                self.ua = np.vstack((self.ua, ua))
                self.va = np.vstack((self.va, va))
                self.elev = np.vstack((self.elev, elev))


        self.time = self.time_JD + 678942 + self.time_second / (24*3600)

    def findFiles(self, filename):
        '''
        Wesley comment: the name needs to be a linux expression to find files
        you want. For multiple station files, this would work
        name = '*station*.nc'

        For just dngrid_0001 and no restart files:
        name = 'dngrid_0*.nc'
        will work
        '''

        name = '*station*.nc'
        self.matches = []
        for root, dirnames, filenames in os.walk(filename):
            for filename in fnmatch.filter(filenames, name):
                self.matches.append(os.path.join(root, filename))

        return sorted(self.matches)


if __name__ == '__main__':
    #filename = '/array2/data3/rkarsten/dncoarse_3D/output2/dn_coarse_station_timeseries.nc'
    #filename = '/array2/data3/rkarsten/dncoarse_3D/output2/dn_coarse_station_timeseries.nc'
    #filename = '/EcoII/EcoEII_server_data_tree/data/simulated/FVCOM/dngrid/june_2013_3D/'
    multi = True
    if multi:
        #filename = '/home/wesley/ncfiles/'
        filename = '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/'
    else:
        filename = '/home/wesley/ncfiles/dn_coarse_station_timeseries.nc'

    data = station(filename)
