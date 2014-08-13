from __future__ import division
import numpy as np
import netCDF4 as nc
import sys
import os
import fnmatch
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr
from shortest_element_path import shortest_element_path
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn
#import scipy.io as sio


class FVCOM:
    '''
    A class for FVCOM data.
    As of right now, only takes a filename as input. It will then load in the
    data (except for timeseries, since loading in the whole time series can be
    too large)

    ax can be defined as a region, i.e. a bounding box.
    An example:
        ax = [min(lon_coord), max(lon_coord), min(lat_coord), max(lat_coord)]
    '''

    def __init__(self, filename, elements=slice(None), ax=[], onlyGrid=False, debug=False):
        self.QC = ['raw data']

        if ax:
            self.ax = ax
        else:
            #self.ax = [min(self.lon), max(self.lon), min(self.lat), max(self.lat)]
            self.ax = elements

        self.debug = debug

        if onlyGrid:
            self.loadGrid(filename)
        else:
            self.isMulti(filename, self.ax)



    def isMulti(self, filename, ax):
        split = filename.split('/')

        if split[-1]:
            self.multi = False
            self.load(filename, ax)
        else:
            self.multi = True
            self.loadMulti(filename, ax)

    def loadMulti(self, filename, ax):
        self.matches = self.findFiles(filename)

        #self.loadGrid(filename)

        self.x = np.array([])
        self.y = np.array([])
        self.xc = np.array([])
        self.yc = np.array([])
        self.lonc = np.array([])
        self.latc = np.array([])
        self.lon = np.array([])
        self.lat = np.array([])
        self.siglay = np.array([])
        self.siglev = np.array([])
        self.h = np.array([])
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
            xc = data.variables['xc'][:]
            yc = data.variables['yc'][:]
            lon = data.variables['lon'][:]
            lat = data.variables['lat'][:]
            lonc = data.variables['lonc'][:]
            latc = data.variables['latc'][:]
            siglay = data.variables['siglay'][:]
            siglev = data.variables['siglev'][:]
            h = data.variables['h'][:]
            time = data.variables['time'][:]

            # WES_COMMENT: I do this since el_region takes the class variable
            # self.latc, etc. Therefore this could be done more effectively,
            # but it shouldn't hurt anything since all of these variables are
            # small and quick
            self.lon = data.variables['lon'][:]
            self.lat = data.variables['lat'][:]
            self.lonc = data.variables['lonc'][:]
            self.latc = data.variables['latc'][:]

            self.el_region()
            self.node_region()

            try:
                u = data.variables['u'][:, :, self.region_e]
                v = data.variables['v'][:, :, self.region_e]
                ww = data.variables['ww'][:, :, self.region_e]
                self.D3 = True
            except KeyError:
                self.D3 = False

            print self.region_e
            print data.variables['ua'].shape
            ua = data.variables['ua'][:, self.region_e]
            va = data.variables['va'][:, self.region_e]
            elev = data.variables['zeta'][:, self.region_n]

            self.x = np.hstack((self.x, x))
            self.y = np.hstack((self.y, y))
            self.xc = np.hstack((self.xc, xc))
            self.yc = np.hstack((self.yc, yc))
            self.lon = np.hstack((self.lon, lon))
            self.lat = np.hstack((self.lat, lat))
            self.lonc = np.hstack((self.lonc, lonc))
            self.latc = np.hstack((self.latc, latc))
            self.h = np.hstack((self.h, h))
            self.time = np.hstack((self.time, time))

            if i == 0:
                self.siglay = siglay
                self.siglev = siglev
                self.ua = ua
                self.va = va
                self.elev = elev
                if self.D3:
                    self.u = u
                    self.v = v
                    self.ww = ww

            else:
                self.siglay = np.vstack((self.siglay, siglay))
                self.siglev = np.vstack((self.siglev, siglev))
                self.ua = np.vstack((self.ua, ua))
                self.va = np.vstack((self.va, va))
                self.elev = np.vstack((self.elev, elev))
                if self.D3:
                    self.u = np.vstack((self.u, u))
                    self.v = np.vstack((self.v, v))
                    self.ww = np.vstack((self.ww, ww))


    def findFiles(self, filename):
        '''
        Wesley comment: the name needs to be a linux expression to find files
        you want. For multiple station files, this would work
        name = '*station*.nc'

        For just dngrid_0001 and no restart files:
        name = 'dngrid_0*.nc'
        will work
        '''

        # WES_COMMENT: This has been hardcoded, and once we have a regular
        # naming convention a hard-coded name will work fine.

        name = 'dngrid_0*.nc'
        name = 'small*.nc'
        self.matches = []
        for root, dirnames, filenames in os.walk(filename):
            for filename in fnmatch.filter(filenames, name):
                self.matches.append(os.path.join(root, filename))

        return sorted(self.matches)



    def el_region(self):

        self.region_e = np.argwhere((self.lonc >= self.ax[0]) &
                                    (self.lonc <= self.ax[1]) &
                                    (self.latc >= self.ax[2]) &
                                    (self.latc <= self.ax[3]))

        self.region_e = self.region_e.flatten()
        self.QC.append('Made region for elements out of {}'.format(self.ax))

        if self.debug:
            print self.region_e

    def node_region(self):

        self.region_n = np.argwhere((self.lon >= self.ax[0]) &
                                    (self.lon <= self.ax[1]) &
                                    (self.lat >= self.ax[2]) &
                                    (self.lat <= self.ax[3]))

        self.region_n = self.region_n.flatten()
        self.QC.append('Made region for nodes out of {}'.format(self.ax))

        if self.debug:
            print self.region_n


    def loadGrid(self, filename):
        self.data = nc.Dataset(filename, 'r')
        self.x = self.data.variables['x'][:]
        self.y = self.data.variables['y'][:]
        self.xc = self.data.variables['xc'][:]
        self.yc = self.data.variables['yc'][:]
        self.lon = self.data.variables['lon'][:]
        self.lat = self.data.variables['lat'][:]
        self.lonc = self.data.variables['lonc'][:]
        self.latc = self.data.variables['latc'][:]
        self.nbe = self.data.variables['nbe'][:]

        self.nv = self.data.variables['nv'][:]
        # Make Trinode available for Python indexing
        self.trinodes = self.nv.T - 1

        # get time and adjust it to matlab datenum
        self.julianTime = self.data.variables['time'][:]
        self.time = self.julianTime + 678942
        self.QC.append('Changed time from Julian to matlab datenum')


    def load(self, filename, ax):

#        self.data = nc.Dataset(filename, 'r')
#        # self.data = sio.netcdf.netcdf_file(filename, 'r')
#        self.x = self.data.variables['x'][:]
#        self.y = self.data.variables['y'][:]
#        self.xc = self.data.variables['xc'][:]
#        self.yc = self.data.variables['yc'][:]
#        self.lon = self.data.variables['lon'][:]
#        self.lat = self.data.variables['lat'][:]
#        self.lonc = self.data.variables['lonc'][:]
#        self.latc = self.data.variables['latc'][:]
        self.loadGrid(filename)
        self.h = self.data.variables['h'][:]
        self.nbe = self.data.variables['nbe'][:]
        self.a1u = self.data.variables['a1u'][:]
        self.a2u = self.data.variables['a2u'][:]
        self.aw0 = self.data.variables['aw0'][:]
        self.awx = self.data.variables['awx'][:]
        self.awy = self.data.variables['awy'][:]
        self.siglay = self.data.variables['siglay'][:]
        self.siglev = self.data.variables['siglev'][:]

        self.nv = self.data.variables['nv'][:]
        # Make Trinode available for Python indexing
        self.trinodes = self.nv.T - 1

        # get time and adjust it to matlab datenum
        self.julianTime = self.data.variables['time'][:]
        self.time = self.julianTime + 678942
        self.QC.append('Changed time from Julian to matlab datenum')

        # Need to use len to get size of dimensions
        self.nele = self.data.dimensions['nele']
        self.node = self.data.dimensions['node']

        # Get regions
        if len(ax) == 4:
            self.el_region()
            self.node_region()
        else:
            print ax
            #print ax.shape
            self.region_e = self.closest_point(ax[0], ax[1])
            self.region_n = self.closest_point(ax[0], ax[1], center=False)
            print self.region_e, self.region_n

        # elev timeseries
        print self.data.variables['zeta'].shape
        self.elev = self.data.variables['zeta'][:, self.region_n]

        try:
            self.ww = self.data.variables['ww'][:, :, self.region_e]
            self.u = self.data.variables['u'][:, :, self.region_e]
            self.v = self.data.variables['v'][:, :, self.region_e]
            self.ua = self.data.variables['ua'][:, self.region_e]
            self.va = self.data.variables['va'][:, self.region_e]
            self.D3 = True

        except KeyError:
            self.ua = self.data.variables['ua'][:, self.region_e]
            self.va = self.data.variables['va'][:, self.region_e]
            self.D3 = False

    def centers(self, elements):
        '''Currently doesn't work with whole grid'''

        size = self.trinodes.T[elements].shape[0]
        size1 = self.elev.shape[0]
        elc = np.zeros((size1, size))
        hc = np.zeros((size))
        for ind, value in enumerate(self.trinodes.T[elements]):
            elc[:, ind] = np.mean(self.elev[:, value-1], axis=1)
            hc[ind] = np.mean(self.h[value-1])

        return elc, hc

    def closest_point(self, pt_lon, pt_lat, center=True):
        # def closest_point(self, pt_lon, pt_lat, lon, lat):
        '''
        Finds the closest exact lon, lat to a lon, lat coordinate.
        Example input:
            closest_point([-65.37], [45.34], lon, lat)

        where lon, lat are from data
        '''

        points = np.array([pt_lon, pt_lat]).T

        # point_list = np.array([lon,lat]).T
        if center:
            point_list = np.array([self.lonc, self.latc]).T
        else:
            point_list = np.array([self.lon, self.lat]).T

        closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                        (point_list[:, 1] - points[:, 1, None])**2)

        closest_point_indexes = np.argmin(closest_dist, axis=1)

        return closest_point_indexes

    def harmonics(self, ind, twodim=True, **kwarg):

        if twodim:
            self.coef = ut_solv(self.time, self.ua[:, ind], self.va[:, ind],
                                self.lat[ind], **kwarg)

            self.QC.append('ut_solv done for velocity')

        else:
            self.coef = ut_solv(self.time, self.ua[:, ind], [],
                                self.lat[ind], **kwarg)

            self.QC.append('ut_solv done for elevation')

    def reconstr(self, time):
        if self.coef['aux']['opt']['twodim']:
            self.U, self.V = ut_reconstr(time, self.coef)
            self.QC.append('ut_reconstr done for velocity')
        else:
            self.ts_recon, _ = ut_reconstr(time, self.coef)
            self.QC.append('ut_reconstr done for elevation')

    def graphGrid(self):
        nv = self.nv.T - 1
        h = self.h
        tri = Tri.Triangulation(self.lon, self.lat, triangles=nv.T-1)

        levels = np.arange(-38, -4, 1)   # depth contours to plot

        fig = plt.figure(figsize=(18, 10))
        plt.rc('font', size='22')
        ax = fig.add_subplot(111, aspect=(1.0/np.cos(np.mean(self.lat)*np.pi/180.0)))
        plt.tricontourf(tri, -h, levels=levels, shading='faceted', cmap=plt.cm.gist_earth)
        plt.triplot(tri)
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.gca().patch.set_facecolor('0.5')
        cbar = plt.colorbar()
        cbar.set_label('Water Depth (m)', rotation=-90, labelpad=30)

        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        ax.xaxis.set_major_formatter(ticks)
        ax.yaxis.set_major_formatter(ticks)
        plt.grid()
        plt.show()


if __name__ == '__main__':

    filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
    #filename = '/home/wesley/ncfiles/'

    ind = [-66.3419, -66.3324, 44.2755, 44.2815]

    test = FVCOM(filename, ax=ind)
    test.harmonics(0, cnstit='auto', notrend=True, nodiagn=True)
    #test.reconstr(test.time)

    #test.closest_point([-66.3385], [44.277])
    #t = shortest_element_path(test.latc,test.lonc,test.lat,test.lon,test.nv,test.h)
    #elements, _ = t.getTargets([[41420,39763],[48484,53441],
    #                            [27241,24226],[21706,17458]])



    # t.graphGrid()
