from __future__ import division
import numpy as np
import netCDF4 as nc
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr
from shortest_element_path import shortest_element_path
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import seaborn
import scipy.io as sio


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

    def __init__(self, filename, ax=[]):
        self.QC = ['raw data']

        self.load(filename)

        if ax:
            self.ax = ax
        else:
            self.ax = [min(self.lon), max(self.lon), min(self.lat), max(self.lat)]


    def el_region(self):

        self.region_e = np.argwhere((self.lonc >= self.ax[0]) &
                                    (self.lonc <= self.ax[1]) &
                                    (self.latc >= self.ax[2]) &
                                    (self.latc <= self.ax[3]))

        print self.region_e

    def node_region(self):

        self.region_n = np.argwhere((self.lon >= self.ax[0]) &
                                    (self.lon <= self.ax[1]) &
                                    (self.lat >= self.ax[2]) &
                                    (self.lat <= self.ax[3]))

        print self.region_e


    def load(self, filename):
        #self.data = nc.Dataset(filename, 'r')
        self.data = sio.netcdf.netcdf_file(filename, 'r')
        self.x = self.data.variables['x'][:]
        self.y = self.data.variables['y'][:]
        self.xc = self.data.variables['xc'][:]
        self.yc = self.data.variables['yc'][:]
        self.lon = self.data.variables['lon'][:]
        self.lat = self.data.variables['lat'][:]
        self.lonc = self.data.variables['lonc'][:]
        self.latc = self.data.variables['latc'][:]
        self.nv = self.data.variables['nv'][:]
        self.h = self.data.variables['h'][:]
        self.nbe = self.data.variables['nbe'][:]
        self.a1u = self.data.variables['a1u'][:]
        self.a2u = self.data.variables['a2u'][:]
        self.aw0 = self.data.variables['aw0'][:]
        self.awx = self.data.variables['awx'][:]
        self.awy = self.data.variables['awy'][:]
        self.trinodes = self.data.variables['nv'][:]
        self.siglay = self.data.variables['siglay'][:]
        self.siglev = self.data.variables['siglev'][:]

        # Need to use len to get size of dimensions
        self.nele = self.data.dimensions['nele']
        self.node = self.data.dimensions['node']

        # elev timeseries
        self.el = self.data.variables['zeta']

        # get time and adjust it to matlab datenum
        self.julianTime = self.data.variables['time'][:]
        self.time = self.julianTime + 678942
        self.QC.append('Changed time from Julian to matlab datenum')

        try:
            self.ww = self.data.variables['ww']
            self.u = self.data.variables['u']
            self.v = self.data.variables['v']
            self.ua = self.data.variables['ua']
            self.va = self.data.variables['va']
            self.D3 = True

        except KeyError:
            self.ua = self.data.variables['ua']
            self.va = self.data.variables['va']
            self.D3 = False

    def centers(self):
        size = self.trinodes.T.shape[0]
        size1 = self.el.shape[0]
        elc = np.zeros((size1, size))
        hc = np.zeros((size))
        for i,v in enumerate(self.trinodes.T):
            elc[:, i] = np.mean(self.el[:, v], axis=1)
            hc[i] = np.mean(self.h[v], axis=1)

        return elc, hc

    def hc(self):
        size = self.trinodes.T.shape[0]
        size1 = self.el.shape[0]
        elc = np.zeros((size1, size))
        for i,v in enumerate(self.trinodes.T):
            elc[:, i] = np.mean(self.el[:, v], axis=1)

    def closest_point(self, pt_lon, pt_lat):
    # def closest_point(self, pt_lon, pt_lat, lon, lat):
        '''
        Finds the closest exact lon, lat to a lon, lat coordinate.
        Example input:
            closest_point([-65.37], [45.34], lon, lat)

        where lon, lat are from data
        '''

        points = np.array([pt_lon, pt_lat]).T

        #point_list = np.array([lon,lat]).T
        point_list = np.array([self.lonc, self.latc]).T

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
        nv = self.nv.T -1
        h = self.h
        tri = Tri.Triangulation(self.lon, self.lat, triangles=nv)

        levels=np.arange(-38,-4,1)   # depth contours to plot

        fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')
        ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(self.lat)*np.pi/180.0)))
        plt.tricontourf(tri, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
        plt.triplot(tri)
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.gca().patch.set_facecolor('0.5')
        cbar=plt.colorbar()
        cbar.set_label('Water Depth (m)', rotation=-90,labelpad=30)

        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        ax.xaxis.set_major_formatter(ticks)
        ax.yaxis.set_major_formatter(ticks)
        plt.grid()
        plt.show()






if __name__ == '__main__':

    filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
    test = FVCOM(filename)
    test.harmonics(0, cnstit='auto', notrend=True, nodiagn=True)
    test.reconstr(test.time)

    t = shortest_element_path(test.latc,test.lonc,test.lat,test.lon,test.nv,test.h)
    elements, _ = t.getTargets([[41420,39763],[48484,53441],
                                [27241,24226],[21706,17458]])



    # t.graphGrid()
