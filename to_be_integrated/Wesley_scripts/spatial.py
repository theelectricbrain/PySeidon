from __future__ import division
import numpy as np
from datetime import datetime
from datetime import timedelta
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from shortest_element_path import shortest_element_path
from fvcomClass import FVCOM
import numexpr as ne
import scipy.io as sio
import matplotlib.tri as Tri


filename = '/EcoII/july_2012/output/dngrid_0001_03.nc'
filename = '/EcoII/june_2013/output/dngrid_0001_week2.nc'
filename = '/EcoII/july_2012/output/dngrid_0001_03.nc'
siglay = np.array([0.98999,0.94999,0.86999,0.74999,0.58999,0.41000,0.25000,0.13000,0.05000,0.01000])

ax = [-66.35, -66.325, 44.261, 44.272]
#data = FVCOM(filename)
data = FVCOM(filename, ax=ax)
data.el_region()

print 'Loading Data'
print 'U'
u = data.u[:, :, data.region_e[:,0]]
print 'V'
v = data.v[:, :, data.region_e[:,0]]
#ww = data.ww[:, :, data.region_e[:,0]]

#vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')
vel = ne.evaluate('sqrt(u**2 + v**2)')

mean_vel = np.mean(vel, axis=0)
h = data.h

levels = np.arange(-40,5,1)
vmax = 0.08
vmin = -0.08

print 'Plotting... '
for i in xrange(mean_vel.shape[0]):
    grid = Tri.Triangulation(data.lon, data.lat, triangles=data.trinodes.T-1)
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(data.lat)*np.pi/180.0)))
    CS = ax.tricontour(grid, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
    #ax.clabel(CS, fontsize=9, inline=1,colors='k',fmt='%1d')
    f = ax.tripcolor(grid, mean_vel[i,:],vmax=vmax,vmin=vmin,cmap=plt.cm.PuOr)
    frame = plt.gca().patch.set_facecolor('0.5')
    cbar = fig.colorbar(f,ax=ax)
    cbar.set_label(r'Velocity', rotation=-90,labelpad=30)
    #plt.scatter(GPBPb[0],GPBPb[1],s=200,color='black')
    #plt.scatter(fs5[0],fs5[1],s=200,color='black')
    #plt.title(str(time[j]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.grid()
    scale = 1
    ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
    ax.xaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_formatter(ticks)
    ax.set_xlim([-66.35,-66.325])
    ax.set_ylim([44.261,44.272])
    plt.show()
