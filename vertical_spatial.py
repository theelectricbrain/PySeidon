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


def time_index(a):
    a = a.reindex(index=a.index.to_datetime())
    return a

def magnitude(a):
    mag = np.sqrt(a['u']**2+a['v']**2+a['w']**2)
    return mag

def theta(a):
    a = np.arctan2(a['v'],a['u'])
    return a

def date2py(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + \
        timedelta(days=matlab_datenum%1) - timedelta(days = 366)

    return python_datetime

def main(filename, eastwest=True):
    #filename = '/array2/data3/rkarsten/july_2012/output/dngrid_0001_02.nc'
    #filename = '/EcoII/june_2013/output/dngrid_0001_week2.nc'
    #filename = '/EcoII/july_2012/output/dngrid_0001_03.nc'

    #eastwest = True

    siglay = np.array([0.98999,0.94999,0.86999,0.74999,0.58999,0.41000,0.25000,0.13000,0.05000,0.01000])


    data = FVCOM(filename)

    if eastwest:
        # East- West
        ind = data.closest_point([-66.3419, -66.3324], [44.2778, 44.2778])
    else:
        # North-South
        ind = data.closest_point([-66.3385, -66.3385], [44.2815, 44.2755])

    short_path = shortest_element_path(data.lonc, data.latc,
                                        data.lon, data.lat,
                                        data.nv, data.h)

    #short_path = shortest_element_path(data.xc, data.yc,
    #                                    data.x, data.y,
    #                                    data.nv, data.h)

    el, _ = short_path.getTargets([ind])
    #short_path.graphGrid()
    #plt.show()
    #saveName = './figures/e-wPath.png'
    #plt.savefig(saveName, bbox_inches=0)
    plt.clf()
    print 'Path Saved'

    t_slice = ['2014-02-02T06:45:00','2014-02-02T07:05:00']
    t_slice = np.array(t_slice,dtype='datetime64[us]')

    t = data.time.shape[0]
    l = []
    for i in range(t):
        date = datetime.fromordinal(int(data.time[i]))+timedelta(days=data.time[i]%1)-timedelta(days=366)
        l.append(date)

    time = np.array(l,dtype='datetime64[us]')
    if t_slice.shape[0] != 1:
        argtime = np.argwhere((time>=t_slice[0])&(time<=t_slice[-1])).flatten()

    #vel = np.sqrt(nc['u'][argtime,:,el]**2+nc['v'][argtime,:,el]**2+nc['ww'][argtime,:,el]**2)
    #vel = np.sqrt(data.u[argtime,:,el]**2+data.u[argtime,:,el]**2+data.ww[argtime,:,el]**2)
    print len(el[0])
    print 'Loading timeseries'
    #vel = np.sqrt(data.u[:, :, el[0]]**2 + data.v[:, :, el[0]]**2 + data.ww[:, :, el[0]]**2)
    print 'U'
    u = data.u[:, :, el[0]]
    print 'V'
    v = data.v[:, :, el[0]]
    #ww = data.ww[:, :, el[0]]

    print 'Calculating'
    #vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')
    vel = ne.evaluate('sqrt(u**2 + v**2)')

    mean_vel = np.mean(vel, axis=0)

    lat = data.latc[el]
    lon = data.lonc[el]

    if eastwest:
        line = lon
    else:
        line = lat

    print data.time[0]
    new = date2py(data.time[0])
    print new
    print vel.shape
    print mean_vel.shape


    print 'Calculating elc and hc'
    size = data.trinodes.T[el].shape[0]
    size1 = data.el.shape[0]
    elc = np.zeros((size1, size))
    hc = np.zeros((size))
    for ind,value in enumerate(data.trinodes.T[el[0]]):
        elc[:, ind] = np.mean(data.el[:, value-1], axis=1)
        hc[ind] = np.mean(data.h[value-1])

    print 'Saving mat file'
    mat = {'u':u, 'v':v, 'latc':lat, 'lonc':lon, 'time':data.time,
        'siglay':data.siglay, 'siglev':data.siglev, 'ua':data.ua[:, el[0]],
        'va':data.va[:, el[0]],
        'elc':elc, 'hc':hc}

    if eastwest:
        sio.savemat('east-west.mat', mat)
    else:
        sio.savemat('north-south.mat', mat)

    vmax = 2.5
    vmin = 0

    print 'Plotting'

    fig,ax = plt.subplots()
    plt.rc('font',size='22')
    levels = np.linspace(0,3.3,34)
    cs = ax.contourf(line,siglay,mean_vel,levels=levels, cmap=plt.cm.jet)
    ax.contour(line,siglay,mean_vel,cs.levels, colors='k')
    cbar = fig.colorbar(cs,ax=ax)
    cbar.set_label(r'Velocity $(m/s)$', rotation=-90,labelpad=30)
    if eastwest:
        ax.set_xlabel('Longitude')
    else:
        ax.set_xlabel('Latitude')
    ax.set_title('vel_mean')
    scale = 1
    ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
    ax.xaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_formatter(ticks)

    if eastwest:
        saveName = './figures/eastwest/mean_vel.png'.format(i)
    else:
        saveName = './figures/northsouth/mean_vel.png'.format(i)

    #plt.show()
    plt.savefig(saveName, bbox_inches=0)
    plt.clf()

    for i in range(vel.shape[0]):
        print i
        fig,ax = plt.subplots()
        plt.rc('font',size='12')
        levels = np.linspace(0,3.3,34)
        cs = ax.contourf(line,siglay,vel[i,:],levels=levels, cmap=plt.cm.jet)
        ax.contour(line,siglay,vel[i,:],cs.levels,colors='k',hold='on')
        cbar = fig.colorbar(cs,ax=ax)
        cbar.set_label(r'Velocity $(m/s)$', rotation=-90,labelpad=30)
        title = '{}'.format(date2py(data.time[i]))
        ax.set_title(title)
        if eastwest:
            ax.set_xlabel('Longitude')
        else:
            ax.set_xlabel('Latitude')
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        ax.xaxis.set_major_formatter(ticks)
        ax.yaxis.set_major_formatter(ticks)

        if eastwest:
            saveName = './figures/eastwest/figure{0:>04d}.png'.format(i)
        else:
            saveName = './figures/northsouth/figure{0:>04d}.png'.format(i)

        #plt.show()
        plt.savefig(saveName, bbox_inches=0)
        plt.clf()

if __name__ == '__main__':

    filename = '/EcoII/july_2012/output/dngrid_0001_03.nc'

    main(filename, eastwest=True)
    main(filename, eastwest=False)
