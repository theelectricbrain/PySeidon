from __future__ import division
import numpy as np
import matplotlib
import netCDF4 as net
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import scipy.io as sio



if __name__ == "__main__":

    GPBPb = [-66.3391,44.2761]
    fs5 = [-66.3381,44.2505]
    fileDirll = '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/'
    fileDir = fileDirll
    filename = 'dngrid_0001_week2.nc'
    nc = net.Dataset(fileDirll+filename).variables
    lat = nc['lat'][:]
    lon = nc['lon'][:]
    print max(lat)
    nc = net.Dataset(fileDir+filename).variables
    time = nc['time'][:]+678942
    time = np.array(time)
    t = time.shape[0]
    l = []
    for i in range(t):
        date = datetime.fromordinal(int(time[i]))+timedelta(days=time[i]%1)-timedelta(days=366)
        l.append(date)
    time = np.array(l,dtype='datetime64[us]')

    t_slice = ['2013-06-22T00:00:00','2013-06-24T00:00:00']
    #%t_slice = [time[0],time[10]]
    t_slice = np.array(t_slice,dtype='datetime64[us]')

    if t_slice.shape[0] != 1:
        argtime = np.argwhere((time>=t_slice[0])&(time<=t_slice[-1])).flatten()
    print argtime

    h = nc['h'][:]
    zeta = nc['zeta'][argtime,:]+h[None,:]
    nv = nc['nv'][:].T-1
    siglay = nc['siglay'][:]
    z = zeta[:,None,:]*siglay[None,:,:]
    dep = np.zeros([argtime.shape[0],siglay.shape[0],nv.shape[0]])
    for i in range(z.shape[0]):
        for j in range(nv.shape[0]):
            el = (z[i,:,nv[j,0]]+z[i,:,nv[j,1]]+z[i,:,nv[j,2]])/3
            dep[i,:,j] = el


    u = nc['u'][argtime,:,:]
    v = nc['v'][argtime,:,:]
    w = nc['ww'][argtime,:,:]


    bot_lvl = 4
    top_lvl = 1
    vel = np.sqrt(u[:]**2+v[:]**2+w[:]**2)
    #topl = np.argwhere((vel[:,top_lvl,:] < vel[:,1,:])).flatten()
    #vel[:,:,topl] = 0
    speed1=vel[:,top_lvl,:]
    speed4=vel[:,bot_lvl,:]
    vel = speed4-speed1
    z = dep[:,bot_lvl,:] - dep[:,top_lvl,:]
    dudz = vel/z

    mdict = {'lon':lon, 'lat':lat, 'nv':nv, 'dudz':dudz, 'speed1':speed1, 'speed4':speed4}
    sio.savemat('dudz.mat', mdict)




    levels = np.arange(-40,5,1)
    vmax = 0.08
    vmin = -0.08

    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    for i in range(vel.shape[0]):
        grid =  Tri.Triangulation(lon,lat,triangles=nv)
        fig = plt.figure()
        plt.rc('font',size='22')
        ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
        CS = ax.tricontour(grid, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
        #ax.clabel(CS, fontsize=9, inline=1,colors='k',fmt='%1d')
        f = ax.tripcolor(grid, dudz[i,:],vmax=vmax,vmin=vmin,cmap=plt.cm.PuOr)
        frame = plt.gca().patch.set_facecolor('0.5')
        cbar = fig.colorbar(f,ax=ax)
        cbar.set_label(r'Vertical Shear', rotation=-90,labelpad=30)
        plt.scatter(GPBPb[0],GPBPb[1],s=200,color='black')
        plt.scatter(fs5[0],fs5[1],s=200,color='black')
        #plt.title(str(time[j]))
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.grid()
        scale = 1
        ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        ax.xaxis.set_major_formatter(ticks)
        ax.yaxis.set_major_formatter(ticks)
        ax.set_xlim([-66.35,-66.325])
        ax.set_ylim([44.268,44.282])
        plt.show()


