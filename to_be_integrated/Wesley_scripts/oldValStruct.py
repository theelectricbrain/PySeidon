from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import cPickle as pickle
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv
import scipy.io as sio
from stationClass import station
from adcpClass import ADCP
from fvcomClass import FVCOM

def mjd2num(x):

    y = x + 678942

    return y


def closest_point(points, lon, lat):

    point_list = np.array([lon,lat]).T

    closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                    (point_list[:, 1] - points[:, 1, None])**2)

    closest_point_indexes = np.argmin(closest_dist, axis=1)

    return closest_point_indexes


def datetime2matlabdn(dt):
    # ordinal = dt.toordinal()
    mdn = dt + timedelta(days=366)
    frac = (dt-datetime(dt.year, dt.month, dt.day, 0, 0, 0)).seconds / \
        (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac



def main(fvFiles, adcpFiles, tideFiles, debug=False):

    fvdebugData = FVCOM(fvdebug)
    saveName = 'validationStruct.p'
    #Name = 'june_2013_3D_station'
    Struct = {}

    for fvFile in fvFiles:

        print fvFile
        fvData = station(fvFile)

        struct = np.array([])
        for adcpFile in adcpFiles:
            print adcpFile
            adcpData = ADCP(adcpFile)
            lonlat = np.array([adcpData.lon[0], adcpData.lat[0]]).T
            #lonlat = np.array([adcpData.x[0], adcpData.y[0]]).T
            #ind = closest_point(lonlat, fvData.lon, fvData.lat)
            newind = closest_point(lonlat, fvdebugData.lonc, fvdebugData.latc)
            #ind = closest_point(lonlat, fvData.x, fvData.y)
            new = np.array([fvdebugData.xc[newind], fvdebugData.yc[newind]])
            ind = closest_point(new.T, fvData.x, fvData.y)

            print ind
            print adcpData.mtime.shape
            print adcpData.ua.shape
            print adcpData.va.shape
            print adcpData.surf.shape

            adcpVelCoef = ut_solv(adcpData.mtime, adcpData.ua,
                            adcpData.va, adcpData.lat[0],
                            cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, coef_int=True)

            adcpElevCoef = ut_solv(adcpData.mtime, adcpData.surf,
                            [], adcpData.lat[0],
                            cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, coef_int=True)

            adcpName = adcpFile.split('/')[-1].split('.')[0]
            #WB_COMMENT: Doesn't currently work
            obs = pd.DataFrame({'u':adcpData.ua, 'v':adcpData.va, 'elev':adcpData.surf})

            print fvData.time.shape
            print fvData.ua[:, ind].shape
            print fvData.va[:, ind].shape
            print fvData.lat[ind].shape

            fvVelCoef = ut_solv(fvData.time, fvData.ua[:, ind].flatten(),
                                fvData.va[:, ind].flatten(),
                        adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                        method='ols', nodiagn=True, linci=True, conf_int=True)

            print fvData.elev[:, ind].shape
            fvElevCoef = ut_solv(fvData.time, fvData.elev[:, ind].flatten(), [],
                        adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                        method='ols', nodiagn=True, linci=True, conf_int=True)

            mod = pd.DataFrame({'ua':fvData.ua[:, ind].flatten(),
                                'va':fvData.va[:, ind].flatten(),
                                'elev':fvData.elev[:, ind].flatten()})


            obs_loc = {'name': adcpName, 'type':'ADCP', 'lat':fvdebugData.lat[newind],
                    'lon':fvdebugData.lon[newind], 'obs_timeseries':obs,
                    'mod_timeseries':mod, 'obs_time':adcpData.mtime,
                    'mod_time':fvData.time, 'vel_obs_harmonics':adcpVelCoef,
                    'elev_obs_harmonics':adcpElevCoef,
                    'vel_mod_harmonics':fvVelCoef, 'elev_mod_harmonics':fvElevCoef}

            struct = np.hstack((struct, obs_loc))

        Struct[Name] = struct

    if debug:
        pickle.dump(Struct, open("structADCP.p", "wb"))

    pickle.dump(Struct, open(saveName, "wb"))
    return Struct

if __name__ == '__main__':

    fvFiles = ['/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/']

    adcpFiles = ['/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPa_avg5.mat',
                 '/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPb_avg5.mat']

    fvdebug = '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/dngrid_0001_week2.nc'

    tideFiles = \
    ['/EcoII/EcoEII_server_data_tree/data/observed/GP/TideGauge/Westport_015892_20140325_1212_Z.mat',
     '/EcoII/EcoEII_server_data_tree/data/observed/DG/TideGauge/DigbyWharf_015893_20140115_2221_Z.mat']

    Struct = main(fvFiles, adcpFiles, tideFiles)
