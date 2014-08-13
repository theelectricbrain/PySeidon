from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import cPickle as pickle
import sys
import os
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv
import scipy.io as sio
from stationClass import station
from adcpClass import ADCP
from fvcomClass import FVCOM
from tidegaugeClass import Tidegauge

def mjd2num(x):

    y = x + 678942

    return y


def closest_point(points, lon, lat):

    point_list = np.array([lon,lat]).T

    print point_list
    print points
    closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                    (point_list[:, 1] - points[:, 1, None])**2)

    print closest_dist
    closest_point_indexes = np.argmin(closest_dist, axis=1)

    return closest_point_indexes


def datetime2matlabdn(dt):
    # ordinal = dt.toordinal()
    mdn = dt + timedelta(days=366)
    frac = (dt-datetime(dt.year, dt.month, dt.day, 0, 0, 0)).seconds / \
        (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac



def main(fvFiles, adcpFiles, tideFiles, isStation=True, ax=[], debug=False):

    #fvdebugData = FVCOM(fvdebug)
    saveName = 'validationStruct.p'
    #Name = 'june_2013_3D_station'
    #Struct = {}

    struct = np.array([])
    for adcpFile in adcpFiles:
        print adcpFile
        adcpData = ADCP(adcpFile)
        lonlat = np.array([adcpData.lon[0], adcpData.lat[0]]).T

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

        #adcpName = adcpFile.split('/')[-1].split('.')[0]

        adcp_obs = {'ua':adcpData.ua,
                    'va':adcpData.va,
                    'elev':adcpData.surf,
                    'u':adcpData.east_vel,
                    'v':adcpData.north_vel,
                    'bins':adcpData.bins}

#        adcp_obs = pd.DataFrame({'ua':adcpData.ua,
#                                 'va':adcpData.va,
#                                 'elev':adcpData.surf,
#                                 'u':adcpData.east_vel,
#                                 'v':adcpData.north_vel})

        for fvFile in fvFiles:

            print fvFile
            saveName = fvFile + 'validationStruct.p'
            if isStation:
                fvData = station(fvFile)
                ind = closest_point(lonlat, fvData.lon, fvData.lat)
            else:
                #ax = np.array([adcpData.lon[0], adcpData.lat[0]]).T
                ax = [[adcpData.lon[0][0]], [adcpData.lat[0][0]]]
                #ax = [adcpData.lon[0][0], adcpData.lat[0][0]]
                fvData = FVCOM(fvFile, ax)
                #print ax
#                lonlat = np.array([[adcpData.lon[0][0],
#                                   adcpData.lat[0][0]]])
#                ind = closest_point(lonlat, fvData.lon, fvData.lat)
#                print ind

#                ind = fvData.closest_point([adcpData.lon[0][0]],
#                                           [adcpData.lat[0][0]])


            # right one
            #ind = closest_point(lonlat, fvData.lon, fvData.lat)

            #lonlat = np.array([adcpData.x[0], adcpData.y[0]]).T
            #newind = closest_point(lonlat, fvdebugData.lonc, fvdebugData.latc)
            #ind = closest_point(lonlat, fvData.x, fvData.y)
            #new = np.array([fvdebugData.xc[newind], fvdebugData.yc[newind]])
            #ind = closest_point(new.T, fvData.x, fvData.y)

            print fvData.time.shape
            print fvData.ua.shape
            print fvData.ua
            #print fvData.ua[:, ind].shape
            #print fvData.va[:, ind].shape
            #print fvData.lat[ind].shape

            if isStation:
                fvVelCoef = ut_solv(fvData.time, fvData.ua[:, ind].flatten(),
                                    fvData.va[:, ind].flatten(),
                            adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, conf_int=True)

                print fvData.elev[:, ind].shape
                fvElevCoef = ut_solv(fvData.time, fvData.elev[:, ind].flatten(), [],
                            adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, conf_int=True)

                mod = {'ua':fvData.ua[:, ind].flatten(),
                        'va':fvData.va[:, ind].flatten(),
                        'elev':fvData.elev[:, ind].flatten(),
                        'u':fvData.u,
                        'v':fvData.v}
            else:
                fvVelCoef = ut_solv(fvData.time, fvData.ua.flatten(),
                                    fvData.va.flatten(),
                            adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, conf_int=True)

                #print fvData.elev[:, ind].shape
                fvElevCoef = ut_solv(fvData.time, fvData.elev.flatten(), [],
                            adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, conf_int=True)

                if fvData.D3:
                    mod = {'ua':fvData.ua.flatten(),
                            'va':fvData.va.flatten(),
                            'elev':fvData.elev.flatten(),
                            'u':fvData.u,
                            'v':fvData.v}
                else:
                    mod = {'ua':fvData.ua.flatten(),
                            'va':fvData.va.flatten(),
                            'elev':fvData.elev.flatten()}



            obs_loc = {'name': adcpFile,
                       'type':'ADCP',
                       'lat':adcpData.lat[0],
                       'lon':adcpData.lon[0],
                       'obs_timeseries':adcp_obs,
                       'mod_timeseries':mod,
                       'obs_time':adcpData.mtime,
                       'mod_time':fvData.time,
                       'vel_obs_harmonics':adcpVelCoef,
                       'elev_obs_harmonics':adcpElevCoef,
                       'vel_mod_harmonics':fvVelCoef,
                       'elev_mod_harmonics':fvElevCoef}
                       #'adcp_bins':adcpData.bins}

#            obs_loc = {'name': adcpName, 'type':'ADCP', 'lat':fvdebugData.lat[newind],
#                    'lon':fvdebugData.lon[newind], 'obs_timeseries':adcp_obs,
#                    'mod_timeseries':mod, 'obs_time':adcpData.mtime,
#                    'mod_time':fvData.time, 'vel_obs_harmonics':adcpVelCoef,
#                    'elev_obs_harmonics':adcpElevCoef,
#                    'vel_mod_harmonics':fvVelCoef, 'elev_mod_harmonics':fvElevCoef}

            struct = np.hstack((struct, obs_loc))


    for tideFile in tideFiles:

        print tideFile

        tideData = Tidegauge(tideFile)
        ut_constits = ['M2','S2','N2','K2','K1','O1','P1','Q1']
        tideData.harmonics(cnstit=ut_constits, notrend=True,
                           rmin=0.95, method='ols', nodiagn=True, linci=True,
                           ordercnstit='frq')

        tide_obs = {'data':tideData.data, 'elev':tideData.elev}

        for fvFile in fvFiles:

            print fvFile

            if isStation:
                fvData = station(fvFile)
                ind = np.argmin(np.sqrt((fvData.lon-tideData.lon)**2+(fvData.lat-tideData.lat)**2))
                #ind = closest_point(lonlat, fvData.lon, fvData.lat)
            else:
                #ax = np.array([adcpData.lon[0], adcpData.lat[0]]).T
                ax = [[tideData.lon], [tideData.lat]]
                fvData = FVCOM(fvFile, ax)

            if isStation:

                print fvData.elev[:, ind].shape
                fvElevCoef = ut_solv(fvData.time, fvData.elev[:, ind].flatten(), [],
                            adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, conf_int=True)

                mod = {'ua':fvData.ua[:, ind].flatten(),
                        'va':fvData.va[:, ind].flatten(),
                        'elev':fvData.elev[:, ind].flatten(),
                        'u':fvData.u,
                        'v':fvData.v}
            else:

                #print fvData.elev[:, ind].shape
                fvElevCoef = ut_solv(fvData.time, fvData.elev.flatten(), [],
                            adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, conf_int=True)

                if fvData.D3:
                    mod = {'ua':fvData.ua.flatten(),
                            'va':fvData.va.flatten(),
                            'elev':fvData.elev.flatten(),
                            'u':fvData.u,
                            'v':fvData.v}
                else:
                    mod = {'ua':fvData.ua.flatten(),
                            'va':fvData.va.flatten(),
                            'elev':fvData.elev.flatten()}



            obs_loc = {'name':tideFile, 'type':'TideGauge',
                       'mod_time':fvData.time,
                       'obs_time':tideData.time,
                       'lon':tideData.lon,
                       'lat':tideData.lat,
                       'elev_obs_harmonics':tideData.coef,
                       'elev_mod_harmonics': fvElevCoef,
                       'obs_timeseries':tide_obs,
                       'mod_timeseries':mod}


            struct = np.hstack((struct, obs_loc))

    saveName = 'validationStruct.p'
    pickle.dump(struct, open(saveName, "wb"))
    return struct

if __name__ == '__main__':

    #fvFiles = ['/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/']
    fvFiles = ['/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/calibration/bottom_roughness/2D/0.0015/output/dngrid_0001.nc',
     '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/calibration/bottom_roughness/2D/0.0020/output/dngrid_0001.nc',
     '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/calibration/bottom_roughness/2D/0.0025/output/dngrid_0001.nc',
     '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/calibration/bottom_roughness/2D/0.002848/output/dngrid_0001.nc',
     '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/calibration/bottom_roughness/2D/0.0030/output/dngrid_0001.nc']

    adcpFiles = ['/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPa_avg5.mat',
                 '/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPb_avg5.mat']

    fvdebug = '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/dngrid_0001_week2.nc'

    tideFiles = \
    ['/EcoII/EcoEII_server_data_tree/data/observed/GP/TideGauge/Westport_015892_20140325_1212_Z.mat',
     '/EcoII/EcoEII_server_data_tree/data/observed/DG/TideGauge/DigbyWharf_015893_20140115_2221_Z.mat']

    #ind = [-66.3419, -66.3324, 44.2755, 44.2815]
    ind = [-66.3419, -66.3324, 44.2755, 44.2815]
    struct = main(fvFiles, adcpFiles, tideFiles, isStation=False, ax=ind)

