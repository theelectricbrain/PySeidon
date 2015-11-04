#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from datetime import datetime
from datetime import timedelta
import fnmatch
import os
import time
from scipy.io import netcdf
from pydap.client import open_url

# Custom error
from pyseidon_error import PyseidonError

def date2py(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + \
        timedelta(days=matlab_datenum%1) - timedelta(days = 366)

    return python_datetime

def op_angles_from_vectors(u, v, debug=False):
    """
    This function takes in vectors in the form (u,v) and compares them in
    order to find the angles of the vectors without any wrap-around issues.
    This is accomplished by finding the smallest difference between angles
    compared at different wrap-around values.
    This appears to work correctly.

    Inputs:
      -u = velocity component along x (West-East) direction, 1D array
      -v = velocity component along y (South-North) direction, 1D array

    Outputs:
      -angle = corresponidng angle in degrees, 1D array

    Notes:
      -Angles are reported in compass coordinates, i.e. 0 and 360 deg.,
       0/360=East, 90=North, 180=West, 270=South
    """
    if debug:
        print 'Computing angles from velocity component...'
        start = time.time()

    phi = np.mod((-1.0*np.arctan2(v,u)) * (180.0/np.pi) + 90.0, 360.0)
    if len(phi.shape)==1:#Assuming the only dimension is time
        #Compute difference between angles
        diff1 = np.abs(phi[:-1]-phi[1:]) #initial difference between angles
        diff2 = np.abs(phi[:-1]-phi[1:]-360.0) #diff when moved down a ring
        diff3 = np.abs(phi[:-1]-phi[1:]+360.0) #diff when moved up a ring

        index1 = np.where((diff2 < diff1) & (diff2 < diff3))[0]
        index2 = np.where((diff3 < diff1) & (diff3 < diff2))[0]

        phi[index1] = np.mod(phi[index1] - 360.0, 360.0)
        phi[index2] = np.mod(phi[index2] + 360.0, 360.0)
    elif len(phi.shape)==2:#Assuming the only dimension is time and sigma level
        #Compute difference between angles
        diff1 = np.abs(phi[:-1,:]-phi[1:,:]) #initial difference between angles
        diff2 = np.abs(phi[:-1,:]-phi[1:,:]-360.0) #diff when moved down a ring
        diff3 = np.abs(phi[:-1,:]-phi[1:,:]+360.0) #diff when moved up a ring

        index1 = np.where((diff2 < diff1) & (diff2 < diff3))[0]
        index2 = np.where((diff3 < diff1) & (diff3 < diff2))[0]

        phi[index1] = phi[index1] - 360.0
        phi[index2] = phi[index2] + 360.0
    else: #Assuming the only dimension is time ,sigma level and element
        #Compute difference between angles
        diff1 = np.abs(phi[:-1,:,:]-phi[1:,:,:]) #initial difference between angles
        diff2 = np.abs(phi[:-1,:,:]-phi[1:,:,:]-360.0) #diff when moved down a ring
        diff3 = np.abs(phi[:-1,:,:]-phi[1:,:,:]+360.0) #diff when moved up a ring

        index1 = np.where((diff2 < diff1) & (diff2 < diff3))[0]
        index2 = np.where((diff3 < diff1) & (diff3 < diff2))[0]

        phi[index1] = phi[index1] - 360.0
        phi[index2] = phi[index2] + 360.0     

    if debug:
        end = time.time()
        print "...processing time: ", (end - start)
    
    return phi

def time_to_index(t_start, t_end, time, debug=False):
    """
    Convert datetime64[us] string in FVCOM index

    Inputs:
      - t_start = start time in datetime
      - t_end = end time in datetime
      - time = array of julian days

    Outputs:
      - argtime = arry of indices
    """
    start = date_to_julian_day(t_start)
    end = date_to_julian_day(t_end)

    t_slice = [start, end]

    argtime = np.argwhere((time>=t_slice[0])&(time<=t_slice[-1])).ravel()
    if debug:
        print 'Argtime: ', argtime
    if argtime == []:
        raise PyseidonError("Wrong time input")
    return argtime

def mattime_to_datetime(mattime, debug=False):
    """Convert matlab time to datetime64[us] """
    date = datetime.fromordinal(int(mattime)) + \
               timedelta(days=mattime%1)-timedelta(days=366)
    time = np.array(date,dtype='datetime64[us]')

    return time

def datetime_to_mattime(dt, debug=False):
    """Convert datetime64[us] to matlab time"""
    mdn = dt + timedelta(days = 366)
    s = (dt.hour * (60.0*60.0)) + (dt.minute * 60.0) + dt.second
    day = 24.0*60.0*60.0
    frac = s/day

    return mdn.toordinal() + frac

def findFiles(filename, name):
    """
    Wesley comment[elements] the name needs to be a linux expression to find files
    you want. For multiple station files, this would work
    name = '*station*.nc'

    For just dngrid_0001 and no restart files:
    name = 'dngrid_0*.nc'
    will work
    """

    name = '*' + name + '*.nc'
    matches = []
    for root, dirnames, filenames in os.walk(filename):
        for filename in fnmatch.filter(filenames, name):
            matches.append(os.path.join(root, filename))
            filenames.remove(filename)
        for filename in fnmatch.filter(filenames, name.lower()):
            matches.append(os.path.join(root, filename))

    return sorted(matches)

def date_to_julian_day(my_date):
    """Returns the Julian day number of a date."""
    # a = (14 - my_date.month)//12
    # y = my_date.year + 4800 - a
    # m = my_date.month + 12*a - 3
    # s = (my_date.hour * (60.0*60.0)) + (my_date.minute * 60.0) + my_date.second
    # day = 24.0*60.0*60.0
    # jtime = my_date.day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045 + s/day
    jtime = datetime_to_mattime(my_date) - 678942.0
    return jtime    
    
def distance(locs,loce):
    """Returns the distance in meters between two locations in long/lat."""
    TPI=111194.92664455874
    y0c = TPI * (loce[1] - locs[1])
    dx_sph = loce[0] - locs[0]
    if (dx_sph > 180.0):
        dx_sph=dx_sph-360.0
    elif (dx_sph < -180.0):
        dx_sph =dx_sph+360.0
    x0c = TPI * np.cos(np.deg2rad(loce[1] + locs[1])*0.5) * dx_sph

    dist=np.linalg.norm([x0c, y0c])

    return dist
