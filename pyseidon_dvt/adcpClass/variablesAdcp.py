#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from numpy.ma import MaskError
import h5py
from pyseidon_dvt.utilities.miscellaneous import mattime_to_datetime

# Custom error
from pyseidon_error import PyseidonError

class _load_adcp:
    """
    **'Variables' subset in ADCP class**

    It contains the following numpy arrays: ::

                      _bins = depth of measurement bins, 1D array, shape=(bins)
                     |_depth = depth, negative from surface down, time serie, 2D array, shape=(time,bins)
                     |_dir_vel = velocity direction time serie, 2D array, shape=(time,bins)
                     |_east_vel = East velocity time serie, 2D array, shape=(time,bins)
                     |_lat = latitude, float, decimal degrees
                     |_lon = lontitude, float, decimal degrees
                     |_mag_signed_vel = signed velocity time serie, 2D array, shape=(time,bins)
                     |_matlabTime = matlab time, 1D array, shape=(time)
     ADCP.Variables._|_north_vel = East velocity time serie, 2D array, shape=(time,bins)
                     |_percent_of_depth = percent of the water column measured by ADCP, float
                     |_surf = pressure at surface timeserie, 1D array, shape=(time)
                     |_el = elevation (m) at surface timeserie, 1D array, shape=(time)
                     |_ua = depth averaged velocity component timeserie, 1D array, shape=(time)
                     |_va = depth averaged velocity component timeserie, 1D array, shape=(time)
                     |_ucross = ???, 1D array, shape=(time)
                     |_ualong = ???, 1D array, shape=(time)

    """
    def __init__(self,cls, History, debug=False):
        if debug:
            print 'Loading variables...'
        # Pointer to History
        setattr(self, '_History', History)


        # TR: fudge factor, squeeze out the 5 top % of the water column
        self.percent_of_depth=0.95
       
        # Convert some mat_struct objects into a dictionaries.
        # This will enable compatible syntax to h5py files.
        if not type(cls.Data) == h5py._hl.files.File:
            try:
                cls.Data['data'] = cls.Data['data'].__dict__
            except KeyError:
                raise PyseidonError("Missing 'data' field from ADCP file.")
            for key in cls.Data['data']:
                if key is not '_fieldnames':
                    # in ther is a mat_struct in the mat_struct, like 'surf'
                    if hasattr(cls.Data['data'][key], '__module__') and \
                       cls.Data['data'][key].__module__ == 'scipy.io.matlab.mio5_params':
                        cls.Data['data'][key] = cls.Data['data'][key].__dict__
                        for kk in cls.Data['data'][key]:
                            if kk is not '_fieldnames':
                                cls.Data['data'][key][kk] = cls.Data['data'][key][kk].T
                    else:
                        cls.Data['data'][key] = cls.Data['data'][key].T
            
            # Convert other fields to dictionaries if they exist.
            if 'pres' in cls.Data:
                cls.Data['pres']=cls.Data['pres'].__dict__
            if 'time' in cls.Data:
                cls.Data['time']=cls.Data['time'].__dict__
            
            
        try:
            self.lat = np.ravel(cls.Data['lat'])[0]
            self.lon = np.ravel(cls.Data['lon'])[0]
        except KeyError:
            if debug:
                print 'Missing lon and/or lat data'
                
        try:
            self.bins = cls.Data['data']['bins'][:].ravel()
        except KeyError:
            if debug:
                print 'Missing bins'   
                             
        try:
            self.north_vel = cls.Data['data']['north_vel'][:].T
            self.east_vel = cls.Data['data']['east_vel'][:].T
            self.v = self.north_vel
            self.u = self.east_vel
            try:
                self.north_vel=np.ma.masked_array(self.north_vel,np.isnan(self.north_vel))
                self.east_vel=np.ma.masked_array(self.east_vel,np.isnan(self.east_vel))
                self.v=self.north_vel
                self.u=self.east_vel
            except MaskError:
                print 'Failed to mask horizontal velocities (north_vel, east_vel)'
        except KeyError:
            if debug:
                print 'Missing horizontal velocities (north_vel, east_vel)'
                                
        try:
            self.vert_vel = cls.Data['data']['vert_vel'][:].T
            try:
                self.vert_vel=np.ma.masked_array(self.vert_vel,np.isnan(self.vert_vel))
            except MaskError:
                print 'Failed to mask vertical velocity (vert_vel)'
        except KeyError:
            if debug:
                print 'Missing vertical velocity (vert_vel)'  
                              
        try:
            self.dir_vel = cls.Data['data']['dir_vel'][:].T
            try:
                self.dir_vel=np.ma.masked_array(self.dir_vel,np.isnan(self.dir_vel))
            except MaskError:
                print 'Failed to mask dir_vel'
        except KeyError:
            if debug:
                print 'Missing dir_vel'  
                                            
        try:
            self.mag_signed_vel = cls.Data['data']['mag_signed_vel'][:].T
            try:
                self.mag_signed_vel=np.ma.masked_array(self.mag_signed_vel,np.isnan(self.mag_signed_vel))
            except MaskError:
                print 'Failed to mask mag_signed_vel'
        except KeyError:
            if debug:
                print 'Missing mag_signed_vel'   
                                           
        try:
            self.pressure = cls.Data['pres']
            self.surf = self.pressure['surf'][:].ravel()
            self.el = self.surf
        except KeyError:
            if debug:
                print 'Missing elevation data (pres.surf)' 
                                             
        try:
            self.time = cls.Data['time']
            self.matlabTime = self.time['mtime'][:].ravel()
        except KeyError:
            if debug:
                print 'Missing time data (time.mtime)' 
                                             
        try:
            self.ucross = cls.Data['data']['Ucross'][:].T
            self.ualong = cls.Data['data']['Ualong'][:].T
        except KeyError:
            if debug:
                print 'Missing along/cross velocities (Ucross, Ualong)'
                
        try:
            self.ua = cls.Data['data']['ua'][:].T
            self.va = cls.Data['data']['va'][:].T
            try:
                self.ua=np.ma.masked_array(self.ua,np.isnan(self.ua))
                self.va=np.ma.masked_array(self.va,np.isnan(self.va))
            except MaskError:
                print 'Failed to mask depth averaged velocities (ua, va)'
        except KeyError:
            if debug:
                print 'Missing depth averaged velocities (ua, va)'


        #-Append message to History field
        try:
            start = mattime_to_datetime(self.matlabTime[0])
            end = mattime_to_datetime(self.matlabTime[-1])
            text = 'Temporal domain from ' + str(start) +\
                    ' to ' + str(end)
            self._History.append(text)
        except AttributeError:
            if debug:
                print 'Missing time variable failed to add history note'

        #Find the depth average of a variable based on percent_of_depth
        #choosen by the user but only if not loaded from file directly.
        # Currently only working for east_vel (u) and north_vel (v)
        if (('ua' not in dir(self)) and ('va' not in dir(self))):
            try:
                #TR: alternative with percent of the depth
                ind = np.argwhere(self.bins < self.percent_of_depth * self.surf[:,np.newaxis])
                #ind = np.argwhere(self.bins < self.surf[:,np.newaxis])
                index = ind[np.r_[ind[1:,0] != ind[:-1,0], True]]
                try:
                    data_ma_u = np.ma.array(self.east_vel,
                                mask=np.arange(self.east_vel.shape[1]) > index[:, 1, np.newaxis])
                    data_ma_u=np.ma.masked_array(data_ma_u,np.isnan(data_ma_u))
                except MaskError:
                    data_ma_u=np.ma.masked_array(self.east_vel,np.isnan(self.east_vel))

                try:
                    data_ma_v = np.ma.array(self.north_vel,
                                mask=np.arange(self.north_vel.shape[1]) > index[:, 1, np.newaxis])
                    data_ma_v=np.ma.masked_array(data_ma_v,np.isnan(data_ma_v))
                except MaskError:
                    data_ma_v=np.ma.masked_array(self.north_vel,np.isnan(self.north_vel))
                
                self.ua = np.array(data_ma_u.mean(axis=1))
                self.va = np.array(data_ma_v.mean(axis=1))
            except AttributeError:
                if debug:
                    print 'Missing atleast one variable required ' + \
                          'to compute depth averaged velocities'
        

        # Compute depth with fvcom convention, negative from surface down
        try:
            self.depth = np.ones(self.north_vel.shape) * np.nan
            for t in range(self.matlabTime.shape[0]):
                #i = np.where(np.isnan(self.north_vel[t,:]))[0][0]
                #z = self.bins[i]
                self.depth[t, :] = self.bins[:] - self.surf[t]
            self.depth[np.where(self.depth>0.0)] = np.nan
        except AttributeError:
            if debug:
                print 'Missing atleast one variable required ' + \
                  'to compute depth with fvcom convention'

        if debug:
            print '...Passed'

        return
