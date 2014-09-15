#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import h5py

class _load_adcp:
    """
'Variables' subset in FVCOM class contains the following numpy arrays:
    """
    def __init__(self,cls, debug=False):
        if debug:
            print 'Loading variables...'
        #Test if load with h5py or scipy
        if not type(cls.Data)==h5py._hl.files.File:
            self.lat = cls.Data['lat']
            self.lon = cls.Data['lon']
            self.bins = cls.Data['data'].bins[:].flatten()
            self.north_vel = cls.Data['data'].north_vel[:]
            self.east_vel = cls.Data['data'].east_vel[:]
            self.vert_vel = cls.Data['data'].vert_vel[:]
            self.dir_vel = cls.Data['data'].dir_vel[:]
            self.mag_signed_vel = cls.Data['data'].mag_signed_vel[:]
            self.pressure = cls.Data['pres']
            self.surf = self.pressure.surf[:].flatten()
            self.time = cls.Data['time']
            self.mtime = self.time.mtime[:].flatten()
            try:
                self.ucross = cls.Data['data'].ucross[:]
                self.ualong = cls.Data['data'].ualong[:]
            except AttributeError:
                pass

        else:
            self.lat = cls.Data['lat'][0][0]
            self.lon = cls.Data['lon'][0][0]
            self.bins = cls.Data['data']['bins'][:].flatten()
            self.north_vel = cls.Data['data']['north_vel'][:].T
            self.east_vel = cls.Data['data']['east_vel'][:].T
            self.vert_vel = cls.Data['data']['vert_vel'][:].T
            self.dir_vel = cls.Data['data']['dir_vel'][:].T
            self.mag_signed_vel = cls.Data['data']['mag_signed_vel'][:].T
            self.pressure = cls.Data['pres']
            self.surf = self.pressure['surf'][:].flatten()
            self.time = cls.Data['time']
            self.mtime = self.time['mtime'][:].flatten()
            try:
                self.ucross = cls.Data['data']['Ucross'][:].T
                self.ualong = cls.Data['data']['Ualong'][:].T
            except KeyError:
                pass

        #Find the depth average of a variable based on percent_of_depth
        #choosen by the user. Currently only working for east_vel (u) and
        #north_vel (v)
        #TR: alaternative with percent of the depth
        #ind = np.argwhere(self.bins < self.percent_of_depth * self.surf[:,None])
        ind = np.argwhere(self.bins < self.surf[:,None])
        index = ind[np.r_[ind[1:,0] != ind[:-1,0], True]]
        data_ma_u = np.ma.array(self.east_vel, mask=np.arange(self.east_vel.shape[1]) > index[:, 1, None])
        data_ma_v = np.ma.array(self.north_vel, mask=np.arange(self.north_vel.shape[1]) > index[:, 1, None])
        self.ua = np.array(data_ma_u.mean(axis=1))
        self.va = np.array(data_ma_v.mean(axis=1))

        if debug:
            print '...Passed'
