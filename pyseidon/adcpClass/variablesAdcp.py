#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np

class _load_adcp:
    """
'Variables' subset in FVCOM class contains the following numpy arrays:
    """
    def __init__(self,cls, debug=False):
        if debug:
            print 'Loading variables...'

        self.lat = cls.Data['lat']
        self.lon = cls.Data['lon']
        self.north_vel = cls.Data['data']['north_vel']
        self.east_vel = cls.Data['data']['east_vel']
        self.vert_vel = cls.Data['data']['vert_vel']
        self.dir_vel = cls.Data['data']['dir_vel']
        self.mag_signed_vel = cls.Data['data']['mag_signed_vel']
        self.ucross = cls.Data['data']['Ucross']
        self.ualong = cls.Data['data']['Ualong']
        self.pressure = cls.Data['pres']
        self.surf = self.pressure['surf']
        self.time = cls.Data['time']
        self.mtime = self.time['mtime']

        if debug:
            print '...Passed'
