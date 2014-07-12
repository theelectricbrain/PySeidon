#!/usr/bin/python2.7
# encoding: utf-8

class _load_adcp:
    """
'Variables' subset in FVCOM class contains the following numpy arrays:
    """
    def _load(self, filename, debug=False):
        if debug:
            print 'Loading variables...'

        self.lat = self.mat['lat']
        self.lon = self.mat['lon']

        self.data = self.mat['data']
        self.north_vel = self.data.north_vel
        self.east_vel = self.data.east_vel
        self.vert_vel = self.data.vert_vel
        self.dir_vel = self.data.dir_vel
        self.mag_signed_vel = self.data.mag_signed_vel
        self.ucross = self.data.ucross
        self.ualong = self.data.ualong

        self.pressure = self.mat['pres']
        self.surf = self.pressure.surf

        self.time = self.mat['time']
        self.mtime = self.time.mtime

        if debug:
            print '...Passed'
