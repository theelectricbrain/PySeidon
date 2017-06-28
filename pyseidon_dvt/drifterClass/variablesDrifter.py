#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
# from numpy.ma import MaskError
# import h5py
from pyseidon_dvt.utilities.miscellaneous import mattime_to_datetime
import sys

class _load_drifter:
    """
    **'Variables' subset in Tidegauge class**
    It contains the following numpy arrays: ::

                         _u = u velocity component (m/s), 1D array (ntime)
                        |_v = v velocity component (m/s), 1D array (ntime)
     Drifter.Variables._|_matlabTime = matlab time, 1D array (ntime)
                        |_lon = longitudes (deg.), 1D array (ntime)
                        |_lat = latitudes (deg.), 1D array (ntime)

    """
    def __init__(self,cls, History, debug=False):
        if debug:
            print 'Loading variables...'
        # Pointer to History
        setattr(self, '_History', History)
        try:
            self.matlabTime = cls.Data['velocity'].vel_time[:]
            #Sorting values with increasing time step
            sortedInd = self.matlabTime.argsort()
            self.matlabTime.sort()
            self.lat = cls.Data['velocity'].vel_lat[sortedInd]
            self.lon = cls.Data['velocity'].vel_lon[sortedInd]
            self.u = cls.Data['velocity'].u[sortedInd]
            self.v = cls.Data['velocity'].v[sortedInd]
        # Luna Ocean Consulting Ltd. new drifter format
        except KeyError:
            self.matlabTime = cls.Data['time']
            self.original = {}
            self.original['lat'] = cls.Data['lat']
            self.original['lon'] = cls.Data['lon']
            self.original['u'] = cls.Data['u']
            self.original['v'] = cls.Data['v']
            self.original['drift_start'] = cls.Data['drift_start']
            self.original['drift_stop'] = cls.Data['drift_stop']
            self.smooth = {}
            self.smooth['lat'] = cls.Data['lat_smooth']
            self.smooth['lon'] = cls.Data['lon_smooth']
            self.smooth['u'] = cls.Data['u_smooth']
            self.smooth['v'] = cls.Data['v_smooth']
            self.smooth['drift_start'] = cls.Data['drift_start_smooth']
            self.smooth['drift_stop'] = cls.Data['drift_stop_smooth']
        except:
            sys.exit('Drifter file format incompatible')

        #-Append message to History field
        try:
            start = mattime_to_datetime(self.matlabTime[0])
            end = mattime_to_datetime(self.matlabTime[-1])
            text = 'Temporal domain from ' + str(start) +\
                    ' to ' + str(end)
            self._History.append(text)
        except ValueError:
            start = mattime_to_datetime(np.nanmin(self.matlabTime))
            end = mattime_to_datetime(np.nanmax(self.matlabTime))
            text = 'Temporal domain from ' + str(start) +\
                    ' to ' + str(end)
            self._History.append(text)

        if debug: print '...Passed'

        return
