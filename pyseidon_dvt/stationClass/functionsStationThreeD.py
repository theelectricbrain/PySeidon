#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

import numpy as np
import numexpr as ne
from pyseidon_dvt.utilities.miscellaneous import *
from pyseidon_dvt.utilities.BP_tools import *
import time

# Custom error
from pyseidon_error import PyseidonError

class FunctionsStationThreeD:
    """
    **'Utils3D' subset of Station class gathers useful functions for 3D runs**
    """
    def __init__(self, variable, grid, plot, History, debug):
        #Inheritance
        self._debug = debug
        self._plot = plot

        #Create pointer to FVCOM class
        setattr(self, '_var', variable)
        setattr(self, '_grid', grid)
        setattr(self, '_History', History)

    def search_index(self, station):
        """Search for the station index"""
        if type(station)==int:
            index = station
        elif type(station).__name__ in ['str', 'ndarray']:
            station = "".join(station).strip().upper()
            for i in range(self._grid.nele):
                if station=="".join(self._grid.name[i]).strip().upper():
                   index=i
        else:
            raise PyseidonError("---Wrong station input---")
        if not 'index' in locals():
            raise PyseidonError("---Wrong station input---")

        return index

    def depth(self, station, debug=False):
        """
        Compute depth at given point

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
          - dep = depth, 2D array (ntime, nlevel)

        Notes:
          - depth convention: 0 = free surface
          - index is used in case one knows already at which
            element depth is requested
        """
        debug = debug or self._debug
        if debug:
            print "Computing depth..."
            start = time.time()

        #Search for the station
        index = self.search_index(station)

        #Compute depth
        h = self._grid.h[index]
        el = self._var.el[:,index]
        zeta = el + h
        siglay = self._grid.siglay[:,index]
        dep = zeta[:, np.newaxis]*siglay[np.newaxis, :]
     
        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        return np.squeeze(dep)

    def verti_shear(self, station, t_start=[], t_end=[],  time_ind=[],
                    bot_lvl=[], top_lvl=[], graph=True, debug=False):
        """
        Compute vertical shear at any given location

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
          - dveldz = vertical shear (1/s), 2D array (time, nlevel - 1)

        Options:
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - time_ind = time indices to work in, list of integers
          - bot_lvl = index of the bottom level to consider, integer
          - top_lvl = index of the top level to consider, integer
          - graph = plots graph if True

        *Notes*
          - use time_ind or t_start and t_end, not both
        """
        debug = debug or self._debug
        if debug:
            print 'Computing vertical shear at point...'

        # Find time interval to work in
        argtime = []
        if not time_ind==[]:
            argtime = time_ind
        elif not t_start==[]:
            if type(t_start)==str:
                start = datetime.datetime.strptime(t_start, '%Y-%m-%d %H:%M:%S')
                end = datetime.datetime.strptime(t_end, '%Y-%m-%d %H:%M:%S')
                argtime = time_to_index(start, end, self._var.julianTime[:], debug=debug)
            else:
                argtime = np.arange(t_start, t_end) 

        #Search for the station
        index = self.search_index(station)

        #Compute depth
        dep = self.depth(station, debug=debug)

        if not argtime==[]:
            depth = dep[argtime,:]
        else:
            depth = dep

        #Sigma levels to consider
        if top_lvl==[]:
            top_lvl = self._grid.nlevel - 1
        if bot_lvl==[]:
            bot_lvl = 0
        sLvl = range(bot_lvl, top_lvl+1)
           
        #Extracting velocity at point
        if not argtime==[]:
            U = self._var.u[argtime,:,index]
            V = self._var.v[argtime,:,index]
        else:
            U = self._var.u[:,:,index]
            V = self._var.v[:,:,index]

        norm = ne.evaluate('sqrt(U**2 + V**2)').squeeze()

        # Compute shear
        dz = depth[:,sLvl[1:]] - depth[:,sLvl[:-1]]
        dvel = norm[:,sLvl[1:]] - norm[:,sLvl[:-1]]           
        dveldz = dvel / dz

        if debug:
            print '...Passed'

        #Plot mean values
        if graph:
            mean_depth = np.mean((depth[:,sLvl[1:]]
                       + depth[:,sLvl[:-1]]) / 2.0, 0)
            mean_dveldz = np.mean(dveldz,0)
            error = np.std(dveldz,axis=0)
            self._plot.plot_xy(mean_dveldz, mean_depth, xerror=error[:],
                               title='Shear profile ',
                               xLabel='Shear (1/s) ', yLabel='Depth (m) ')

        return np.squeeze(dveldz)

    def velo_norm(self, station, t_start=[], t_end=[], time_ind=[],
                  graph=True, debug=False):
        """
        Compute the velocity norm at any given location

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
          - velo_norm = velocity norm, 2D array (time, level)

        Options:
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - time_ind = time indices to work in, list of integers
          - graph = plots vertical profile averaged over time if True

        *Notes*
          - use time_ind or t_start and t_end, not both
        """
        debug = debug or self._debug
        if debug:
            print 'Computing velocity norm at point...'
       
        # Find time interval to work in
        argtime = []
        if not time_ind==[]:
            argtime = time_ind
        elif not t_start==[]:
            if type(t_start)==str:
                start = datetime.datetime.strptime(t_start, '%Y-%m-%d %H:%M:%S')
                end = datetime.datetime.strptime(t_end, '%Y-%m-%d %H:%M:%S')
                argtime = time_to_index(start, end, self._var.julianTime[:], debug=debug)
            else:
                argtime = np.arange(t_start, t_end)

        #Search for the station
        index = self.search_index(station)

        #Computing velocity norm
        try:
            if not argtime==[]:          
                U = self._var.u[argtime, :, index]
                V = self._var.v[argtime, :, index]
                W = self._var.w[argtime, :, index]
                velo_norm = ne.evaluate('sqrt(U**2 + V**2 + W**2)').squeeze()
            else:            
                U = self._var.u[:, :, index]
                V = self._var.v[:, :, index]
                W = self._var.w[:, :, index]
                velo_norm = ne.evaluate('sqrt(U**2 + V**2 + W**2)').squeeze()
        except AttributeError:
            if not argtime==[]:          
                U = self._var.u[argtime, :, index]
                V = self._var.v[argtime, :, index]
                velo_norm = ne.evaluate('sqrt(U**2 + V**2)').squeeze()
            else:            
                U = self._var.u[:, :, index]
                V = self._var.v[:, :, index]
                velo_norm = ne.evaluate('sqrt(U**2 + V**2)').squeeze()

        if debug:
            print '...passed'

        #Plot mean values
        if graph:
            depth = np.mean(self.depth(station),axis=0)
            vel = np.mean(velo_norm,axis=0)
            error = np.std(velo_norm,axis=0)
            self._plot.plot_xy(vel, depth, xerror=error[:],
                               title='Velocity norm profile ',
                               xLabel='Velocity (m/s) ', yLabel='Depth (m) ')
      

        return velo_norm

    def flow_dir(self, station, t_start=[], t_end=[], time_ind=[], 
                          vertical=True, debug=False):
        """
        Compute flow directions and associated norm at any given location.

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
          - flowDir = flowDir at (pt_lon, pt_lat), 2D array (ntime, nlevel)

        Options:
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - time_ind = time indices to work in, list of integers
          - vertical = True, compute flowDir for each vertical level

        *Notes*
          - directions between -180 and 180 deg., i.e. 0=East, 90=North,
            +/-180=West, -90=South
          - use time_ind or t_start and t_end, not both
        """
        debug = debug or self._debug
        if debug:
            print 'Computing flow directions at point...'

        #Search for the station
        index = self.search_index(station)

        # Find time interval to work in
        argtime = []
        if not time_ind==[]:
            argtime = time_ind
        elif not t_start==[]:
            if type(t_start)==str:
                start = datetime.datetime.strptime(t_start, '%Y-%m-%d %H:%M:%S')
                end = datetime.datetime.strptime(t_end, '%Y-%m-%d %H:%M:%S')
                argtime = time_to_index(start, end, self._var.julianTime[:], debug=debug)
            else:
                argtime = np.arange(t_start, t_end)
        
        #Choose the right pair of velocity components
        if not argtime==[]:
            if self._var._3D and vertical:
                u = self._var.u[argtime,:,index]
                v = self._var.v[argtime,:,index]
            else:
                u = self._var.ua[argtime,index]
                v = self._var.va[argtime,index]    

        #Compute directions
        if debug: print 'Computing arctan2 and norm...'
        dirFlow = np.rad2deg(np.arctan2(v,u))
       
        if debug:
                print '...Passed'

        return np.squeeze(dirFlow)

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
