#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import sys
import numexpr as ne
from datetime import datetime
from datetime import timedelta
from interpolation_utils import *
from miscellaneous import *
import time

#TR comment: This all routine needs to be tested and debugged
class FunctionsFvcomThreeD:
    """'UtilsThreeD' subset of FVCOM class gathers useful functions for 3D runs"""
    def __init__(self, variable, grid, plot, util, QC, debug):
        #Inheritance
        self._debug = debug
        self._var = variable
        self._grid = grid
        self._plot = plot
        self._QC = QC
        self._util = util
        self.interpolation_at_point = self._util.interpolation_at_point
        self.hori_velo_norm = self._util.hori_velo_norm

        #Create pointer to FVCOM class
        variable = self._var
        grid = self._grid
        QC = self._QC

    def depth(self, debug=False):
        """
        Compute new grid variable 'depth' -> FVCOM.Grid.depth

        Notes:
        -----
          - Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            print "Computing depth..."
            start = time.time()

        #Compute depth
        #Different interpolation method if partial dat
        if hasattr(self._grid, '_region_e'):
            h = self._grid.h[self._grid_region_e]
        else:
            h = self._grid.h
        zeta = self._var.el[:,:] + h[None,:]
        siglay = self._grid.siglay[:]
        z = zeta[:,None,:]*siglay[None,:,:]
        dep = np.zeros([self._grid.ntime,self._grid.nlevel,self._grid.node])
        for i in range(self._grid.node):
            dep[:,:,i] = interpolation_at_point(z[:,:,i],
                                                self._grid.lonc[i],
                                                self._grid.latc[i], debug=debug)

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        # Add metadata entry
        self._grid.depth = dep
        self._QC.append('depth computed')
        print '-Flow directions added to FVCOM.Variables.-'

    def depth_at_point(self, pt_lon, pt_lat, debug=False):
        """
        Compute depth at given point

        Inputs:
        ------
          - pt_lon = longitude in degrees to find
          - pt_lat = latitude in degrees to find
        """
        debug = debug or self._debug
        if debug:
            print "Computing depth"
            start = time.time()

        #Finding index        
        index = closest_point([pt_lon], [pt_lat],
                              self._grid.lonc,
                              self._grid.latc, debug=debug)[0]

        if not hasattr(self._grid, 'depth'):
            #Compute depth
            h = self._grid.h[index]
            zeta = self._var.el[:,index] + h
            siglay = self._grid.siglay[:]
            z = zeta[:,None]*siglay[None,:,index]
            dep = np.zeros([self._grid.ntime,self._grid.nlevel])
            dep = self.interpolation_at_point(z, pt_lon,
                                                 pt_lat,
                                                 debug=debug)
        else:
            dep = self._grid.depth[:,:,index]         
        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        return dep

    def verti_shear(self, debug=False):
        """
        Compute new variable 'vertical shear' -> FVCOM.Variables.verti_shear

        Inputs:
        ------
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'),
            or time index (integer)
          - t_end = end time, as string ('yyyy-mm-ddThh:mm:ss'),
            or time index (integer)
        Keywords:
        --------
          - bot_lvl = index of the bottom level to consider, integer
          - top_lvl = index of the top level to consider, integer
        Notes:
        -----
          - Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            print 'Computing vertical shear...'
              
        #Compute depth if necessary
        if not hasattr(self._grid, 'depth'):        
           depth = self.depth(debug=debug)
        depth = self._grid.depth

        # Checking if horizontal velocity norm already exists
        if not hasattr(self._var, 'hori_velo_norm'):
            self.hori_velo_norm()
        vel = self._var.hori_velo_norm

        #Sigma levels to consider
        top_lvl = (self._grid.nlevel) - 1
        bot_lvl = 0
        sLvl = range(bot_lvl, top_lvl+1)

        # Compute shear
        dz = depth[:,sLvl[1:],:] - depth[:,sLvl[:-1],:]
        dvel = vel[:,sLvl[1:],:] - vel[:,sLvl[:-1],:]           
        dveldz = dvel / dz

        #Custom return
        self._var.verti_shear = dveldz 
            
        # Add metadata entry
        self._QC.append('vertical shear computed')
        print '-Vertical shear added to FVCOM.Variables.-'

        if debug:
            print '...Passed'

    def verti_shear_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[],  time_ind=[],
                             bot_lvl=[], top_level= [], graph=True, debug=False):
        """
        Compute vertical shear -> FVCOM.Variables.verti_shear
        Inputs:
        ------
          - pt_lon = longitude in degrees to find
          - pt_lat = latitude in degrees to find
        Keywords:
        --------
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - t_end = end time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - time_ind = time indexes to work in, list of integers
          - bot_lvl = index of the bottom level to consider, integer
          - top_lvl = index of the top level to consider, integer
          - graph = plot graph if True
        """
        debug = debug or self._debug
        if debug:
            print 'Computing vertical shear at point...'

        # Find time interval to work in
        argtime = []
        if time_ind:
            argtime = time_ind
        elif t_start:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                argtime = arange(t_start, t_end) 

        #Compute depth
        dep = self.depth_at_point(pt_lon, pt_lat, debug=debug)
        if argtime:
            depth = dep[argtime,:]
        else:
            depth = dep

        #Sigma levels to consider
        if not top_lvl:
            top_lvl = (self._grid.nlevel) - 1
        if not bot_lvl:
            bot_lvl = 0
        sLvl = range(bot_lvl, top_lvl+1)


        # Checking if vertical shear already exists
        if not hasattr(self._var, 'verti_shear'): 
             
            #Extracting velocity at point
            if argtime:
                u = self._var.u[argtime,:,:]
                v = self._var.v[argtime,:,:]
            else:
                u = self._var.u
                v = self._var.v

            #Extraction at point
            if debug:
                print 'Extraction of u and v at point...'
            U = self.interpolation_at_point(u, pt_lon, pt_lat,
                                            debug=debug)  
            V = self.interpolation_at_point(v, pt_lon, pt_lat,
                                            debug=debug)
            norm = ne.evaluate('sqrt(U**2 + V**2)')     

            # Compute shear
            dz = depth[:,sLvl[1:]] - depth[:,sLvl[:-1]]
            dvel = norm[:,sLvl[1:]] - norm[:,sLvl[:-1]]           
            dveldz = dvel / dz
        else:
            dveldz = interpolation_at_point(self._var.verti_shear,
                                            pt_lon, pt_lat, debug=debug)

        if debug:
            print '...Passed'

        #Plot mean values
        if graph:
            mean_depth = np.mean((depth[:,sLvl[1:]]
                       + depth[:,sLvl[:-1]]) / 2.0, 0)
            mean_dveldz = np.mean(dveldz,0)
            self._plot.plot_xy(mean_dveldz, mean_depth, title='Shear profile ',
                               xLabel='Shear (1/s) ', yLabel='Depth (m) ')

        return dveldz             

    def velo_norm(self, debug=False):
        """
        Compute new variable 'velocity norm' -> FVCOM.Variables.velo_norm

        Notes:
        -----
          -Can take time over the full domain
        """
        if debug or self._debug:
            print 'Computing velocity norm...'

        #Computing velocity norm
        u = self._var.u[:, :, :]
        v = self._var.v[:, :, :]
        ww = self._var.ww[:, :, :]
        vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')

        #Custom return    
        self._var.velo_norm = vel 
       
        # Add metadata entry
        self._QC.append('Velocity norm computed')
        print '-Velocity norm added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def velo_norm_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[], time_ind=[],
                           debug=False):
        """
        Compute vertical shear -> FVCOM.Variables.verti_shear
        Inputs:
        ------
          - pt_lon = longitude in degrees to find
          - pt_lat = latitude in degrees to find
        Keywords:
        --------
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - t_end = end time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - time_ind = time indexes to work in, list of integers
        """
         debug = debug or self._debug
        if debug:
            print 'Computing velocity norm at point...'
       
        # Find time interval to work in
        argtime = []
        if time_ind:
            argtime = time_ind
        elif t_start:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                argtime = arange(t_start, t_end)

        #computing velocity norm
        if argtime:
            if not hasattr(self._var, 'velo_norm'):             
                u = self._var.u[argtime, :, :]
                v = self._var.v[argtime, :, :]
                ww = self._var.ww[argtime, :, :]
                vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')
            else:
                vel = self._var.velo_norm[argtime, :, :]
        else:
            if not hasattr(self._var, 'velo_norm'):             
                u = self._var.u
                v = self._var.v
                ww = self._var.ww
                vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')
            else:
                vel = self._var.velo_norm

       #Interpolation
       velo_norm = self.interpolation_at_point(vel, pt_lon, pt_lat,
                                               debug=debug)

        if debug:
            print '...passed'

       return velo_norm 


    def flow_dir_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[], time_ind=[],
                          vertical=True, debug=False):
        """
        Flow directions and associated norm at any give location.

        Inputs:
        ------
          - pt_lon = longitude in degrees to find
          - pt_lat = latitude in degrees to find
        Outputs:
        -------
          - flowDir = flowDir at (pt_lon, pt_lat)
          - norm = velocity norm at (pt_lon, pt_lat)
        Keywords:
        -------
          - time_ind = time indexes to work in, list of integers
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - t_end = end time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - vertical = True, compute flowDir for each vertical level
        Notes:
        -----
          - directions between 0 and 360 deg.
        """
        debug = debug or self._debug
        if debug:
            print 'Computing flow directions at point...'

        # Find time interval to work in
        argtime = []
        if time_ind:
            argtime = time_ind
        elif t_start:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                argtime = arange(t_start, t_end)
        
        #Checking
        if not hasattr(self._var, 'dir_flow'):
            #Choose the right pair of velocity components
            if argtime:
                if self._var._3D and vertical:
                    u = self._var.u[argtime,:,:]
                    v = self._var.v[argtime,:,:]
                else:
                    u = self._var.ua[argtime,:]
                    v = self._var.va[argtime,:]
            else:
                if self._var._3D and vertical:
                   u = self._var.u
                   v = self._var.v
                else:
                   u = self._var.ua
                   v = self._var.va

            #Extraction at point
            if debug:
                print 'Extraction of u and v at point...'
            U = self._util.interpolation_at_point(u, pt_lon, pt_lat,
                                                  debug=debug)  
            V = self._util.interpolation_at_point(v, pt_lon, pt_lat,
                                                  debug=debug)       

            #Compute directions
            if debug:
                print 'Computing arctan2 and norm...'
            dirFlow = np.arctan2(U,V)
            dirFlow = np.mod(((np.pi/2.0) - dirFlow) * (180.0 / np.pi), 360.0)

        else:
            if argtime:
                dir_flow = self._var.dir_flow[argtime,:,:]
                dirFlow = self._util.interpolation_at_point(dir_flow, pt_lon, pt_lat,
                                                            debug=debug)   
            else:
                dirFlow = self._util.interpolation_at_point(self._var.dir_flow,
                                                            pt_lon, pt_lat, debug=debug) 
         
        if debug:
                print '...Passed'

        return dirFlow, norm

    def flow_dir(self, debug=False):
        """"
        Compute new variable 'flow directions' -> FVCOM.Variables.flow_dir
        Notes:
        -----
          - directions between 0 and 360 deg.
          - Can take time over the full domain
        """
        if debug or self._debug:
            print 'Computing flow directions...'

        u = self._var.u
        v = self._var.v

        dirFlow = (np.pi/2.0) - np.arctan2(u,v)
        dirFlow = dirFlow * (180.0 / np.pi)

        #Custom return    
        self._var.dir_flow = dirFlow 

        # Add metadata entry
        self._QC.append('flow directions computed')
        print '-Flow directions added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

