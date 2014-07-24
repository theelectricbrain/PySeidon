#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import sys
import numexpr as ne
from datetime import datetime
from datetime import timedelta
from interpolation_utils import *
import time

#TR comment: This all routine needs to be tested and debugged
class FunctionsFvcomThreeD:
    """'UtilsThreeD' subset of FVCOM class gathers useful functions for 3D runs"""
    def __init__(self, cls):
        #Inheritance
        self._debug = cls._debug
        self._var = cls.Variables
        self._grid = cls.Grid
        self._plot = cls.Plots
        self._QC = cls.QC
        self._util = cls.Utils
        self.interpolation_at_point = self._util.interpolation_at_point
        #Create pointer to FVCOM class
        cls.Variables = self._var
        cls.Grid = self._grid
        cls.QC = self._QC

    def depth(self, debug=False):
        """
        Compute depth -> FVCOM.Grid.depth
        Notes:
        -----
           Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            print "Computing depth"
            start = time.time()

        #Compute depth
        h = self._grid.h
        zeta = self._var.el[:,:] + h[None,:]
        nv = self._grid.trinodes
        siglay = self._grid.siglay[:]
        z = zeta[:,None,:]*siglay[None,:,:]
        dep = np.zeros([zeta.shape[0],siglay.shape[0],nv.shape[0]])
        #for i in range(z.shape[0]):
        #    for j in range(nv.shape[0]):
        #        dep[i,:,j] = (z[i,:,nv[j,0]]+z[i,:,nv[j,1]]+
        #                      z[i,:,nv[j,2]])/3
        #TR alternative: over the whole thing
        #dep = (z[:,:,nv[:,0]] + z[:,:,nv[:,1]] + z[:,:,nv[:,2]]) / 3
        #end = time.time()
        #TR comment: I have doubt on this interp approach
        #print "Computation time method1: ", (end - start)            
        #TR alternative2: using the interp function
        for i in range(nv.shape[0]):
            dep[:,:,i] = interpN_at_pt(z, self._grid.xc[i],
                                       self._grid.yc[i],
                                       self._grid.xc,
                                       self._grid.yc, i,
                                       self._grid.trinodes,
                                       self._grid.aw0,
                                       self._grid.awx,
                                       self._grid.awy)
        #end = time.time()
        #print "Computation time method1: ", (end - start)  
        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        self._grid.depth = dep

    def verti_shear(self, t_start, t_end,
                    bot_lvl=[], top_level= [], debug=False):
        """
        Compute vertical shear -> FVCOM.Variables.verti_shear
        Inputs:
        ------
          t_start = start time, datetime64[us] type
          t_end = end time, datetime64[us] type
        Keywords:
        --------
          bot_lvl = index of the bottom level to consider, integer
          top_lvl = index of the top level to consider, integer
        Notes:
        -----
           Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            print 'Computing vertical shear...'
        if self._var._3D:
            # Find simulation time contains in [t_start, t_end]
            time = self._var.matlabTime
            t = time.shape[0]
            l = []
            for i in range(t):
                date = datetime.fromordinal(int(time[i])) + \
                       timedelta(days=time[i]%1)-timedelta(days=366)
                l.append(date)
            time = np.array(l,dtype='datetime64[us]')
            t_slice = [t_start, t_end]
            t_slice = np.array(t_slice,dtype='datetime64[us]')

            if t_slice.shape[0] != 1:
                argtime = np.argwhere((time>=t_slice[0])&
                                      (time<=t_slice[-1])).flatten()
            if debug or self._debug:
                print argtime
            
            #Compute depth
            dep = self.depth(debug=debug)
            # Checking if horizontal velocity norm already exists
            if not hasattr(self._var, 'hori_velo_norm'):
                self.hori_velo_norm()
            #Extract horizontal velocity norm contained in t_slice
            vel = self._var.hori_velo_norm[argtime,:,:]

            #Sigma levels to consider
            if not top_lvl:
                top_lvl = (dep.shape[1]) - 1
            if not bot_lvl:
                bot_lvl = 0
            sLvl = range(bot_lvl, top_lvl+1)

            # Compute shear
            dz = dep[:,sLvl[1:],:] - dep[:,sLvl[:-1],:]
            dvel = vel[:,sLvl[1:],:] - vel[:,sLvl[:-1],:]           
            dveldz = dvel / dz

            #Custopm return
            self._var.verti_shear = dveldz 
            
            # Add metadata entry
            self._QC.append('vertical shear computed')
            print '-Vertical shear added to FVCOM.Variables.-'

        else:
            print "This function is only available for 3D FVCOM runs"

        if debug:
            print '...Passed'

    def velo_norm(self, debug=False):
        """Compute velocity norm -> FVCOM.Variables.velo_norm"""
        if debug or self._debug:
            print 'Computing velocity norm...'

        #TR_comment: not sure we should compute norm only in bounding box
        #u = self._var.u[:, :, self._grid.region_e[:]]
        #v = self._var.v[:, :, self._grid.region_e[:]]
        #ww = self._var.ww[:, :, self._grid.region_e[:]]
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

    def flow_dir_at_point(self, pt_lon, pt_lat,
                          vertical=False, debug=False):
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
          -vertical = True, compute flowDir for each vertical level
        Notes:
        -----
          directions between 0 and 360 deg.
        """
        debug = debug or self._debug
        if debug:
            print 'Computing flow directions at point...'

        #Choose the right pair of velocity components
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
        norm = ne.evaluate('sqrt(U**2 + V**2)')
        if debug:
            print '...Passed'

        return dirFlow, norm

    def flow_dir(self, debug=False):
        """"
        Compute flow directions over the whole grid -> FVCOM.Variables.flow_dir
        Notes:
        -----
          - directions between 0 and 360 deg.
          - This is very very SLOW !!!
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

