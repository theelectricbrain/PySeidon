#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np
import numexpr as ne
from datetime import datetime
from datetime import timedelta

class FunctionsFvcom:
    ''''Utils' subset of FVCOM class gathers useful functions""" '''
    def __init__(self, cls):
        self._debug = cls._debug
        self._var = cls.Variables
        self._grid = cls.Grid
        self._QC = cls.QC
        #Create pointer to FVCOM class
        cls.Variables = self._var
        cls.Grid = self._grid
        cls.QC = self._QC
  
    def elc(self, debug=False):
        '''Compute elevation at center points -> FVCOM.Variables.elc'''
        if debug or self._debug:
            print 'Computing elc...'

        size = self._grid.trinodes.T.shape[0]
        size1 = self._var.el.shape[0]
        elc = np.zeros((size1, size))
        for i,v in enumerate(self._grid.trinodes.T):
            elc[:, i] = np.mean(self._var.el[:, v], axis=1)

        if debug or self._debug:
            print '...Passed'

        # Add metadata entry
        self._QC.append('central elevation computed')
        print '-Central elevation added to FVCOM.Variables.-'

        #Custom return 
        self._var.elc = elc


    def hc(self, debug=False):
        '''Compute bathymetry at center points -> FVCOM.Variables.elc'''
        if debug or self._debug:
            print 'Computing hc...'

        size = self._grid.trinodes.T.shape[0]
        size1 = self._var.el.shape[0]
        elc = np.zeros((size1, size))
        for i,v in enumerate(self._grid.trinodes.T):
            hc[i] = np.mean(self._grid.h[v], axis=1)

        if debug or self._debug:
            print '...Passed'

        # Add metadata entry
        self._QC.append('central bathy computed')
        print '-Central bathy added to FVCOM.Variables.-'

        #Custom return 
        self._var.hc = hc       

    def closest_point(self, pt_lon, pt_lat,debug=False):
        '''
        Finds the closest exact lon, lat centre indexes of an FVCOM class
        to given lon, lat coordinates.

        Inputs:
          - pt_lon = list of longitudes in degrees
          - pt_lat = list of latitudes in degrees
        Outputs:
          - closest_point_indexes = numpy array of grid indexes
        '''
        if debug or self._debug:
            print 'Computing closest_point_indexes...'

        points = np.array([pt_lon, pt_lat]).T

        #point_list = np.array([lon,lat]).T
        point_list = np.array([self._var.lonc, self._var.latc]).T

        closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                        (point_list[:, 1] - points[:, 1, None])**2)

        closest_point_indexes = np.argmin(closest_dist, axis=1)

        if debug or self._debug:
            print '...Passed'

        return closest_point_indexes

    def ele_region(self, debug=False):
        '''Return element indexes included in bounding box, aka ax'''
        
        if debug or self._debug:
            print 'Computing region_e...'

        region_e = np.argwhere((self._grid.lonc >= self._grid.ax[0]) &
                                     (self._grid.lonc <= self._grid.ax[1]) &
                                     (self._grid.latc >= self._grid.ax[2]) &
                                     (self._grid.latc <= self._grid.ax[3]))

        # Create new grid variable through pointer
        self._grid.region_e = region_e[:,0]

        if debug or self._debug:
            print '...Passed'

        return region_e

    def node_region(self, debug=False):
        '''Return node indexes included in bounding box, aka ax'''
        if debug or self._debug:
            print 'Computing region_n...'

        region_n = np.argwhere((self._grid.lon >= self._grid.ax[0]) &
                                     (self._grid.lon <= self._grid.ax[1]) &
                                     (self._grid.lat >= self._grid.ax[2]) &
                                     (self._grid.lat <= self._grid.ax[3]))

        # Create new grid var through pointer
        self._grid.region_n = region_n[:,0]

        if debug or self._debug:
            print '...Passed'

        return region_n

    def bounding_box(self, ax=[]):
        """
        Define bounding box and reset the box by default.
        Input ex:
        --------
          .bounding_box(ax=[min lon, max lon, min lat, max lat])
        """
        # reset through pointer
        if ax:
            self._grid.ax = ax
        else:
            self._grid.ax = [min(self._grid.lon), max(self._grid.lon),
                             min(self._grid.lat), max(self._grid.lat)]
        self.node_region()
        self.ele_region()

        # Add metadata entry
        text = 'bounding box =' + str(self._grid.ax)
        self._QC.append(text)
        print '-Now working in bounding box-'      

    def velo_norm(self, debug=False):
        """Compute velocity norm -> FVCOM.Variables.velo_norm"""
        if debug or self._debug:
            print 'Computing velocity norm...'
        if self._var._3D:
            #TR_comment: not sure we should compute norm only in bounding box
            #u = self._var.u[:, :, self._grid.region_e[:]]
            #v = self._var.v[:, :, self._grid.region_e[:]]
            #ww = self._var.ww[:, :, self._grid.region_e[:]]
            u = self._var.u[:, :, :]
            v = self._var.v[:, :, :]
            ww = self._var.ww[:, :, :]
            vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')
        else:
            #TR_comment: not sure we should compute norm only in bounding box
            #u = self._var.ua[:, self._grid.region_e[:]]
            #v = self._var.va[:, self._grid.region_e[:]]
            u = self._var.ua[:, :]
            v = self._var.va[:, :]
            vel = ne.evaluate('sqrt(u**2 + v**2)')  

        #Custom return    
        self._var.velo_norm = vel 
          
        # Add metadata entry
        self._QC.append('velocity norm computed')
        print '-Velocity norm added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    # Functions available only for 3D runs
    if self._var._3D:
        def verti_shear(self, t_start, t_end, bot_lvl=[], top_level= [], debug=False):
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
            """
            if debug or self._debug:
                print 'Computing vertical shear...'

            # Find simulation time contains in [t_start, t_end]
            time = self._var.matlabTime
            t = time.shape[0]
            l = []
            for i in range(t):
                date = datetime.fromordinal(int(time[i]))
                     + timedelta(days=time[i]%1)-timedelta(days=366)
                l.append(date)
            time = np.array(l,dtype='datetime64[us]')
            t_slice = [t_start, t_end]
            t_slice = np.array(t_slice,dtype='datetime64[us]')

            if t_slice.shape[0] != 1:
                argtime = np.argwhere((time>=t_slice[0])&(time<=t_slice[-1])).flatten()
            if debug or self._debug:
                print argtime
            
            #Compute depth
            h = self._grid.h
            zeta = self.Variables.el[argtime,:] + h[None,:]
            nv = self._grid.nv[:].T-1
            siglay = self._grid.siglay[:]
            z = zeta[:,None,:]*siglay[None,:,:]
            dep = np.zeros([argtime.shape[0],siglay.shape[0],nv.shape[0]])
            for i in range(z.shape[0]):
                for j in range(nv.shape[0]):
                    el = (z[i,:,nv[j,0]]+z[i,:,nv[j,1]]+z[i,:,nv[j,2]])/3
                    dep[i,:,j] = el

            # Checking if velocity norm already exists
            if not hasattr(self._var, 'velo_norm'):
                self.velo_norm()
            #Extract velocity norm contained in t_slice
            vel = self._var.velo_norm[argtime,:,:]

            #Sigma levels to consider
            if not top:
                top = (dep.shape[1]) - 1
            if not bot:
                bot = 0
            sLvl = range(bot, top+1)

            # Compute shear
            dz = dep[:,sLvl[1:],:] - dep[:,sLvl[:-1],:]
            dvel = vel[:,sLvl[1:],:] - vel[:,sLvl[:-1],:]           
            dveldz = dvel / dz

            #Custopm return
            self._var.verti_shear = dveldz 
            
            # Add metadata entry
            self._QC.append('vertical shear computed')
            print '-Vertical shear added to FVCOM.Variables.-'


            if debug or self._debug:
                print '...Passed'
            
