#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import sys
import numexpr as ne
from datetime import datetime
from datetime import timedelta
from interpolation_utils import closest_point, interp_at_point

class FunctionsFvcom:
    ''''Utils' subset of FVCOM class gathers useful functions""" '''
    def __init__(self, cls):
        self._debug = cls._debug
        self._var = cls.Variables
        self._grid = cls.Grid
        self._plot = cls.Plots
        self._QC = cls.QC
        #Create pointer to FVCOM class
        cls.Variables = self._var
        cls.Grid = self._grid
        cls.QC = self._QC
  
    def _centers(self, var, debug=False):
        '''Compute bathy and elevation at center points -> FVCOM.Grid.hc, elc'''
        if debug or self._debug:
            print 'Computing central bathy...'

        size = self.trinodes.T[elements].shape[0]
        size1 = self.elev.shape[0]
        elc = np.zeros((size1, size))
        hc = np.zeros((size))
        for ind, value in enumerate(self.trinodes.T[elements]):
            elc[:, ind] = np.mean(self.elev[:, value-1], axis=1)
            hc[ind] = np.mean(self.h[value-1])

        #Custom return    
        self._grid.hc = hc
        self._var.elc = elc
          
        # Add metadata entry
        self._QC.append('bathymetry at center points computed')
        self._QC.append('elevation at center points computed')
        print '-Central bathy and elevation added to FVCOM.Grid.-'

        if debug or self._debug:
            print '...Passed'   

    def _ele_region(self, debug=False):
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

    def _node_region(self, debug=False):
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

    def bounding_box(self, ax=[], quiet=False):
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
        self._node_region()
        self._ele_region()

        # Add metadata entry
        if not quiet:
            text = 'bounding box =' + str(self._grid.ax)
            self._QC.append(text)
            print '-Now working in bounding box-'      

    def hori_velo_norm(self, debug=False):
        """Compute horizontal velocity norm -> FVCOM.Variables.hori_velo_norm"""
        if debug or self._debug:
            print 'Computing horizontal velocity norm...'
        if self._var._3D:
            #TR_comment: not sure we should compute norm only in bounding box
            #u = self._var.u[:, :, self._grid.region_e[:]]
            #v = self._var.v[:, :, self._grid.region_e[:]]
            u = self._var.u[:, :, :]
            v = self._var.v[:, :, :]
            vel = ne.evaluate('sqrt(u**2 + v**2)')
        else:
            #TR_comment: not sure we should compute norm only in bounding box
            #u = self._var.ua[:, self._grid.region_e[:]]
            #v = self._var.va[:, self._grid.region_e[:]]
            u = self._var.ua[:, :]
            v = self._var.va[:, :]
            vel = ne.evaluate('sqrt(u**2 + v**2)')  

        #Custom return    
        self._var.hori_velo_norm = vel 
          
        # Add metadata entry
        self._QC.append('horizontal velocity norm computed')
        print '-Horizontal velocity norm added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    #TR: flow_dir over the whole grid far too slow !!!
    def flow_dir(self, debug=False):
        """"
        Compute flow directions over the whole grid -> FVCOM.Variables.flow_dir
        Notes:
        -----
          - directions between 0 and 360 deg.
          - This is very SLOW !!!
        """
        if debug or self._debug:
            print 'Computing flow directions...'

        u = self._var.ua
        v = self._var.va
        dirFlow = (np.pi/2.0) - np.arctan2(u,v)
        dirFlow = dirFlow * (180.0 / np.pi)

        #Custom return    
        self._var.dir_flow = dirFlow 

        # Add metadata entry
        self._QC.append('depth averaged flow directions computed')
        print '-Depth averaged flow directions added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def flow_dir_at_point(self, pt_lon, pt_lat, vertical=False, debug=False):
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

       u = self._var.ua
       v = self._var.va

        #Extraction at point
        if debug:
            print 'Extraction of u and v at point...'
        U = self.interpolation_at_point(u, pt_lon, pt_lat,
                                        debug=debug)  
        V = self.interpolation_at_point(v, pt_lon, pt_lat,
                                        debug=debug)       

        #Compute directions
        if debug:
            print 'Computing arctan2 and norm...'
        dirFlow = np.arctan2(U,V)
        dirFlow = np.mod(((np.pi/2.0) - dirFlow) * (180.0 / np.pi), 360.0)
        norm = ne.evaluate('sqrt(U**2 + V**2)')
        if debug:
            print '...Passed'

        #Rose diagram
        self._plot.rose_diagram(dirFlow, norm)

        return dirFlow, norm

    def ebb_flood_split(self, debug=False):
        """"
        Compute principal flow directions -> FVCOM.Variables.principal_axes
        Notes: directions between 0 and 360 deg.
        """
        if debug or self._debug:
            print 'Computing principal flow directions...'

        #TR : still to be defined

        # Add metadata entry
        self._QC.append('principal flow directions computed')
        print '-Principal Flow directions added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def principal_axes(self, lon, lat, dirFlow=[], debug=False):
        """"
        Compute principal flow directions and associated variances for (x, y) point
        Inputs:
        ------
          lon = longitude in deg., float
          lat = latitude in deg., float
        Outputs:
        -------
           pr_axes = 2 principal flow axes, 2 element list
           pr_ax_var = associated variances, 2 element list            
        Notes:
        -----
          directions between 0 and 360 deg.
        """
        if debug or self._debug:
            print 'Computing principal flow directions...'

        #TR: Still to be developed

        if debug or self._debug:
            print '...Passed'

    def interpolation_at_point(self, var, pt_lon, pt_lat, debug=False):
        """
        Interpol any given variables at any give location.
        Inputs:
        ------
          - var = variable, numpy array, dim=(time, nele or node)
          - pt_lon = longitude in degrees to find
          - pt_lat = latitude in degrees to find
        Outputs:
        -------
           - varInterp = var interpolate at (pt_lon, pt_lat)
        """
        debug = (debug or self._debug)
        if debug:
            print 'Interpoling at point...'
        lon = self._grid.lonc
        lat = self._grid.latc
        index = closest_point([pt_lon], [pt_lat], lon, lat, debug=debug)[0]

        lon = self._grid.lon
        lat = self._grid.lat    
        trinodes = self._grid.trinodes

        varInterp = interp_at_point(var, pt_lon, pt_lat, lon, lat, index,
                                    trinodes , debug=(self._debug or debug))
        return varInterp

    def exceedance(self, var, time, pt_lon=[], pt_lat=[], debug=False):
        """
        This function calculate the excedence curve of a var(time).
        Inputs:
        ------
          - time: time in seconds, 1D array of n elements
          - var: given quantity, 1 or 2D array of n elements, i.e (time) or (time,ele)
          - pt_lon, pt_lat: coordinates, necessary if var = 2D
        Outputs:
        -------
           - Exceedance: list of % of occurences
           - Ranges: list of signal amplitude bins
        Keywords:
        --------
           - graph: True->plots curve; False->does not
        Notes:
        -----
          This method is not suitable for SSE
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing exceedance...'

        #Distinguish between 1D and 2D var
        if len(var.shape)!=2:
            if not pt_lon or pt_lat:
                print 'Lon, lat coordinates are needed'
                sys.exit()  
            signal = self.interpolation_at_point(var, pt_lon, pt_lat, debug=debug)
        else:
            signal=var
        
        Max = max(signal)	
        dy = (Max/50.0)
        Ranges = np.arange(0,(Max + dy), dy)
        Exceedance = np.zeros(Ranges.shape[0])
        Period = time[-1] - time[0]

        N = len(signal)
        M = len(Ranges)

        for i in range(M):
            r = Ranges[i]
            for j in range(N-1):
                if signal[j] > r:
                    Exceedance[i] = Exceedance[i] + (time[j+1] - time[j])

        Exceedance = (Exceedance * 100) / Period

        if debug:
            print '...Passed'
       
        #Plot
        self._plot.plot_xy(Exceedance, Ranges, yLabel='Amplitudes',
                           xLabel='Exceedance probability in %')

        return Exceedance, Ranges
