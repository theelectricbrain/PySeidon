#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import sys
import numexpr as ne
from datetime import datetime
from datetime import timedelta
from interpolation_utils import * #closest_point, interp_at_point, interpN_at_pt
import time

class FunctionsFvcom:
    """'Utils' subset of FVCOM class gathers useful functions"""
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
        dirFlow = np.arctan2(u,v)
        dirFlow = np.mod((dirFlow + np.pi) * (180.0 / np.pi), 360.0)

        #Custom return    
        self._var.dir_flow = dirFlow 

        # Add metadata entry
        self._QC.append('depth averaged flow directions computed')
        print '-Depth averaged flow directions added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def flow_dir_at_point(self, pt_lon, pt_lat, exceedance=False,
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
        --------
          - excedance = True, compute associated exceedance curve
          - vertical = True, compute flowDir for each vertical level
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
        #TR: not quite sure here, seems to change from location to location
        #dirFlow = np.mod(((np.pi/2.0) - dirFlow) * (180.0 / np.pi), 360.0)
        #dirFlow = ((np.pi/2.0) - dirFlow) * (180.0 / np.pi)
        #dirFlow = dirFlow * (180.0 / np.pi)
        #dirFlow = np.pi - dirFlow) * (180.0 / np.pi)
        dirFlow = np.mod((dirFlow + np.pi) * (180.0 / np.pi), 360.0)
        norm = ne.evaluate('sqrt(U**2 + V**2)')
        if debug:
            print '...Passed'

        #Rose diagram
        self._plot.rose_diagram(dirFlow, norm)
        if exceedance:
            self.exceedance(norm, self._var.julianTime)

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
        lonc = self._grid.lonc
        latc = self._grid.latc
        xc = self._grid.xc
        yc = self._grid.yc
        lon = self._grid.lon
        lat = self._grid.lat
        trinodes = self._grid.trinodes

        # Find indexes of the closest element
        index = closest_point([pt_lon], [pt_lat], lonc, latc, debug=debug)[0]
        # Conversion (lon, lat) to (x, y)
        pt_x = interp_at_point(self._grid.x, pt_lon, pt_lat, lon, lat, index,
                                                       trinodes , debug=debug)
        pt_y = interp_at_point(self._grid.y, pt_lon, pt_lat, lon, lat, index,
                                                       trinodes , debug=debug)


        #change in function of the data you dealing with
        if any(i == self._grid.node for i in var.shape):
            if debug:
                start = time.time() 
            varInterp = interpN_at_pt(var, pt_x, pt_y, xc, yc, index, trinodes,
                                      self._grid.aw0, self._grid.awx, self._grid.awy,
                                      debug=debug)
            if debug:
                end = time.time()
                print "Time elapsed for Richard's method: ", (end - start) 
            #if debug:
            #    start = time.time()      
            #TR_alternative: works only for h and el at the mo i.e. node variables
            #test = interp_at_point(var, pt_lon, pt_lat, lon, lat, index,
            #                       trinodes , debug=debug)
            #TR_comment: seems to give same accuracy and computational time
            #            than Richard's approach
            #if debug:
            #    end = time.time()
            #    print "Time elapsed for Wesley's method: ", (end - start)
        else:
            triele = self._grid.triele
            if debug:
                start = time.time() 
            varInterp = interpE_at_pt(var, pt_x, pt_y, xc, yc, index, triele,
                                      trinodes, self._grid.a1u, self._grid.a2u,
                                      debug=debug)
            if debug:
                end = time.time()
                print "Time elapsed for Richard's method: ", (end - start)         

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
        if len(var.shape)==2:
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

    def vorticity(self, time_index=[], debug=False):
        """
        Compute the depth averaged vorticity for a time period:
        Dva/Dx - Dua/Dy -> vort

        Keyword:
        ------
          -time_index: time indexes, 1D array
        Notes:
        -----
          - Very lon if run over every time step
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing vorticity...'
            start = time.time()
        #Define time_index
        if time_index:
           t = np.asarray(time_index)
        else:
           t = arange(self._var.va.shape[0])
           print "It will compute for every time step. This may take a while !"
        #Surrounding elements
        n1 = self._grid.triele[:,0]
        n2 = self._grid.triele[:,1]
        n3 = self._grid.triele[:,2]
        #TR comment: not quiet sure what this step does
        n1[np.where(n1==-1)[0]] = self._grid.trinodes.shape[1] + 1
        n2[np.where(n2==-1)[0]] = self._grid.trinodes.shape[1] + 1
        n3[np.where(n3==-1)[0]] = self._grid.trinodes.shape[1] + 1
        if debug:
            end = time.time()
            print "Check element=-1, computation time in (s): ", (end - start)
            print "start np.multiply" 

        x0 = self._grid.xc
        y0 = self._grid.yc
        
        dvdx = np.zeros((t.shape[0], self._var.va.shape[1]))
        dudy = np.zeros((t.shape[0],self._var.va.shape[1]))

        j=0
        for i in t:
            dvdx[j,:] = np.multiply(self._grid.a1u[0,:], self._var.va[i,:]) \
                      + np.multiply(self._grid.a1u[1,:], self._var.va[i,n1]) \
                      + np.multiply(self._grid.a1u[2,:], self._var.va[i,n2]) \
                      + np.multiply(self._grid.a1u[3,:], self._var.va[i,n3])
            dudy[j,:] = np.multiply(self._grid.a2u[0,:], self._var.ua[i,:]) \
                      + np.multiply(self._grid.a2u[1,:], self._var.ua[i,n1]) \
                      + np.multiply(self._grid.a2u[2,:], self._var.ua[i,n2]) \
                      + np.multiply(self._grid.a2u[3,:], self._var.ua[i,n3])
            j+=1
        if debug:
            print "loop number ", i

        vort = dvdx - dudy

        # Add metadata entry
        #self._var.depth_av_vorticity = vort
        #self._QC.append('depth averaged vorticity computed')
        #print '-Depth averaged vorticity added to FVCOM.Variables.-'
        return vort

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start) 
