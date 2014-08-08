#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from scipy import linalg as LA
import sys
import numexpr as ne
from datetime import datetime
from datetime import timedelta
from interpolation_utils import *
from miscellaneous import *
from BP_tools import *
import time

class FunctionsFvcom:
    """
    Description:
    -----------
    'Utils' subset of FVCOM class gathers
    useful functions for 2D runs
    """
    def __init__(self, variable, grid, plot, QC, debug):
        self._debug = debug
        self._var = variable
        self._grid = grid
        self._plot = plot
        self._QC = QC
        #Create pointer to FVCOM class
        variable = self._var
        grid = self._grid
        QC = self._QC

    #TR comment: I don't think I need this anymore  
    def _centers(self, var, debug=False):
        """
        Compute bathy and elevation at center points -> FVCOM.Grid.hc, elc

        Notes:
        -----
          - Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            print 'Computing central bathy...'

        #Interpolation at centers
        size = self._grid.nele
        size1 = self._grid.ntime
        elc = np.zeros((size1, size))
        hc = np.zeros((size))
        #TR comment: I am dubeous about the interpolation method here
        for ind, value in enumerate(self._grid.trinodes):
            elc[:, ind] = np.mean(self._var.el[:, value-1], axis=1)
            hc[ind] = np.mean(self._grid.h[value-1])

        #Custom return    
        self._grid.hc = hc
        self._var.elc = elc
          
        # Add metadata entry
        self._QC.append('bathymetry at center points computed')
        self._QC.append('elevation at center points computed')
        print '-Central bathy and elevation added to FVCOM.Grid.-'

        if debug:
            print '...Passed'   

    def hori_velo_norm(self, debug=False):
        """
        Compute new variable 'horizontal velocity norm' -> FVCOM.Variables.hori_velo_norm

        Notes:
        -----
          - Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            print 'Computing horizontal velocity norm...'

        try:
            u = self._var.ua[:, :]
            v = self._var.va[:, :]
            vel = ne.evaluate('sqrt(u**2 + v**2)')  
        except MemoryError:
            print '---Data too large for machine memory---'
            print 'Tip: use ax or tx during class initialisation'
            print '---  to use partial data'
            raise

        #Custom return    
        self._var.hori_velo_norm = vel 
          
        # Add metadata entry
        self._QC.append('horizontal velocity norm computed')
        print '-Horizontal velocity norm added to FVCOM.Variables.-'

        if debug:
            print '...Passed'

    def flow_dir(self, debug=False):
        """"
        Compute depth averaged flow directions over the whole grid
        -> FVCOM.Variables.depth_av_flow_dir

        Notes:
        -----
          - directions between 0 and 360 deg., i.e. 0/360=East, 90=North,
            180=West, 270=South
          - Can take time over the full domain
        """
        if debug or self._debug:
            print 'Computing flow directions...'

        try:
            u = self._var.ua
            v = self._var.va
            dirFlow = np.rad2deg(np.arctan2(V,U))
            #Adapt to Rose diagram
            #ind = np.where(dirFlow<0)[0]
            #dirFlow[ind] = 360.0 + dirFlow[ind]
            #TR: not quite sure here, seems to change from location to location
            #express principal axis in compass
            dirFlow = np.mod(90 - dirFlow, 360.0)
        except MemoryError:
            print '---Data too large for machine memory---'
            print 'Tip: use ax or tx during class initialisation'
            print '---  to use partial data'
            raise

        #Custom return    
        self._var.depth_av_flow_dir = dirFlow 

        # Add metadata entry
        self._QC.append('depth averaged flow directions computed')
        print '-Depth averaged flow directions added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def flow_dir_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[], time_ind=[],
                          exceedance=False, vertical=False, debug=False):
        """
        Flow directions and associated norm at any give location.

        Inputs:
        ------
          - pt_lon = longitude in decimal degrees East to find
          - pt_lat = latitude in decimal degrees North to find
        Outputs:
        -------
           - flowDir = flowDir at (pt_lon, pt_lat)
           - norm = velocity norm at (pt_lon, pt_lat)
        Keywords:
        --------
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - t_end = end time, as string ('yyyy-mm-ddThh:mm:ss'), or time index (integer)
          - time_ind = time indexes to work in, list of integers
          - excedance = True, compute associated exceedance curve
          - vertical = True, compute flowDir for each vertical level
        Notes:
        -----
          -directions between 0 and 360 deg., i.e. 0/360=East, 90=North,
           180=West, 270=South 
        """
        debug = debug or self._debug
        if debug:
            print 'Computing flow directions at point...'

        # Find time interval to work in
        argtime = []
        if not time_ind==[]:
            argtime = time_ind
        elif not t_start==[]:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                argtime = arange(t_start, t_end)

        #Choose the right pair of velocity components
        if not argtime==[]:
            u = self._var.ua[argtime,:]
            v = self._var.va[argtime,:]
        else:
            u = self._var.ua
            v = self._var.va

        #Extraction at point
        # Finding closest point
        index = closest_point([pt_lon], [pt_lat],
                              self._grid.lonc,
                              self._grid.latc, debug=debug)[0]
        if debug:
            print 'Extraction of u and v at point...'
        U = self.interpolation_at_point(u, pt_lon, pt_lat, index=index,
                                        debug=debug)  
        V = self.interpolation_at_point(v, pt_lon, pt_lat, index=index,
                                        debug=debug)       
 
        #Checking if dir_flow already computed
        if not hasattr(self._var, 'depth_av_flow_dir'):
            #Compute directions
            if debug:
                print 'Computing arctan2 and norm...'
            dirFlow = np.rad2deg(np.arctan2(V,U))#
            #Adapt to Rose diagram
            #TR: not quite sure here, seems to change from location to location
            dirFlow = np.mod(90.0 - dirFlow, 360.0)
        else:
            if not argtime==[]:
                dir_flow = self._var.dir_flow_hori[argtime,:,:]
                dirFlow = self.interpolation_at_point(self._var.depth_av_flow_dir,
                                                      pt_lon, pt_lat,
                                                      index=index, debug=debug)   
            else:
                dirFlow = self.interpolation_at_point(self._var.depth_av_flow_dir,
                                                      pt_lon, pt_lat,
                                                      index=index, debug=debug) 

        #Compute velocity norm
        norm = ne.evaluate('sqrt(U**2 + V**2)')
        if debug:
            print '...Passed'
        #Rose diagram
        self._plot.rose_diagram(dirFlow, norm)
        if exceedance:
            self.exceedance(norm)

        return dirFlow, norm

    def ebb_flood_split_at_point(self, pt_lon, pt_lat,
                                 t_start=[], t_end=[], time_ind=[], debug=False):
        """
        Compute time indexes for ebb and flood but also the 
        principal flow directions and associated variances for (lon, lat) point

        Inputs:
        ------
          - pt_lon = longitude in decimal degrees East to find
          - pt_lat = latitude in decimal degrees North to find
        Outputs:
        -------
          - floodIndex = flood time index, 1D array of integers
          - ebbIndex = ebb time index, 1D array of integers
          - pr_axis = principal flow ax1s, number in degrees
          - pr_ax_var = associated variance
        Keywords:
        --------
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'),
            or time index (integer)
          - t_end = end time, as string ('yyyy-mm-ddThh:mm:ss'),
            or time index(integer)
          - time_ind = time indexes to work in, 1D array of integers         
        Notes:
        -----
          - may take time to compute if time period too long
          - directions between 0 and 360 deg., i.e. 0/360=East, 90=North,
            180=West, 270=South
          - use time_ind or t_start and t_end, not both
          - assume that flood is aligned with principal direction
        """
        debug = debug or self._debug
        if debug:
            start = time.time()
            print 'Computing principal flow directions...'

        # Find time interval to work in
        argtime = []
        if not time_ind==[]:
            argtime = time_ind
        elif not t_start==[]:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end,
                                        self._var.matlabTime,
                                        debug=debug)
            else:
                argtime = arange(t_start, t_end)

        #Choose the right pair of velocity components
        if not argtime==[]:
            u = self._var.ua[argtime,:]
            v = self._var.va[argtime,:]
        else:
            u = self._var.ua
            v = self._var.va

        #Extraction at point
        # Finding closest point
        index = closest_point([pt_lon], [pt_lat],
                              self._grid.lonc,
                              self._grid.latc, debug=debug)[0]
        if debug:
            print 'Extraction of u and v at point...'
        U = self.interpolation_at_point(u, pt_lon, pt_lat, index=index,
                                        debug=debug)  
        V = self.interpolation_at_point(v, pt_lon, pt_lat, index=index,
                                        debug=debug)  

        #WB version of BP's principal axis
        if debug:
            print 'Computin principal axis at point...'
        pr_axis, pr_ax_var = principal_axis(U, V)

        #ebb/flood split
        if debug:
            print 'Splitting ebb and flood at point...'
        # reverse 0-360 deg convention
        ra = (-pr_axis - 90.0) * np.pi /180.0
        if ra>np.pi:
            ra = ra - (2.0*np.pi)
        elif ra<-np.pi:
            ra = ra + (2.0*np.pi)    
        dirFlow = np.arctan2(V,U)
        #Define bins of angles
        if ra == 0.0:
            binP = [0.0, np.pi/2.0]
            binP = [0.0, -np.pi/2.0]
        elif ra > 0.0:
            if ra == np.pi:
                binP = [np.pi/2.0 , np.pi]
                binM = [-np.pi, -np.pi/2.0 ]        
            elif ra < (np.pi/2.0):
                binP = [0.0, ra + (np.pi/2.0)]
                binM = [-((np.pi/2.0)-ra), 0.0]
            else:
                binP = [ra - (np.pi/2.0), np.pi]
                binM = [-np.pi, -np.pi + (ra-(np.pi/2.0))]
        else:
            if ra == -np.pi:
                binP = [np.pi/2.0 , np.pi]
                binM = [-np.pi, -np.pi/2.0]
            elif ra > -(np.pi/2.0):
                binP = [0.0, ra + (np.pi/2.0)]
                binM = [ ((-np.pi/2.0)+ra), 0.0]
            else:
                binP = [np.pi - (ra+(np.pi/2.0)) , np.pi]
                binM = [-np.pi, ra + (np.pi/2.0)]
                                
        test = (((dirFlow > binP[0]) * (dirFlow < binP[1])) +
                ((dirFlow > binM[0]) * (dirFlow < binM[1])))
        floodIndex = np.where(test == True)[0]
        ebbIndex = np.where(test == False)[0]

        #TR fit with Rose diagram angle convention
        pr_axis = pr_axis - 90.0
        if pr_axis<0.0:
            pr_axis[ind] = pr_axis[ind] + 360   

        if debug:
            end = time.time()
            print "...processing time: ", (end - start)

        return floodIndex, ebbIndex, pr_axis, pr_ax_var

    def interpolation_at_point(self, var, pt_lon, pt_lat, index=[], debug=False):
        """
        Interpol any given variables at any give location.

        Inputs:
        ------
          - var = variable, numpy array, dim=(time, nele or node)
          - pt_lon = longitude in decimal degrees East to find
          - pt_lat = latitude in decimal degrees North to find
        Outputs:
        -------
           - varInterp = var interpolate at (pt_lon, pt_lat)
        Keywords:
        --------
          - index = element index, integer
        Notes:
        -----
          - use index if closest element already known
        """
        debug = (debug or self._debug)
        if debug:
            print 'Interpolaling at point...'
        lonc = self._grid.lonc
        latc = self._grid.latc
        xc = self._grid.xc
        yc = self._grid.yc
        lon = self._grid.lon
        lat = self._grid.lat
        trinodes = self._grid.trinodes

        if index==[]:
            # Find indexes of the closest element
            index = closest_point([pt_lon], [pt_lat], lonc, latc, debug=debug)[0]
        # Conversion (lon, lat) to (x, y)
        pt_x = interp_at_point(self._grid.x, pt_lon, pt_lat, lon, lat,
                               index=index, trinodes=trinodes, debug=debug)
        pt_y = interp_at_point(self._grid.y, pt_lon, pt_lat, lon, lat,
                               index=index, trinodes=trinodes, debug=debug)

        #change in function of the data you dealing with
        if any(i == self._grid.node for i in var.shape):
            if debug:
                start = time.time() 
            varInterp = interpN_at_pt(var, pt_x, pt_y, xc, yc, index, trinodes,
                                      self._grid.aw0, self._grid.awx, self._grid.awy,
                                      debug=debug)
            if debug:
                end = time.time()
                print "Processing time: ", (end - start) 
        else:
            triele = self._grid.triele
            if debug:
                start = time.time() 
            varInterp = interpE_at_pt(var, pt_x, pt_y, xc, yc, index, triele,
                                      trinodes, self._grid.a1u, self._grid.a2u,
                                      debug=debug)
            if debug:
                end = time.time()
                print "Processing time: ", (end - start)         

        return varInterp

    def exceedance(self, var, pt_lon=[], pt_lat=[], debug=False):
        """
        This function calculate the excedence curve of a var(time).

        Inputs:
        ------
          - var = given quantity, 1 or 2D array of n elements, i.e (time) or (time,ele)
          - pt_lon, pt_lat = coordinates, necessary if var = 2D
        Outputs:
        -------
          - Exceedance = list of % of occurences
          - Ranges = list of signal amplitude bins
        Keywords:
        --------
          - graph: True->plots curve; False->does not
        Notes:
        -----
          - This method is not suitable for SSE
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing exceedance...'

        #Distinguish between 1D and 2D var
        if len(var.shape)>1:
            if not pt_lon or pt_lat:
                print 'Lon, lat coordinates are needed'
                sys.exit()
            signal = self.interpolation_at_point(var, pt_lon, pt_lat, debug=debug)
        else:
            signal=var
        
        Max = max(signal)	
        dy = (Max/30.0)
        Ranges = np.arange(0,(Max + dy), dy)
        Exceedance = np.zeros(Ranges.shape[0])
        dt = self._var.julianTime[1] - self._var.julianTime[0]
        Period = var.shape[0] * dt
        time = np.arange(0.0, Period, dt)

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

    def vorticity(self, debug=False):
        """
        Create a new variable 'depth averaged vorticity'
        -> FVCOM.Variables.depth_av_vorticity
     
        Notes:
        -----
          - Can take time over the full domain
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing vorticity...'
            start = time.time()

        t = arange(self._grid.ntime)  

        #Surrounding elements
        n1 = self._grid.triele[:,0]
        n2 = self._grid.triele[:,1]
        n3 = self._grid.triele[:,2]
        #TR comment: not quiet sure what this step does
        n1[np.where(n1==0)[0]] = self._grid.trinodes.shape[1]
        n2[np.where(n2==0)[0]] = self._grid.trinodes.shape[1]
        n3[np.where(n3==0)[0]] = self._grid.trinodes.shape[1]
        if debug:
            end = time.time()
            print "Check element=0, computation time in (s): ", (end - start)
            print "start np.multiply" 

        x0 = self._grid.xc
        y0 = self._grid.yc
        
        dvdx = np.zeros((self._grid.ntime,self._grid.nele))
        dudy = np.zeros((self._grid.ntime,self._grid.nele))

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
        self._var.depth_av_vorticity = vort
        self._History.append('depth averaged vorticity computed')
        print '-Depth averaged vorticity added to FVCOM.Variables.-'

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start) 

    def vorticity_over_period(self, time_ind=[], t_start=[], t_end=[], debug=False):
        """
        Compute the depth averaged vorticity for a time period.
     
        Outputs:
        -------
          - vort = horizontal vorticity (1/s), 2D array (time, nele)
        Keywords:
        -------
          - time_ind = time indexes to work in, list of integers
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'),
            or time index (integer)
          - t_end = end time, as string ('yyyy-mm-ddThh:mm:ss'),
            or time index (integer)
        Notes:
        -----
          - Can take time over the full domain
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing vorticity...'
            start = time.time()

        # Find time interval to work in
        t = []
        if not time_ind==[]:
            t = time_ind
        elif not t_start==[]:
            if type(t_start)==str:
                t = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                t = arange(t_start, t_end)
        else:
            t = arange(self._grid.ntime)  

        #Checking if vorticity already computed
        if not hasattr(self._var, 'depth_av_vorticity'): 
            #Surrounding elements
            n1 = self._grid.triele[:,0]
            n2 = self._grid.triele[:,1]
            n3 = self._grid.triele[:,2]
            #TR comment: not quiet sure what this step does
            n1[np.where(n1==0)[0]] = self._grid.trinodes.shape[1]
            n2[np.where(n2==0)[0]] = self._grid.trinodes.shape[1]
            n3[np.where(n3==0)[0]] = self._grid.trinodes.shape[1]
            if debug:
                end = time.time()
                print "Check element=0, computation time in (s): ", (end - start)
                print "start np.multiply" 

            x0 = self._grid.xc
            y0 = self._grid.yc
        
            dvdx = np.zeros((t.shape[0],self._grid.nele))
            dudy = np.zeros((t.shape[0],self._grid.nele))

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
        else:
            vort = self._var.depth_av_vorticity[t[:], :]

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start) 
        return vort

