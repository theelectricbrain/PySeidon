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
    'Util2D' subset of FVCOM class gathers
    useful functions for 2D and 3D runs
    """
    def __init__(self, variable, grid, plot, History, debug):
        self._debug = debug
        self._var = variable
        self._grid = grid
        self._plot = plot
        self._History = History
        #Create pointer to FVCOM class
        variable = self._var
        grid = self._grid
        History = self._History

    #TR comment: I don't think I need this anymore  
    def _centers(self, var, debug=False):
        """
        Create new variable 'bathy and elevation at center points' (m)
        -> FVCOM.Grid.hc, elc

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
        self._History.append('bathymetry at center points computed')
        self._History.append('elevation at center points computed')
        print '-Central bathy and elevation added to FVCOM.Grid.-'

        if debug:
            print '...Passed'   

    def hori_velo_norm(self, debug=False):
        """
        Compute new variable 'horizontal velocity norm' (m/s)
        -> FVCOM.Variables.hori_velo_norm

        Notes:
        -----
          - Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            print 'Computing horizontal velocity norm...'

        try:
            u = self._var.ua
            v = self._var.va
            vel = ne.evaluate('sqrt(u**2 + v**2)')  
        except MemoryError:
            print '---Data too large for machine memory---'
            print 'Tip: use ax or tx during class initialisation'
            print '---  to use partial data'
            raise

        #Custom return    
        self._var.hori_velo_norm = vel 
          
        # Add metadata entry
        self._History.append('horizontal velocity norm computed')
        print '-Horizontal velocity norm added to FVCOM.Variables.-'

        if debug:
            print '...Passed'

    def flow_dir(self, debug=False):
        """"
        Create new variable 'depth averaged flow directions' (deg.)
        -> FVCOM.Variables.depth_av_flow_dir

        Notes:
        -----
          - directions between -180 and 180 deg., i.e. 0=East, 90=North,
            +/-180=West, -90=South
          - Can take time over the full domain
        """
        if debug or self._debug:
            print 'Computing flow directions...'

        try:
            u = self._var.ua
            v = self._var.va
            dirFlow = np.rad2deg(np.arctan2(v,u))

        except MemoryError:
            print '---Data too large for machine memory---'
            print 'Tip: use ax or tx during class initialisation'
            print '---  to use partial data'
            raise

        #Custom return    
        self._var.depth_av_flow_dir = dirFlow 

        # Add metadata entry
        self._History.append('depth averaged flow directions computed')
        print '-Depth averaged flow directions added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def flow_dir_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[], time_ind=[],
                          exceedance=False, debug=False):
        """
        Flow directions and associated norm at any give location.

        Inputs:
        ------
          - pt_lon = longitude in decimal degrees East to find, float number 
          - pt_lat = latitude in decimal degrees North to find, float number 

        Outputs:
        -------
           - flowDir = flowDir at (pt_lon, pt_lat), 1D array
           - norm = velocity norm at (pt_lon, pt_lat), 1D array

        Keywords:
        --------
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'),
                      or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer
          - time_ind = time indices to work in, list of integers
          - excedance = True, compute associated exceedance curve

        Notes:
        -----
          - directions between -180 and 180 deg., i.e. 0=East, 90=North,
            +/-180=West, -90=South
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

        else:
            dirFlow = self.interpolation_at_point(self._var.depth_av_flow_dir,
                                                  pt_lon, pt_lat,
                                                  index=index, debug=debug) 

        #Compute velocity norm
        norm = ne.evaluate('sqrt(U**2 + V**2)')

        #use only the time indices of interest
        if not argtime==[]:
            dirFlow = dirFlow[argtime[:]]
            norm = norm[argtime[:]] 

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
        Compute time indices for ebb and flood but also the 
        principal flow directions and associated variances for (lon, lat) point

        Inputs:
        ------
          - pt_lon = longitude in decimal degrees East to find, float number 
          - pt_lat = latitude in decimal degrees North to find,float number 

        Outputs:
        -------
          - floodIndex = flood time index, 1D array of integers
          - ebbIndex = ebb time index, 1D array of integers
          - pr_axis = principal flow ax1s, float number in degrees
          - pr_ax_var = associated variance, float number

        Keywords:
        --------
          - t_start = start time, as a string ('yyyy-mm-ddThh:mm:ss'),
                      or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer
          - time_ind = time indices to work in, 1D array of integers 
        
        Notes:
        -----
          - may take time to compute if time period too long
          - directions between -180 and 180 deg., i.e. 0=East, 90=North,
            +/-180=West, -90=South
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

        #use only the time indices of interest
        if not argtime==[]:
            U = U[argtime[:]]
            V = V[argtime[:]] 

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
        #pr_axis = pr_axis - 90.0
        #if pr_axis<0.0:
        #    pr_axis[ind] = pr_axis[ind] + 360   

        if debug:
            end = time.time()
            print "...processing time: ", (end - start)

        return floodIndex, ebbIndex, pr_axis, pr_ax_var

    def interpolation_at_point(self, var, pt_lon, pt_lat, index=[], debug=False):
        """
        Interpol any given variables at any give location.

        Inputs:
        ------
          - var = any FVCOM grid data or variable, numpy array
          - pt_lon = longitude in decimal degrees East to find, float number
          - pt_lat = latitude in decimal degrees North to find, float number

        Outputs:
        -------
           - varInterp = var interpolated at (pt_lon, pt_lat)

        Keywords:
        --------
          - index = element index, integer. Use only if closest element index
                    is already known

        Notes:
        -----
          - use index if closest element already known
        """
        debug = (debug or self._debug)
        if debug:
            print 'Interpolaling at point...'
        lonc = self._grid.lonc[:]
        latc = self._grid.latc[:]
        xc = self._grid.xc[:]
        yc = self._grid.yc[:]
        lon = self._grid.lon[:]
        lat = self._grid.lat[:]
        trinodes = self._grid.trinodes[:]

        if index==[]:
            # Find indices of the closest element
            index = closest_point([pt_lon], [pt_lat], lonc, latc, debug=debug)[0]
        # Conversion (lon, lat) to (x, y)
        pt_x = interp_at_point(self._grid.x, pt_lon, pt_lat, lon, lat,
                               index=index, trinodes=trinodes, debug=debug)
        pt_y = interp_at_point(self._grid.y, pt_lon, pt_lat, lon, lat,
                               index=index, trinodes=trinodes, debug=debug)
        #change in function of the data you dealing with
        if any(i == self._grid.nnode for i in var.shape):
            if debug:
                start = time.time() 
            varInterp = interpN_at_pt(var, pt_x, pt_y, xc, yc, index, trinodes,
                                      self._grid.aw0, self._grid.awx, self._grid.awy,
                                      debug=debug)
            if debug:
                end = time.time()
                print "Processing time: ", (end - start) 
        else:
            triele = self._grid.triele[:]
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

        Keywords:
        --------
          - pt_lon, pt_lat = coordinates, float numbers.
                             Necessary if var = 2D (i.e. [time, nnode or nele]
          - graph: True->plots curve; False->does not

        Outputs:
        -------
          - Exceedance = list of % of occurences, 1D array
          - Ranges = list of signal amplitude bins, 1D array

        Notes:
        -----
          - This method is not suitable for SSE
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing exceedance...'

        #Distinguish between 1D and 2D var
        if len(var.shape)>1:
            if pt_lon==[] or pt_lat==[]:
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
        Create a new variable 'depth averaged vorticity (1/s)'
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
        #TR quick fix: due to error with pydap.proxy.ArrayProxy
        #              not able to cop with numpy.int
        n1 = int(n1)
        n2 = int(n2)
        n3 = int(n3)

        if debug:
            end = time.time()
            print "Check element=0, computation time in (s): ", (end - start)
            print "start np.multiply" 
        
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
          - time_ind = time indices to work in, list of integers
          - t_start = start time, as a string ('yyyy-mm-ddThh:mm:ss'),
                     or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer
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
            self.vorticity() 

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
            #TR quick fix: due to error with pydap.proxy.ArrayProxy
            #              not able to cop with numpy.int
            n1 = int(n1)
            n2 = int(n2)
            n3 = int(n3)

            if debug:
                end = time.time()
                print "Check element=0, computation time in (s): ", (end - start)
                print "start np.multiply" 
        
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

    def depth(self, debug=False):
        """
        Create new grid variable 'depth' (m)
        -> FVCOM.Grid.depth

        Notes:
        -----
          - depth convention: 0 = free surface
          - Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            start = time.time()

        print "Computing depth..."
        #Compute depth      
        size = self._grid.nele
        size1 = self._grid.ntime
        size2 = self._grid.nlevel
        elc = np.zeros((size1, size))
        hc = np.zeros((size))

        try:
            for ind, value in enumerate(self._grid.trinodes):
                elc[:, ind] = np.mean(self._var.el[:, value], axis=1)
                hc[ind] = np.mean(self._grid.h[value])

            dep = self._var.el[:,:] + h[None,:]
        except MemoryError:
            print '---Data too large for machine memory---'
            print 'Tip: use ax or tx during class initialisation'
            print '---  to use partial data'
            raise

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        # Add metadata entry
        self._grid.depth2D = dep
        self._History.append('depth 2D computed')
        print '-Depth 2D added to FVCOM.Variables.-'

    def depth_at_point(self, pt_lon, pt_lat, index=[], debug=False):
        """
        Compute depth at given point

        Inputs:
        ------
          - pt_lon = longitude in decimal degrees East, float number
          - pt_lat = latitude in decimal degrees North, float number

        Outputs:
        -------
          - dep = depth, 2D array (ntime, nlevel)

        Keywords:
        --------
          - index = element index, interger

        Notes:
        -----
          - depth convention: 0 = free surface
          - index is used in case one knows already at which
            element depth is requested
        """
        debug = debug or self._debug
        if debug:
            print "Computing depth..."
            start = time.time()

        #Finding index
        if index==[]:      
            index = closest_point([pt_lon], [pt_lat],
                                  self._grid.lonc,
                                  self._grid.latc, debug=debug)[0]

        if not hasattr(self._grid, 'depth2D'):
            #Compute depth
            h = self.interpolation_at_point(self._grid.h, pt_lon, pt_lat,
                                            index=index, debug=debug)
            el = self.interpolation_at_point(self._var.el, pt_lon, pt_lat,
                                             index=index, debug=debug)

            dep = el + h
        else:
            dep = self.interpolation_at_point(self._grid.depth2D, pt_lon, pt_lat,
                                            index=index, debug=debug)
        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        return dep

    def depth_averaged_power_assessment(self, debug=False):
        """
        Create a new variable 'depth averaged power density' (W/m2)
        -> FVCOM.Variables.depth_av_power_density

        Description:
        -----------
        The power density (pd) is then calculated as follows:
            pd = 0.5*1025*(u**3)

        Notes:
        -----
          - This may take some time to compute depending on the size
            of the data set
        """
        debug = (debug or self._debug)
        if debug: print "Computing depth averaged power density..."

        if not hasattr(self._var, 'hori_velo_norm'):
            if debug: print "Computing hori velo norm..."
            self.hori_velo_norm(debug=debug)
        if debug: print "Computing powers of hori velo norm..."
        u = self._var.hori_velo_norm
        if debug: print "Computing pd..."
        pd = ne.evaluate('0.5*1025.0*(u**3)')

        # Add metadata entry
        self._var.depth_av_power_density = pd
        self._History.append('depth averaged power density computed')
        print '-Depth averaged power density to FVCOM.Variables.-' 

    def depth_averaged_power_assessment(self, cut_in=1.0, cut_out=4.5, tsr=4.3, 
                                        a4=0.002, a3=-0.03, a2=0.1,
                                        a1=-0.1, a0=0.8,
                                        b2=-0.02, b1=0.2, b0=-0.005, debug=False):
        """
        Create a new variable 'depth averaged power assessment' (W/m2)
        -> FVCOM.Variables.depth_av_power_assessment

        Description:
        -----------
        This function performs tidal turbine power assessment by accounting for
        cut-in and cut-out speed, power curve (pc):
            pc = a4*(u**4) + a3*(u**3) + a2*(u**2) + a1*u + a0
           (where u is the flow speed)

        and device controled power coefficient (dcpc):
            dcpc =  b2*(tsr**2) + b1*tsr + b0

        The power density (pd) is then calculated as follows:
            pd = pc*dcpc*(1/2)*1025*(u**3)

        Keywords:
        --------
          - cut_in = cut-in speed in m/s, float number
          - cut_out = cut-out speed in m/s, float number
          - tsr = tip speed ratio, float number
          - a4 = pc curve parameter, float number
          - a3 = pc curve parameter, float number       
          - a2 = pc curve parameter, float number
          - a1 = pc curve parameter, float number
          - a0 = pc curve parameter, float number    
          - b2 = dcpc curve parameter, float number
          - b1 = dcpc curve parameter, float number
          - b0 = dcpc curve parameter, float number

        Notes:
        -----
          - This may take some time to compute depending on the size
            of the data set
        """
        debug = (debug or self._debug)
        if debug: print "Computing depth averaged power density..."

        if not hasattr(self._var, 'hori_velo_norm'):
            if debug: print "Computing hori velo norm..."
            self.hori_velo_norm(debug=debug)
        if debug: print "Computing powers of hori velo norm..."
        u = self._var.hori_velo_norm
        if debug: print "Computing pc and dcpc..."
        pc = ne.evaluate('a4*(u**4) + a3*(u**3) + a2*(u**2) + a1*u + a0')
        dcpc = ne.evaluate('b2*(tsr**2) + b1*tsr + b0')
        if debug: print "Computing pd..."
        pd = ne.evaluate('pc*dcpc*0.5*1025.0*(u**3)')

        if debug: print "finding cut-in and out..."
        u = cut_in
        pcin = ne.evaluate('a4*(u**4) + a3*(u**3) + a2*(u**2) + a1*u + a0')
        dcpcin = ne.evaluate('b2*(tsr**2) + b1*tsr + b0')
        pdin = ne.evaluate('pcin*dcpcin*0.5*1025.0*(u**3)')
        #TR comment huge bottleneck here
        #ind = np.where(pd<pdin)[0]
        #if not ind.shape[0]==0:
        #    pd[ind] = 0.0
        for i in range(pd.shape[0]):
            for j in range(pd.shape[1]):
                if pd[i,j] < pdin:
                   pd[i,j] = 0.0 

        u = cut_out
        pcout = ne.evaluate('a4*(u**4) + a3*(u**3) + a2*(u**2) + a1*u + a0')
        dcpcout = ne.evaluate('b2*(tsr**2) + b1*tsr + b0')
        pdout = ne.evaluate('pcout*dcpcout*0.5*1025.0*(u**3)')
        #TR comment huge bottleneck here
        #ind = np.where(pd>pdout)[0]
        #if not ind.shape[0]==0:
        #    pd[ind] = pdout
        for i in range(pd.shape[0]):
            for j in range(pd.shape[1]):
                if pd[i,j] > pdout:
                   pd[i,j] = pdout     

        # Add metadata entry
        self._var.depth_av_power_assessment = pd
        self._History.append('depth averaged power assessment computed')
        print '-Depth averaged power assessment to FVCOM.Variables.-'        
        

                     
