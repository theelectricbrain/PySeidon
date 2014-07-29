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
import time

class FunctionsFvcom:
    """'Utils' subset of FVCOM class gathers useful functions"""
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
        if self._var._3D:
            u = self._var.u[:, :, :]
            v = self._var.v[:, :, :]
            vel = ne.evaluate('sqrt(u**2 + v**2)')
        else:
            u = self._var.ua[:, :]
            v = self._var.va[:, :]
            vel = ne.evaluate('sqrt(u**2 + v**2)')  

        #Custom return    
        self._var.hori_velo_norm = vel 
          
        # Add metadata entry
        self._QC.append('horizontal velocity norm computed')
        print '-Horizontal velocity norm added to FVCOM.Variables.-'

        if debug:
            print '...Passed'

    def flow_dir(self, debug=False):
        """"
        Compute flow directions over the whole grid -> FVCOM.Variables.flow_dir

        Notes:
        -----
          - directions between 0 and 360 deg., i.e. 0/360=East, 90=North,
            180=West, 270=South
          - Can take time over the full domain
        """
        if debug or self._debug:
            print 'Computing flow directions...'

        u = self._var.ua
        v = self._var.va
        dirFlow = np.rad2deg(np.arctan2(V,U))
        #Adapt to Rose diagram
        #ind = np.where(dirFlow<0)[0]
        #dirFlow[ind] = 360.0 + dirFlow[ind]
        #TR: not quite sure here, seems to change from location to location
        #express principal axis in compass
        dirFlow = np.mod(90 - dirFlow, 360.0)

        #Custom return    
        self._var.dir_flow_hori = dirFlow 

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
          - pt_lon = longitude in degrees to find
          - pt_lat = latitude in degrees to find
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
        if time_ind:
            argtime = time_ind
        elif t_start:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                argtime = arange(t_start, t_end)

        #Choose the right pair of velocity components
        if argtime:
            u = self._var.ua[argtime,:]
            v = self._var.va[argtime,:]
        else:
            u = self._var.ua
            v = self._var.va

        #Extraction at point
        if debug:
            print 'Extraction of u and v at point...'
        U = self.interpolation_at_point(u, pt_lon, pt_lat,
                                        debug=debug)  
        V = self.interpolation_at_point(v, pt_lon, pt_lat,
                                        debug=debug)       
 
        #Checking if dir_flow already computed
        if not hasattr(self._var, 'dir_flow_hori'):
            #Compute directions
            if debug:
                print 'Computing arctan2 and norm...'
            dirFlow = np.rad2deg(np.arctan2(V,U))#
            #Adapt to Rose diagram
            #TR: not quite sure here, seems to change from location to location
            dirFlow = np.mod(90.0 - dirFlow, 360.0)
        else:
            if argtime:
                dir_flow = self._var.dir_flow_hori[argtime,:,:]
                dirFlow = self._util.interpolation_at_point(dir_flow, pt_lon, pt_lat,
                                                            debug=debug)   
            else:
                dirFlow = self._util.interpolation_at_point(self._var.dir_flow_hori,
                                                            pt_lon, pt_lat, debug=debug) 

        #Compute velocity norm
        norm = ne.evaluate('sqrt(U**2 + V**2)')
        if debug:
            print '...Passed'

        #Rose diagram
        self._plot.rose_diagram(dirFlow, norm)
        if exceedance:
            self.exceedance(norm, time_ind=argtime)

        return dirFlow, norm

    def ebb_flood_split_at_point(self, pt_lon, pt_lat,
                                 t_start=[], t_end=[], time_ind=[], debug=False):
        """"
        Compute time indexes for ebb and flood but also the 
        principal flow directions and associated variances for (lon, lat) point

        Inputs:
        ------
          - lon = longitude in deg., float
          - lat = latitude in deg., float
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
        """
        debug = debug or self._debug
        if debug:
            print 'Computing principal flow directions...'

        # Find time interval to work in
        argtime = []
        if time_ind:
            argtime = time_ind
        elif t_start:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                argtime = arange(t_start, t_end)

        #Choose the right pair of velocity components
        if argtime:
            u = self._var.ua[argtime,:]
            v = self._var.va[argtime,:]
        else:
            u = self._var.ua
            v = self._var.va

        #Extraction at point
        if debug:
            print 'Extraction of u and v at point...'
        U = self.interpolation_at_point(u, pt_lon, pt_lat,
                                        debug=debug)  
        V = self.interpolation_at_point(v, pt_lon, pt_lat,
                                        debug=debug)  

        #Build deviation matrix
        coord = np.matrix([U,V])
        center = np.mean(coord,1)
        coord = coord - center

        #Pierre Poulain's approach
        #inertia = np.dot(coord.transpose(), coord)
        #if debug:
        #    start = time.time()
        #    print 'Computing eigen values and vectors...'
        #e_values, e_vectors = np.linalg.eig(inertia)
        #if debug:
        #    end = time.time()
        #    print "...processing time: ", (end - start)
        #
        #if debug:
        #    print "(unordered) eigen values:"
        #    print e_values
        #    print "(unordered) eigen vectors:"
        #    print e_vectors

        ## order eigen values (and eigen vectors)
        ## axis1 is the principal axis with the biggest eigen value (eval1)
        ## axis2 is the principal axis with the second biggest eigen value (eval2)
        #for i in xrange(len(e_values)):
        #    # find biggest eigen value
        #    if e_values[i] == max(e_values):
        #        eval1 = e_values[i]
        #        axis1 = e_vectors[:,i]
        #    # find smallest eigen value
        #    else:
        #        eval2 = e_values[i]
        #        axis2 = e_vectors[:,i]

        #pr_axis = np.rad2deg(np.arctan2(axis1[1], axis1[0]))
        #pr_axis = np.mod(-90 -pr_axis, 360.0)
        #pr_ax_var = eval1/np.sum(e_values)

        #TR comment: THis Brian Polagye's and Kristen thyng approach
        ##Covariance
        print "Computing covariance..."
        if debug:
            start = time.time() 
        R = np.cov(coord)    
        if debug:
            end = time.time()
            print "...processing time: ", (end - start)

        ##Eigen values and vectors
        print "Computing eigen values and vectors..."
        if debug:
            start = time.time()
        eigVal, eigVec = np.linalg.eig(R)
        if debug:
            end = time.time()
            print "...processing time: ", (end - start)
            print eigVal, eigVec

        #TR performance test with Scipy
        #if debug:
        #    start = time.time()
        #eigValS, eigVecS = LA.eig(R)
        #if debug:
        #    end = time.time()
        #    print "...processing time with Scipy: ", (end - start)
        #    print eigValS, eigVecS

        #BP comments: Sort eignvalues in descending order
        #BP comments: so that major axis is given by first eigenvecto
        eigValF = np.sort(eigVal)[::-1]
        ilambda = np.argsort(eigVal)[::-1]
        #TR performance test with Scipy
        #eigValFS = np.sort(eigValS)[::-1]
        #ilambdaS = np.argsort(eigValS)[::-1]
        #BP comments: Reconstruct the eigenvalue matrix
        eigVal = np.diag(eigValF)
        #TR performance test with Scipy
        #eigValS = np.diag(eigValFS)
        ##BP comments: reorder the eigenvectors
        eigVec = eigVec[ilambda,:]
        ra = np.arctan2(eigVec[0,1], eigVec[1,1])
        #TR performance test with Scipy
        #eigVecS = eigVec[ilambdaS,:]
        #raS = np.arctan2(eigVecS[0,1], eigVecS[1,1])        
        #express principal axis in compass
        pr_axis = np.mod(90.0 - np.rad2deg(ra), 360.0)
        pr_ax_var = (eigVal[0]/np.trace(eigVal))[0]
        #TR performance test with Scipy
        #pr_axisS = np.mod(90.0 - np.rad2deg(raS), 360.0)
        #pr_ax_varS = (eigValS[0]/np.trace(eigValS))[0]
        #print 'Results with Numpy: ', pr_axis, pr_ax_var
        #print 'Results with Scipy: ', pr_axisS, pr_ax_varS

        #Computing ebb and flood indexes:
        print "Computing flood and ebb axis..."
        if debug:
            start = time.time()
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

        if debug:
            end = time.time()
            print "...processing time: ", (end - start)

        return floodIndex, ebbIndex, pr_axis, pr_ax_var

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
            print 'Interpolaling at point...'
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
        #Different interpolation method if partial data
        if not hasattr(self._grid, '_region_e'):
            pt_x = interp_at_point(self._grid.x, pt_lon, pt_lat, lon, lat,
                                   index=index, trinodes=trinodes, debug=debug)
            pt_y = interp_at_point(self._grid.y, pt_lon, pt_lat, lon, lat,
                                   index=index, trinodes=trinodes, debug=debug)

        #change in function of the data you dealing with
        if any(i == self._grid.node for i in var.shape):
            #Different interpolation method if partial data
            if hasattr(self._grid, '_region_e'):
                if debug:
                    start = time.time() 
                varInterp = interp_at_point(var, pt_lon, pt_lat, lon, lat,
                                            tri=trinodes , debug=debug)
                if debug:
                    end = time.time()
                    print "Processing time: ", (end - start) 
            else:
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
            #Different interpolation method if partial data
            if hasattr(self._grid, '_region_e'):
                if debug:
                    start = time.time() 
                varInterp = interp_at_point(var, pt_lon, pt_lat, lonc, latc,
                                            tri=triele , debug=debug)
                if debug:
                    end = time.time()
                    print "Processing time: ", (end - start) 
            else:
                if debug:
                    start = time.time() 
                varInterp = interpE_at_pt(var, pt_x, pt_y, xc, yc, index, triele,
                                          trinodes, self._grid.a1u, self._grid.a2u,
                                          debug=debug)
                if debug:
                    end = time.time()
                    print "Processing time: ", (end - start)         

        return varInterp

    def exceedance(self, var, time_ind, pt_lon=[], pt_lat=[], debug=False):
        """
        This function calculate the excedence curve of a var(time).

        Inputs:
        ------
          - time_ind = time in seconds, 1D array of n elements
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
        dy = (Max/30.0)
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

    def vorticity(self, time_ind=[], t_start=[], t_end=[], debug=False):
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
        if time_ind:
            t = time_ind
        elif t_start:
            if type(t_start)==str:
                t = time_to_index(t_start, t_end, self._var.matlabTime, debug=debug)
            else:
                t = arange(t_start, t_end)
        else:
            t = arange(self._grid.ntime)  

        #Surrounding elements
        #TR comment: I am not sure this one would work with new regioning protocol
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
        #self._var.depth_av_vorticity = vort
        #self._QC.append('depth averaged vorticity computed')
        #print '-Depth averaged vorticity added to FVCOM.Variables.-'

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start) 
        return vort
