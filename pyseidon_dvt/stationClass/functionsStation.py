#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import sys
import numexpr as ne
from pyseidon_dvt.utilities.miscellaneous import *
from pyseidon_dvt.utilities.BP_tools import *
from utide import solve, reconstruct
import time

# Custom error
from pyseidon_error import PyseidonError

class FunctionsStation:
    """
    **'Util2D' subset of Station class gathers useful functions for 2D and 3D runs**
    """
    def __init__(self, variable, grid, plot, History, debug):
        self._debug = debug
        self._plot = plot
        #Create pointer to Station class
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

    def flow_dir(self,  station, t_start=[], t_end=[], time_ind=[],
                 exceedance=False, debug=False):
        """
        Flow directions and associated norm at any give location.

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
           - flowDir = flowDir at station, 1D array
           - norm = velocity norm at station, 1D array

        Options:
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - time_ind = time indices to work in, list of integers
          - excedance = True, compute associated exceedance curve

        *Notes*
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
                argtime = np.arange(t_start, t_end)

        #Search for the station
        index = self.search_index(station)

        #Choose the right pair of velocity components
        if not argtime==[]:
            U = self._var.ua[argtime,index]
            V = self._var.va[argtime,index]
        else:
            U = self._var.ua[:,index]
            V = self._var.va[:,index]
   
        #Compute directions
        if debug:
            print 'Computing arctan2 and norm...'
        dirFlow = np.rad2deg(np.arctan2(V,U))

        #Compute velocity norm
        norm = ne.evaluate('sqrt(U**2 + V**2)').squeeze()
        if debug:
            print '...Passed'
        #Rose diagram
        self._plot.rose_diagram(dirFlow, norm)
        if exceedance:
            self.exceedance(norm)

        return dirFlow, norm

    def ebb_flood_split(self, station,
                        t_start=[], t_end=[], time_ind=[], debug=False):
        """
        Compute time indices for ebb and flood but also the 
        principal flow directions and associated variances for (lon, lat) point

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
          - floodIndex = flood time index, 1D array of integers
          - ebbIndex = ebb time index, 1D array of integers
          - pr_axis = principal flow ax1s, float number in degrees
          - pr_ax_var = associated variance, float number

        Options:
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - time_ind = time indices to work in, list of integers
        
        *Notes*
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
                argtime = np.arange(t_start, t_end)

        #Search for the station
        index = self.search_index(station)

        #Choose the right pair of velocity components
        if not argtime==[]:
            U = self._var.ua[argtime,index]
            V = self._var.va[argtime,index]
        else:
            U = self._var.ua[:,index]
            V = self._var.va[:,index]

        #WB version of BP's principal axis
        if debug: print 'Computin principal axis at point...'
        pr_axis, pr_ax_var = principal_axis(U, V)

        #ebb/flood split
        if debug: print 'Splitting ebb and flood at point...'
        ###TR: version 1
        ## reverse 0-360 deg convention
        #ra = (-pr_axis - 90.0) * np.pi /180.0
        #if ra>np.pi:
        #    ra = ra - (2.0*np.pi)
        #elif ra<-np.pi:
        #    ra = ra + (2.0*np.pi)    
        #dirFlow = np.arctan2(V,U)
        ##Define bins of angles
        #if ra == 0.0:
        #    binP = [0.0, np.pi/2.0]
        #    binP = [0.0, -np.pi/2.0]
        #elif ra > 0.0:
        #    if ra == np.pi:
        #        binP = [np.pi/2.0 , np.pi]
        #        binM = [-np.pi, -np.pi/2.0 ]        
        #    elif ra < (np.pi/2.0):
        #        binP = [0.0, ra + (np.pi/2.0)]
        #        binM = [-((np.pi/2.0)-ra), 0.0]
        #    else:
        #        binP = [ra - (np.pi/2.0), np.pi]
        #        binM = [-np.pi, -np.pi + (ra-(np.pi/2.0))]
        #else:
        #    if ra == -np.pi:
        #        binP = [np.pi/2.0 , np.pi]
        #        binM = [-np.pi, -np.pi/2.0]
        #    elif ra > -(np.pi/2.0):
        #        binP = [0.0, ra + (np.pi/2.0)]
        #        binM = [ ((-np.pi/2.0)+ra), 0.0]
        #    else:
        #        binP = [np.pi - (ra+(np.pi/2.0)) , np.pi]
        #        binM = [-np.pi, ra + (np.pi/2.0)]
        #                        
        #test = (((dirFlow > binP[0]) * (dirFlow < binP[1])) +
        #        ((dirFlow > binM[0]) * (dirFlow < binM[1])))
        #floodIndex = np.where(test == True)[0]
        #ebbIndex = np.where(test == False)[0]

        #Defines interval
        flood_heading = np.array([-90, 90]) + pr_axis
        dir_all = np.rad2deg(np.arctan2(V,U))
        ind = np.where(dir_all<0)
        dir_all[ind] = dir_all[ind] + 360     

        # sign speed - eliminating wrap-around
        dir_PA = dir_all - pr_axis
        dir_PA[dir_PA < -90] += 360
        dir_PA[dir_PA > 270] -= 360

        #general direction of flood passed as input argument
        floodIndex = np.where((dir_PA >= -90) & (dir_PA<90))
        ebbIndex = np.arange(dir_PA.shape[0])
        ebbIndex = np.delete(ebbIndex,floodIndex[:]) 

        if debug:
            end = time.time()
            print "...processing time: ", (end - start)

        return floodIndex, ebbIndex, pr_axis, pr_ax_var

    def exceedance(self, var, station=[], debug=False):
        """
        This function calculate the excedence curve of a var(time).

        Inputs:
          - var = given quantity, 1 or 2D array of n elements, i.e (time) or (time,ele)

        Options:
          - station = either station index (interger) or name (string)
                      Necessary if var = 2D (i.e. [time, nnode or nele]
          - graph: True->plots curve; False->does not

        Outputs:
          - Exceedance = list of % of occurences, 1D array
          - Ranges = list of signal amplitude bins, 1D array

        *Notes*
          - This method is not suitable for SSE
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing exceedance...'

        #Distinguish between 1D and 2D var
        if len(var.shape)>1:
            if station==[]:
                print 'Lon, lat coordinates are needed'
                sys.exit()
            #Search for the station
            index = self.search_index(station)
            signal = var[:,index] 
        else:
            signal=var
        
        Max = max(signal)	
        dy = (Max/30.0)
        Ranges = np.arange(0,(Max + dy), dy)
        Exceedance = np.zeros(Ranges.shape[0])
        dt = self._var.julianTime[1] - self._var.julianTime[0]
        if dt==0:
            dt = self._var.secondTime[1] - self._var.secondTime[0]
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

    def depth(self, station, debug=False):
        """
        Compute depth at given point

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
          - dep = depth, 2D array (ntime, nlevel)

        *Notes*
          - depth convention: 0 = free surface
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
        dep = el + h

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        return dep

    def speed_histogram(self, station, t_start=[], t_end=[], time_ind=[], debug=False):
        """
        This function plots the histogram of occurrences for the signed
        flow speed at any given point.

        Inputs:
          - station = either station index (interger) or name (string)

        Options:
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - time_ind = time indices to work in, list of integers
        
        *Notes*
          - use time_ind or t_start and t_end, not both
        """
        debug = debug or self._debug
        if debug:
            start = time.time()
            print 'Computing speed histogram...'

        pI, nI, pa, pav = self.ebb_flood_split(station,
                          t_start=t_start, t_end=t_end, time_ind=time_ind,
                          debug=debug)
        dirFlow, norm = self.flow_dir(station,
                          t_start=t_start, t_end=t_end, time_ind=time_ind,
                          exceedance=False, debug=debug)
        norm[nI] = -1.0 * norm[nI]

        #compute bins
        #minBound = norm.min()
        #maxBound = norm.max()
        #step = round((maxBound-minBound/51.0),1)
        #bins = np.arange(minBound,maxBound,step)

        #plot histogram
        self._plot.Histogram(norm,
                             xLabel='Signed flow speed (m/s)',
                             yLabel='Occurrences (%)')
   
        if debug:
            end = time.time()
            print "...processing time: ", (end - start)


    def Harmonic_analysis_at_point(self, station,
                                   time_ind=[], t_start=[], t_end=[],
                                   elevation=True, velocity=False,
                                   debug=False, **kwarg):
        """
        This function performs a harmonic analysis on the sea surface elevation
        time series or the velocity components timeseries.

        Inputs:
          - station = either station index (interger) or name (string)

        Outputs:
          - harmo = harmonic coefficients, dictionary

        Options:
          - t_start = start time, as string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'), or time index as an integer
          - time_ind = time indices to work in, list of integers
          - elevation=True means that 'solve' will be done for elevation.
          - velocity=True means that 'solve' will be done for velocity.

        Utide's options:
        Options are the same as for 'solve', which are shown below with
        their default values:
            conf_int=True; cnstit='auto'; notrend=0; prefilt=[]; nodsatlint=0;
            nodsatnone=0; gwchlint=0; gwchnone=0; infer=[]; inferaprx=0;
            rmin=1; method='cauchy'; tunrdn=1; linci=0; white=0; nrlzn=200;
            lsfrqosmp=1; nodiagn=0; diagnplots=0; diagnminsnr=2;
            ordercnstit=[]; runtimedisp='yyy'

        *Notes*
        For more detailed information about 'solve', please see
        https://github.com/wesleybowman/UTide

        """
        debug = (debug or self._debug)

        #Search for the station
        index = self.search_index(station)

        argtime = []
        if not time_ind==[]:
            argtime = time_ind
        elif not t_start==[]:
            if type(t_start)==str:
                argtime = time_to_index(t_start, t_end,
                                        self._var.matlabTime,
                                        debug=debug)
            else:
                argtime = np.arange(t_start, t_end)

        if velocity == elevation:
            raise PyseidonError("---Can only process either velocities or elevation. Change options---")
        
        if velocity:
            time = self._var.matlabTime[:]
            u = self._var.ua[:,index]
            v = self._var.va[:,index]

            if not argtime==[]:
                time = time[argtime[:]]
                u = u[argtime[:]]
                v = v[argtime[:]]

            lat = self._grid.lat[index]
            harmo = solve(time, u, v, lat, **kwarg)

        if elevation:
            time = self._var.matlabTime[:]
            el = self._var.el[:,index]

            if not argtime==[]:
                time = time[argtime[:]]
                el = el[argtime[:]]

            lat = self._grid.lat[index]
            harmo = solve(time, el, None, lat, **kwarg)
            #Write meta-data only if computed over all the elements

        return harmo

    def Harmonic_reconstruction(self, harmo, time_ind=slice(None), debug=False, **kwarg):
        """
        This function reconstructs the velocity components or the surface elevation
        from harmonic coefficients.
        Harmonic_reconstruction calls 'reconstruct'. This function assumes harmonics
        ('solve') has already been executed.

        Inputs:
          - Harmo = harmonic coefficient from harmo_analysis
          - time_ind = time indices to process, list of integers
        
        Output:
          - Reconstruct = reconstructed signal, dictionary

        Options:
        Options are the same as for 'reconstruct', which are shown below with
        their default values:
            cnstit = [], minsnr = 2, minpe = 0

        *Notes*
        For more detailed information about 'reconstruct', please see
        https://github.com/wesleybowman/UTide

        """
        debug = (debug or self._debug)
        time = self._var.matlabTime[time_ind]
        Reconstruct = reconstruct(time, harmo)


        return Reconstruct  
