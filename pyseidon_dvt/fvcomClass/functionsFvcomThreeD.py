#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numexpr as ne
import datetime
from scipy.interpolate import interp1d
from pyseidon_dvt.utilities.interpolation_utils import *
from pyseidon_dvt.utilities.miscellaneous import *
from pyseidon_dvt.utilities.BP_tools import *
from pyseidon_dvt.utilities.shortest_element_path import *
import time
import matplotlib.pyplot as plt
from pydap.exceptions import ServerError

#TR comment: This all routine needs to be tested and debugged
class FunctionsFvcomThreeD:
    """
    **'Utils3D' subset of FVCOM class gathers useful methods and functions for 3D runs**
    """
    def __init__(self, variable, grid, plot, util, History, debug):
        #Inheritance
        self._debug = debug
        self._plot = plot
        self._util = util
        self.interpolation_at_point = self._util.interpolation_at_point
        self.index_finder = self._util.index_finder
        self.hori_velo_norm = self._util.hori_velo_norm

        #Create pointer to FVCOM class
        setattr(self, '_var', variable)
        setattr(self, '_grid', grid)
        setattr(self, '_History', History)

        return

    def depth(self, debug=False):
        """
        This method computes new grid variable: 'depth' (m)
        -> FVCOM.Grid.depth

        *Notes*
          - depth convention: 0 = free surface
          - Can take time over the full domain
        """
        debug = debug or self._debug
        if debug:
            start = time.time()
            print "Computing depth..."

        try:
            elc = interpN(self._var.el[:], self._grid.trinodes[:], self._grid.aw0[:], debug=debug)
            hc = interpN(self._grid.h[:], self._grid.trinodes[:], self._grid.aw0[:], debug=debug)
            siglay = interpN(self._grid.siglay[:], self._grid.trinodes[:], self._grid.aw0[:], debug=debug)
            zeta = elc[:,:] + hc[None,:]
            dep = zeta[:,None,:]*siglay[None,:,:]

        except MemoryError:
             print '---Data too large for machine memory---'
             print 'Tip: use ax or tx during class initialisation'
             print '---  to use partial data'
             raise

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        # TR: need to find vectorized alternative
        #Compute depth
        # size = self._grid.nele
        # size1 = self._grid.ntime
        # size2 = self._grid.nlevel
        #
        # elc = np.zeros((size1, size))
        # hc = np.zeros((size)
        # dep = np.zeros((size1, size2, size))
        # for ind in range(size):
        #     d = self.depth_at_point(self._grid.lonc[ind],self._grid.latc[ind],debug=debug)
        #     dep[:, :, ind] = d[:,:]

        # TR: does not work with netCDF4 lib
        # elc = np.zeros((size1, size))
        # hc = np.zeros((size))
        # siglay = np.zeros((size2, size))
        #
        # try:
        #     for ind, value in enumerate(self._grid.trinodes[:]):
        #         elc[:, ind] = np.mean(self._var.el[:, value], axis=1)
        #         hc[ind] = np.mean(self._grid.h[value])
        #         siglay[:,ind] = np.mean(self._grid.siglay[:,value],1)
        #
        #     #zeta = self._var.el[:,:] + self._grid.h[None,:]
        #     zeta = elc[:,:] + hc[None,:]
        #     dep = zeta[:,None,:]*siglay[None,:,:]

        # Add metadata entry
        setattr(self._grid, 'depth', dep)
        self._History.append('depth computed')
        print '-Depth added to FVCOM.Grid.-'

    def depth_at_point(self, pt_lon, pt_lat, index=[], debug=False):
        """
        This function computes depth at any given point.

        Inputs:
          - pt_lon = longitude in decimal degrees East, float number
          - pt_lat = latitude in decimal degrees North, float number

        Outputs:
          - dep = depth, 2D array (ntime, nlevel)

        Options:
          - index = element index, interger. Use only if closest element
                    index is already known

        *Notes*
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
            index = self.index_finder(pt_lon, pt_lat, debug=False)

        if not hasattr(self._grid, 'depth'):
            #Compute depth
            h = self.interpolation_at_point(self._grid.h, pt_lon, pt_lat,
                                            index=index, debug=debug)
            el = self.interpolation_at_point(self._var.el, pt_lon, pt_lat,
                                             index=index, debug=debug)
            siglay = self.interpolation_at_point(self._grid.siglay, pt_lon, pt_lat,
                                                 index=index, debug=debug)
            zeta = el + h
            dep = zeta[:,None]*siglay[None,:]
        else:
            dep = self.interpolation_at_point(self._grid.depth[:],
                                              pt_lon, pt_lat, index=index,
                                              debug=debug)          
        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start)

        return dep

    def interp_at_depth(self, var, depth, ind=[], debug=False):
        """
        This function interpolates any given FVCOM.Variables field
        onto a specified depth plan

        Inputs:
          - var = 3 dimensional (time, sigma level, element) variable, array
          - depth = interpolation depth (float in meters), if negative = from
                    sea surface downwards, if positive = from sea bottom upwards
        Options:
          - ind = array of closest indexes to depth, 2D array (ntime, nele)

        Output:
          - interpVar = 2 dimensional (time, element) variable, masked array
          - ind = array of closest indexes to depth, 2D array (ntime, nele)
        """
        debug = debug or self._debug
        if debug: print 'Interpolating at '+str(depth)+' meter depth...'

        #checking if depth field already calculated
        if not hasattr(self._grid, 'depth'):
            self.depth()
        Depth = self._grid.depth[:]#otherwise to slow with netcdf4 lib
        if depth > 0.0:  # Changing vertical axis convention
            for tt in range(Depth.shape[0]):
                for ii in range(Depth.shape[2]):
                    mini = np.min(np.squeeze(Depth[tt, :, ii]))
                    Depth[tt, :, ii] = Depth[tt, :, ii] - mini
        dep = Depth[:] - depth
        #Finding closest values to specified depth
        if ind==[]:
            if debug: print 'Finding closest indexes to depth...'
            #mask negative value
            dep = np.ma.masked_where(dep<0.0, dep)
            #find min argument in masked array
            ind = dep.argmin(axis=1)
            ind=ind.astype(float)
            #set to nan to shallow elements
            ind[ind==dep.shape[1]-1.0] = np.nan

            #ind = np.zeros((dep.shape[0],dep.shape[2]))
            #for i in range(dep.shape[0]):
            #    for k in range(dep.shape[2]):
            #        test = dep[i,:,k]
            #        if not test[test>0.0].shape==test.shape:
            #            ind[i,k] = test[test>0.0].argmin()
            #        else:
            #            ind[i,k] = np.nan
                    
        inddown = ind + 1

        if debug: print 'Computing weights...'
        ##weight matrix & interp
        #interpVar = np.ones((var.shape[0], var.shape[2]))*np.nan
        #for i in range(ind.shape[0]):
        #    for j in range(ind.shape[1]):
        #        iU = ind[i,j]
        #        iD = inddown[i,j]
        #        if not np.isnan(iU):
        #            iU = int(iU)
        #            iD = int(iD)      
        #            length = np.abs(self._grid.depth[i,iU,j]\
        #                          - self._grid.depth[i,iD,j])
        #            wU = np.abs(depth - self._grid.depth[i,iU,j])/length
        #            wD = np.abs(depth - self._grid.depth[i,iD,j])/length
        #            interpVar[i,j] = (wU * var[i,iU,j]) + (wD * var[i,iD,j])
        #        else:
        #            interpVar[i,j] = np.nan
        #if debug: print 'Computing nan mask...'
        #interpVar = np.ma.masked_array(interpVar,np.isnan(interpVar))

        ##Streamlining
        I = []; J = []; U = []; D = []
        for i in range(ind.shape[0]):
            for j in range(ind.shape[1]):
                iU = ind[i,j]
                iD = inddown[i,j]
                I.append(i)
                J.append(j)
                U.append(iU)
                D.append(iD)
        if debug: print 'Convert lists to arrays...'
        I=np.asarray(I); J=np.asarray(J); U=np.asarray(U); D=np.asarray(D)
        if debug: print 'Find nan indices...'
        nanI = np.ones(U.shape)
        nanU = np.where(np.isnan(U))
        nanD = np.where(np.isnan(D))
        nanI[nanU] = np.nan
        nanI[nanD] = np.nan
        if debug: print 'convert to integer...'
        U[nanU] = 0
        D[nanD] = 0
        I = I.astype(int)
        J = J.astype(int)
        U = U.astype(int)
        D = D.astype(int)

        if type(var).__name__=='Variable': #Fix for netcdf4 lib
            Var = var[:]
        else:
            Var = var
        # dUp = Depth[I,U,J]
        # varUp = Var[I,U,J]
        # dDo = Depth[I,D,J]
        # varDo = Var[I,D,J]
        # TR: append to speed up caching
        if debug: print 'Caching...'
        for i in I:
            dUp = Depth[i,U,J]
            varUp = Var[i,U,J]
            dDo = Depth[i,D,J]
            varDo = Var[i,D,J]
        if debug: print 'Compute weights...'
        lengths = np.abs(dUp - dDo)
        wU = np.abs(depth - dUp)/lengths
        wD = np.abs(depth - dDo)/lengths
        if debug: print 'interpolation...'
        interpVar = nanI * ((wU * varUp) + (wD * varDo))
        if debug: print 'reshaping...'
        interpVar = np.reshape(interpVar, (Var.shape[0], Var.shape[2]))
        if debug: print 'masking...'        
        interpVar = np.ma.masked_array(interpVar,np.isnan(interpVar))
        if debug: print '...Passed'

        return interpVar, ind

    def verti_shear(self, debug=False):
        """
        This method computes a new variable: 'vertical shear' (1/s)
        -> FVCOM.Variables.verti_shear

        *Notes*
          - Can take time over the full doma
        """
        debug = debug or self._debug
        if debug:
            print 'Computing vertical shear...'
              
        #Compute depth if necessary
        if not hasattr(self._grid, 'depth'):        
           depth = self.depth(debug=debug)
        depth = self._grid.depth[:]

        # Checking if horizontal velocity norm already exists
        if not hasattr(self._var, 'velo_norm'):
            self.velo_norm()
        vel = self._var.velo_norm[:]

        try:
            #Sigma levels to consider
            top_lvl = self._grid.nlevel - 1
            bot_lvl = 0
            sLvl = range(bot_lvl, top_lvl+1)

            # Compute shear
            dz = depth[:,sLvl[1:],:] - depth[:,sLvl[:-1],:]
            dvel = vel[:,sLvl[1:],:] - vel[:,sLvl[:-1],:]           
            dveldz = dvel / dz
        except MemoryError:
            print '---Data too large for machine memory---'
            print 'Tip: use ax or tx during class initialisation'
            print '---  to use partial data'
            raise

        #Custom return
        setattr(self._var, 'verti_shear', dveldz)
            
        # Add metadata entry
        self._History.append('vertical shear computed')
        print '-Vertical shear added to FVCOM.Variables.-'

        if debug:
            print '...Passed'

    def verti_shear_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[],  time_ind=[],
                             bot_lvl=[], top_lvl=[], graph=True, dump=False, debug=False):
        """
        This function computes vertical shear at any given point.

        Inputs:
          - pt_lon = longitude in decimal degrees East, float number
          - pt_lat = latitude in decimal degrees North, float number

        Outputs:
          - dveldz = vertical shear (1/s), 2D array (time, nlevel - 1)

        Options:
          - t_start = start time, as a string ('yyyy-mm-ddThh:mm:ss'),
                      or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer
          - time_ind = time indices to work in, list of integers
          - bot_lvl = index of the bottom level to consider, integer
          - top_lvl = index of the top level to consider, integer
          - graph = plots graph if True
          - dump = boolean, dump profile data in csv file

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

        # Finding closest point
        index = self.index_finder(pt_lon, pt_lat, debug=False)

        #Compute depth
        depth = self.depth_at_point(pt_lon, pt_lat, index=index, debug=debug)       

        #Sigma levels to consider
        if top_lvl==[]:
            top_lvl = self._grid.nlevel - 1
        if bot_lvl==[]:
            bot_lvl = 0
        sLvl = range(bot_lvl, top_lvl+1)


        # Checking if vertical shear already exists
        if not hasattr(self._var, 'verti_shear'):
            if type(self._var.u).__name__=='Variable': 
                u = self._var.u[:]
                v = self._var.v[:]
            else:
                u = self._var.u
                v = self._var.v

            #Extraction at point
            if debug:
                print 'Extraction of u and v at point...'
            U = self.interpolation_at_point(u, pt_lon, pt_lat,
                                            index=index, debug=debug)  
            V = self.interpolation_at_point(v, pt_lon, pt_lat,
                                            index=index, debug=debug)
            norm = ne.evaluate('sqrt(U**2 + V**2)').squeeze()

            # Compute shear
            dz = depth[:,sLvl[1:]] - depth[:,sLvl[:-1]]
            dvel = norm[:,sLvl[1:]] - norm[:,sLvl[:-1]]           
            dveldz = dvel / dz
        else:
            if type(self._var.verti_shear).__name__=='Variable':
                shear = self._var.verti_shear[:]
            else:
                shear = self._var.verti_shear
            dveldz = self.interpolation_at_point(self._var.verti_shear,
                                                 pt_lon, pt_lat,
                                                 index=index, debug=debug)

        if debug:
            print '...Passed'
        #use time indices of interest
        if not argtime==[]:
            dveldz = dveldz[argtime,:]
            depth = depth[argtime,:]

        #Plot mean values
        if graph:
            mean_depth = np.mean((depth[:,sLvl[1:]]
                       + depth[:,sLvl[:-1]]) / 2.0, 0)
            mean_dveldz = np.mean(dveldz,0)
            error = np.std(dveldz,axis=0)/2.0
            self._plot.plot_xy(mean_dveldz, mean_depth, xerror=error[:],
                               title='Shear profile ',
                               xLabel='Shear (1/s) ', yLabel='Depth (m) ',
                               dump=dump)

        return dveldz             

    def velo_norm(self, debug=False):
        """
        This method computes a new variable: 'velocity norm' (m/s)
        -> FVCOM.Variables.velo_norm

        *Notes*
          -Can take time over the full domain
        """
        if debug or self._debug:
            print 'Computing velocity norm...'
        #Check if w if there
        try:
            try:
                #Computing velocity norm
                u = self._var.u[:]
                v = self._var.v[:]
                w = self._var.w[:]
                vel = ne.evaluate('sqrt(u**2 + v**2 + w**2)').squeeze()
            except (MemoryError, ServerError) as e:
                print '---Data too large for machine memory or server---'
                print 'Tip: Save data on your machine first'
                print 'Tip: use ax or tx during class initialisation'
                print '---  to use partial data'
                raise
        except AttributeError:
            try:
                #Computing velocity norm
                u = self._var.u[:]
                v = self._var.v[:]
                vel = ne.evaluate('sqrt(u**2 + v**2)').squeeze()
            except (MemoryError, ServerError) as e:
                print '---Data too large for machine memory or server---'
                print 'Tip: Save data on your machine first'
                print 'Tip: use ax or tx during class initialisation'
                print '---  to use partial data'
                raise

        #Custom return    
        setattr(self._var, 'velo_norm', vel)
       
        # Add metadata entry
        self._History.append('Velocity norm computed')
        print '-Velocity norm added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def velo_norm_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[], time_ind=[],
                           graph=True, dump=False, debug=False):
        """
        This function computes the velocity norm at any given point.

        Inputs:
          - pt_lon = longitude in decimal degrees East, float number
          - pt_lat = latitude in decimal degrees North, float number

        Outputs:
          - velo_norm = velocity norm, 2D array (time, level)

        Options:
          - t_start = start time, as a string ('yyyy-mm-ddThh:mm:ss'),
                      or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer
          - time_ind = time indices to work in, list of integers
          - graph = boolean, plots or not veritcal profile
          - dump = boolean, dump profile data in csv file

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

        try:
            if not hasattr(self._var, 'velo_norm'):
                if type(self._var.u).__name__=='Variable': #Fix for netcdf4 lib            
                    u = self._var.u[:]
                    v = self._var.v[:]
                    w = self._var.w[:]
                else:
                    u = self._var.u
                    v = self._var.v
                    w = self._var.w
            else:
                vel = self._var.velo_norm
        except AttributeError:
            if not hasattr(self._var, 'velo_norm'): 
                if type(self._var.u).__name__=='Variable': #Fix for netcdf4 lib            
                    u = self._var.u[:]
                    v = self._var.v[:]
                else:
                    u = self._var.u
                    v = self._var.v
            else:
                if type(self._var.velo_norm).__name__=='Variable': #Fix for netcdf4 lib:
                    vel = self._var.velo_norm[:]
                else:
                    vel = self._var.velo_norm

        # Finding closest point
        index = self.index_finder(pt_lon, pt_lat, debug=False)

        #Computing horizontal velocity norm
        if debug:
            print 'Extraction of u, v and w at point...'
        if not hasattr(self._var, 'velo_norm'): 
            U = self.interpolation_at_point(u, pt_lon, pt_lat,
                                            index=index, debug=debug)  
            V = self.interpolation_at_point(v, pt_lon, pt_lat,
                                            index=index, debug=debug)
            if 'w' in locals():
                W = self.interpolation_at_point(w, pt_lon, pt_lat,
                                                index=index, debug=debug)
                velo_norm = ne.evaluate('sqrt(U**2 + V**2 + W**2)').squeeze()
            else:
                velo_norm = ne.evaluate('sqrt(U**2 + V**2)').squeeze()
        else:
            velo_norm = self.interpolation_at_point(vel, pt_lon, pt_lat,
                                                    index=index, debug=debug)
        if debug:
            print '...passed'

        #use only the time indices of interest
        if not argtime==[]:
            velo_norm = velo_norm[argtime[:],:,:]

        #Plot mean values
        if graph:
            depth = self.depth_at_point(pt_lon, pt_lat, index=index)
            mean_depth = np.mean(depth, 0)
            mean_vel = np.mean(velo_norm,0)
            error = np.std(velo_norm,axis=0)/2.0
            self._plot.plot_xy(mean_vel, mean_depth, xerror=error[:],
                               title='Flow speed vertical  ',
                               xLabel='Flow speed (1/s) ', yLabel='Depth (m) ',
                               dump=dump)

        return velo_norm 


    def flow_dir_at_point(self, pt_lon, pt_lat, t_start=[], t_end=[], time_ind=[], 
                          vertical=True, debug=False):
        """
        This function computes flow directions and associated norm
        at any given location.

        Inputs:
          - pt_lon = longitude in decimal degrees East to find
          - pt_lat = latitude in decimal degrees North to find

        Outputs:
          - flowDir = flowDir at (pt_lon, pt_lat), 2D array (ntime, nlevel)
          - norm = velocity norm at (pt_lon, pt_lat), 2D array (ntime, nlevel)

        Options:
          - t_start = start time, as a string ('yyyy-mm-ddThh:mm:ss'),
                      or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer
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

        # Finding closest point
        index = self.index_finder(pt_lon, pt_lat, debug=False)

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
        if self._var._3D and vertical:
            if type(self._var.u).__name__=='Variable': #Fix for netcdf4 lib
                u = self._var.u[:]
                v = self._var.v[:]
            else:
                u = self._var.u
                v = self._var.v
        else:
            if type(self._var.ua).__name__=='Variable': #Fix for netcdf4 lib:
                u = self._var.ua[:]
                v = self._var.va[:]
            else:
                u = self._var.ua
                v = self._var.va

        #Extraction at point
        if debug:
            print 'Extraction of u and v at point...'
        U = self._util.interpolation_at_point(u, pt_lon, pt_lat,
                                              index=index, debug=debug)  
        V = self._util.interpolation_at_point(v, pt_lon, pt_lat,
                                              index=index, debug=debug)

        #Compute velocity norm
        norm = ne.evaluate('sqrt(U**2 + V**2)').squeeze()

        #Compute directions
        if debug:
            print 'Computing arctan2...'
        dirFlow = np.rad2deg(np.arctan2(V,U))

        if debug: print '...Passed'
        #use only the time indices of interest
        if not argtime==[]:
            dirFlow = dirFlow[argtime[:],:]

        return dirFlow, norm

    def flow_dir(self, debug=False):
        """"
        This method computes a new variable: 'flow directions' (deg.)
        -> FVCOM.Variables.flow_dir

        *Notes*
          - directions between -180 and 180 deg., i.e. 0=East, 90=North,
            +/-180=West, -90=South
          - Can take time over the full domain
        """
        if debug or self._debug:
            print 'Computing flow directions...'

        try:
            u = self._var.u[:]
            v = self._var.v[:]
            dirFlow = np.rad2deg(np.arctan2(v,u))
        except (MemoryError, ServerError) as e:
            print '---Data too large for machine memory or server---'
            print 'Tip: Save data on your machine'
            print 'Tip: use ax or tx during class initialisation'
            print '---  to use partial data'
            raise

        #Custom return    
        setattr(self._var, 'flow_dir', dirFlow)

        # Add metadata entry
        self._History.append('flow directions computed')
        print '-Flow directions added to FVCOM.Variables.-'

        if debug or self._debug:
            print '...Passed'

    def vorticity(self, debug=False):
        """
        This method creates a new variable: 'depth averaged vorticity' (1/s)
        -> FVCOM.Variables.vorticity
     
        *Notes*
          - Can take time over the full domain
        """
        debug = (debug or self._debug)
        if debug:
            print 'Computing vorticity...'
            start = time.time()

        t = np.arange(self._grid.ntime)  

        #Surrounding elements
        n1 = self._grid.triele[:,0]
        n2 = self._grid.triele[:,1]
        n3 = self._grid.triele[:,2]
        #No need anymore
        ##change end bound indices 
        #test = self._grid.triele.shape[0]
        #n1[np.where(n1==test)[0]] = 0
        #n2[np.where(n2==test)[0]] = 0
        #n3[np.where(n3==test)[0]] = 0

        #TR quick fix: due to error with pydap.proxy.ArrayProxy
        #              not able to cop with numpy.int
        N1 = []
        N2 = []
        N3 = []

        N1[:] = n1[:]
        N2[:] = n2[:]
        N3[:] = n3[:]

        if debug:
            end = time.time()
            print "Check element=0, computation time in (s): ", (end - start)
            print "start np.multiply" 

        x0 = self._grid.xc
        y0 = self._grid.yc
        
        dvdx = np.zeros((self._grid.ntime,self._grid.nlevel,self._grid.nele))
        dudy = np.zeros((self._grid.ntime,self._grid.nlevel,self._grid.nele))
        nele = self._grid.nele

        j=0
        for i in t:
            dvdx[j,:,:] = np.multiply(self._grid.a1u[0,:].reshape(1,nele),
                          self._var.v[i,:,:]) \
                        + np.multiply(self._grid.a1u[1,:].reshape(nele,1),
                          self._var.v[i,:,N1]).T \
                        + np.multiply(self._grid.a1u[2,:].reshape(nele,1),
                          self._var.v[i,:,N2]).T \
                        + np.multiply(self._grid.a1u[3,:].reshape(nele,1),
                          self._var.v[i,:,N3]).T
            dudy[j,:,:] = np.multiply(self._grid.a2u[0,:].reshape(1,nele),
                          self._var.u[i,:,:]) \
                        + np.multiply(self._grid.a2u[1,:].reshape(nele,1),
                          self._var.u[i,:,N1]).T \
                        + np.multiply(self._grid.a2u[2,:].reshape(nele,1),
                          self._var.u[i,:,N2]).T \
                        + np.multiply(self._grid.a2u[3,:].reshape(nele,1),
                          self._var.u[i,:,N3]).T
            j+=1
        if debug:
            print "loop number ", i

        vort = dvdx - dudy

        # Add metadata entry
        setattr(self._var, 'vorticity', vort)
        self._History.append('vorticity computed')
        print '-Vorticity added to FVCOM.Variables.-'

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start) 

    def vorticity_over_period(self, time_ind=[], t_start=[], t_end=[], debug=False):
        """
        This function computes the vorticity for a time period.
     
        Outputs:
          - vort = horizontal vorticity (1/s), 2D array (time, nele)

        Options:
          - time_ind = time indices to work in, list of integers
          - t_start = start time, as a string ('yyyy-mm-ddThh:mm:ss'),
                      or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer
        *Notes*
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
                start = datetime.datetime.strptime(t_start, '%Y-%m-%d %H:%M:%S')
                end = datetime.datetime.strptime(t_end, '%Y-%m-%d %H:%M:%S')
                t = time_to_index(start, end, self._var.julianTime[:], debug=debug)
            else:
                t = np.arange(t_start, t_end)
        else:
            t = np.arange(self._grid.ntime)

        #Checking if vorticity already computed
        if not hasattr(self._var, 'vorticity'): 
            #Surrounding elements
            n1 = self._grid.triele[:,0]
            n2 = self._grid.triele[:,1]
            n3 = self._grid.triele[:,2]
            #No need anymore
            ##change end bound indices 
            #test = self._grid.triele.shape[0]
            #n1[np.where(n1==test)[0]] = 0
            #n2[np.where(n2==test)[0]] = 0
            #n3[np.where(n3==test)[0]] = 0
            #TR quick fix: due to error with pydap.proxy.ArrayProxy
            #              not able to cop with numpy.int
            N1 = []
            N2 = []
            N3 = []

            N1[:] = n1[:]
            N2[:] = n2[:]
            N3[:] = n3[:]

            if debug:
                end = time.time()
                print "Check element=0, computation time in (s): ", (end - start)
                print "start np.multiply" 

            x0 = self._grid.xc[:]
            y0 = self._grid.yc[:]
        
            dvdx = np.zeros((t.shape[0],self._grid.nlevel,self._grid.nele))
            dudy = np.zeros((t.shape[0],self._grid.nlevel,self._grid.nele))

            j=0
            for i in t:
                dvdx[j,:,:] = np.multiply(self._grid.a1u[0,:], self._var.v[i,:,:]) \
                          + np.multiply(self._grid.a1u[1,:], self._var.v[i,:,N1]) \
                          + np.multiply(self._grid.a1u[2,:], self._var.v[i,:,N2]) \
                          + np.multiply(self._grid.a1u[3,:], self._var.v[i,:,N3])
                dudy[j,:,:] = np.multiply(self._grid.a2u[0,:], self._var.u[i,:,:]) \
                          + np.multiply(self._grid.a2u[1,:], self._var.u[i,:,N1]) \
                          + np.multiply(self._grid.a2u[2,:], self._var.u[i,:,N2]) \
                          + np.multiply(self._grid.a2u[3,:], self._var.u[i,:,N3])
                j+=1
            if debug:
                print "loop number ", i

            vort = dvdx - dudy
        else:
            vort = self._var.vorticity[t[:],:,:]

        if debug:
            end = time.time()
            print "Computation time in (s): ", (end - start) 
        return vort

    def power_density(self, debug=False):
        """
        This method creates a new variable: 'power density' (W/m2)
        -> FVCOM.Variables.power_density

        The power density (pd) is then calculated as follows:
            pd = 0.5*1025*(u**3)

        *Notes*
          - This may take some time to compute depending on the size
            of the data set
        """
        debug = (debug or self._debug)
        if debug: print "Computing power density..."

        if not hasattr(self._var, 'velo_norm'):
            self.velo_norm(debug=debug)
        if debug: print "Computing power density variable..."
        u = self._var.velo_norm[:]
        pd = ne.evaluate('0.5*1025.0*(u**3)').squeeze()
        #pd = 0.5*1025.0*np.power(self._var.hori_velo_norm[:],3.0)  # TR: very slow
        #pd = 0.5*1025.0*self._var.hori_velo_norm[:]*self._var.hori_velo_norm[:]*self._var.hori_velo_norm[:]

        # Add metadata entry
        setattr(self._var, 'power_density', pd)
        self._History.append('power density computed')
        print '-Power density to FVCOM.Variables.-' 

    def power_assessment_at_depth(self, depth, power_mat, rated_speed, 
                                        cut_in=1.0, cut_out=4.5, debug=False):
        """
        This function computes power assessment (W/m2) at given depth.

        Description:
        This function performs tidal turbine power assessment by accounting for
        cut-in and cut-out speed, power curve/function (pc):
            Cp = pc(u)
           (where u is the flow speed)

        The power density (pd) is then calculated as follows:
            pd = Cp*(1/2)*1025*(u**3)

        Inputs:
          - depth = given depth from the surface, float
          - power_mat = power matrix (u,Cp(u)), 2D array (2,n),
                        u being power_mat[0,:] and Ct(u) being power_mat[1,:]
          - rated_speed = rated speed speed in m/s, float number


        Output:
          - pa = power assessment in (W/m2), 2D masked array (ntime, nele)

        Options:
          - cut_in = cut-in speed in m/s, float number
          - cut_out = cut-out speed in m/s, float number

        *Notes*
          - This may take some time to compute depending on the size
            of the data set
        """
        debug = (debug or self._debug)
        if debug: print "Computing depth averaged power density..."

        if not hasattr(self._var, 'power_density'):
            self.power_density(debug=debug)

        if debug: print "Initialising power curve..."
        Cp = interp1d(power_mat[0,:],power_mat[1,:])

        u, ind = self.interp_at_depth(self._var.velo_norm[:], depth, debug=debug)
        pd, ind2 = self.interp_at_depth(self._var.power_density[:], depth,
                                                   ind=ind, debug=debug)

        pa = Cp(u)*pd

        if debug: print "finding cut-in..."
        #TR comment huge bottleneck here
        #ind = np.where(pd<pdin)[0]
        #if not ind.shape[0]==0:
        #    pd[ind] = 0.0
        #for i in range(pa.shape[0]):
        #    for j in range(pa.shape[1]):
        #        if (u[i,j] < cut_in) or (u[i,j] > cut_out):
        #           pa[i,j] = 0.0
        
        inM = np.ma.masked_where(u<cut_in, u).mask
        outM = np.ma.masked_where(u>cut_out, u).mask
        ioM = inM * outM * u.mask
        pa=np.ma.masked_where(ioM, pa)

        if debug: print "finding rated speed..."
        parated = Cp(rated_speed)*0.5*1025.0*(rated_speed**3.0)
        #TR comment huge bottleneck here
        #ind = np.where(pd>pdout)[0]
        #if not ind.shape[0]==0:
        #    pd[ind] = pdout
        #for i in range(pa.shape[0]):
        #    for j in range(pa.shape[1]):
        #        if u[i,j] > rated_speed:
        #           pa[i,j] = parated
        pa[u>rated_speed] = parated     

        return pa 

    def _vertical_slice(self, var, start_pt, end_pt,
                        time_ind=[], t_start=[], t_end=[],
                        title='Title', cmax=[], cmin=[], debug=False):
        """
        Draw vertical slice in var along the shortest path between
        start_point, end_pt.
 
        Inputs:
          - var = 2D dimensional (sigma level, element) variable, array
          - start_pt = starting point, [longitude, latitude]
          - end_pt = ending point, [longitude, latitude]

        Options:
          - time_ind = reference time indices for surface elevation, list of integer
          - t_start = start time, as a string ('yyyy-mm-ddThh:mm:ss'),
                      or time index as an integer
          - t_end = end time, as a string ('yyyy-mm-ddThh:mm:ss'),
                    or time index as an integer

        Keywords for plot:
          - title = plot title, string
          - cmin = minimum limit colorbar
          - cmax = maximum limit colorbar
        """
        debug = debug or self._debug
        if not self._var._3D:
            print "Error: Only available for 3D runs."
            raise
        else: 
            lons = [start_pt[0], end_pt[0]]
            lats = [start_pt[1], end_pt[1]]
            #Finding the closest elements to start and end points
            index = closest_points(lons, lats,
                              self._grid.lonc,
                              self._grid.latc, debug=debug)
    
            #Finding the shortest path between start and end points
            if debug : print "Computing shortest path..."
            short_path = shortest_element_path(self._grid.lonc[:],
                                               self._grid.latc[:],
                                               self._grid.lon[:],
                                               self._grid.lat[:],
                                               self._grid.trinodes[:],
                                               self._grid.h[:], debug=debug)
            el, _ = short_path.getTargets([index])
            # Plot shortest path
            short_path.graphGrid(plot=True)

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
 
            #Extract along line
            ele=np.asarray(el[:])[0,:]
            varP = var[:,ele]
            # Depth along line
            if debug : print "Computing depth..."
            depth = np.zeros((self._grid.ntime, self._grid.nlevel, ele.shape[0]))
            I=0
            for ind in ele:
                value = self._grid.trinodes[ind]
                h = np.mean(self._grid.h[value])
                zeta = np.mean(self._var.el[:,value],1) + h
                siglay = np.mean(self._grid.siglay[:,value],1)
                depth[:,:,I] =  zeta[:,None]*siglay[None,:]
                I+=1
            # Average depth over time
            if not argtime==[]:
                depth = np.mean(depth[argtime,:,:], 0)
            else:
                depth = np.mean(depth, 0)
              
            # Compute distance along line
            x = self._grid.xc[ele]
            y = self._grid.yc[ele]
            # Pythagore + cumulative path 
            line = np.zeros(depth.shape)
            dl = np.sqrt(np.square(x[1:]-x[:-1]) + np.square(y[1:]-y[:-1]))
            for i in range(1,dl.shape[0]):
                dl[i] = dl[i] + dl[i-1]
            line[:,1:] = dl[:]
           
            #turn into gridded
            #print 'Compute gridded data'
            #nx, ny = 100, 100
            #xi = np.linspace(x.min(), x.max(), nx)
            #yi = np.linspace(y.min(), y.max(), ny)

            #Plot features
            #setting limits and levels of colormap
            if cmax==[]:
                cmax = varP[:].max()
            if cmin==[]:
                cmin = varP[:].min()
            step = (cmax-cmin) / 20.0
            levels=np.arange(cmin, (cmax+step), step)
            #plt.clf()
            fig = plt.figure(figsize=(18,10))
            plt.rc('font',size='22')
            ax = fig.add_subplot(111) #,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
            #levels = np.linspace(0,3.3,34)
            #cs = ax.contourf(line,depth,varP,levels=levels, cmap=plt.cm.jet)
            cs = ax.contourf(line,depth,varP,levels=levels, vmax=cmax,vmin=cmin,
                              cmap=plt.get_cmap('jet'))
            cbar = fig.colorbar(cs)
            #cbar.set_label(title, rotation=-90,labelpad=30)
            ax.contour(line,depth,varP,cs.levels) #, linewidths=0.5,colors='k')
            #ax.set_title()
            plt.title(title)
            #scale = 1
            #ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
            #ax.xaxis.set_major_formatter(ticks)
            #ax.yaxis.set_major_formatter(ticks)
            plt.xlabel('Distance along line (m)')
            plt.ylabel('Depth (m)')
