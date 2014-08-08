#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
from jdcal import gcal2jd
import numpy as np
import matplotlib.tri as Tri
from regioner import *

class _load_var:
    """
'Variables' subset in FVCOM class contains the numpy arrays:
-----------------------------------------------------------
  Some variables are directly passed on from FVCOM output
  (i.e. el, julianTime, w, u, v, ua, va) and possess in-build
  set of descriptors, ex:
         _long_name = 'Vertically Averaged x-velocity'
        |_units = 'meters s-1'
    ua._|_grid = 'fvcom_grid'
        |_...
  Some others shall be generated as methods are being called, ex:
    hori_velo_norm = horizontal velocity norm
    velo_norm = velocity norm
    verti_shear = vertical shear
    vorticity...            
    """
    def __init__(self, data, grid, tx, History, debug=False):
        self._debug = debug
        #Pointer to History
        self._History = History
        History = self._History

        #Check if time period defined
        self.julianTime = data.variables['time']      
        if tx:
            #Time period           
            region_t = self._t_region(tx, debug=debug)
            #Quick reshape
            region_t = region_t.T[0,:]
            self._region_time = region_t
            # get time and adjust it to matlab datenum
            self.julianTime = data.variables['time'][region_t]
            self.matlabTime = self.julianTime + 678942
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]

            #Check if bounding box has been defined
            if grid._ax==[]:
                if debug:
                    print 'Caching variables...'
                # elev timeseries
                #TR: accelerating caching but have a memory cost
                #temp = data.variables['zeta'][:]
                #self.el = temp[region_t,:]
                #TR Comment: using scipy netcdf we can use .data
                #self.el = data.variables['zeta'][region_t,:]
                self.el = data.variables['zeta'].data[region_t,:]            
                try:
                    #TR: accelerating caching but have a memory cost
                    #temp = data.variables['ww'][:]
                    #self.w = temp[region_t,:,:]
                    #temp = data.variables['u'][:]
                    #self.u = temp[region_t,:,:]
                    #temp = data.variables['v'][:]
                    #self.v = temp[region_t,:,:]
                    #temp = data.variables['ua'][:]
                    #self.ua = temp[region_t,:]
                    #temp = data.variables['va'][:]
                    #self.va = temp[region_t,:]
                    #TR Comment: using scipy netcdf we can use .data
                    #self.w = data.variables['ww'][region_t,:,:]
                    #self.u = data.variables['u'][region_t,:,:]
                    #self.v = data.variables['v'][region_t,:,:]
                    #self.ua = data.variables['ua'][region_t,:]
                    #self.va = data.variables['va'][region_t,:]
                    self.w = data.variables['ww'].data[region_t,:,:]
                    self.u = data.variables['u'].data[region_t,:,:]
                    self.v = data.variables['v'].data[region_t,:,:]
                    self.ua = data.variables['ua'].data[region_t,:]
                    self.va = data.variables['va'].data[region_t,:]
                    # invisible variables
                    self._3D = True
                except KeyError:
                    #temp = data.variables['ua'][:]
                    #self.ua = temp[region_t,:]
                    #temp = data.variables['va'][:]
                    #self.va = temp[region_t,:]
                    #TR Comment: using scipy netcdf we can use .data
                    #self.ua = data.variables['ua'][region_t,:]
                    #self.va = data.variables['va'][region_t,:]
                    self.ua = data.variables['ua'].data[region_t,:]
                    self.va = data.variables['va'].data[region_t,:]
                    self._3D = False
            else:
                if debug:
                    print 'Caching variables...'
                #TR comment: very slow...gonna need optimisation down the line  
                #Bounding box
                region_e = grid._element_index
                region_n = grid._node_index
                #Redefine variables in bounding box & time period
                # elev timeseries
                #TR: accelerating caching but have a memory cost
                #temp = data.variables['zeta'][:]
                #self.el = temp[region_t,region_n]
                #TR Comment: using scipy netcdf we can use .data
                #self.el = data.variables['zeta'][region_t,region_n]
                self.el = data.variables['zeta'].data[region_t,region_n] 
                try:
                    #TR: accelerating caching but have a memory cost
                    #temp = data.variables['ww'][:]
                    #self.w = temp[region_t,:,region_e]
                    #temp = data.variables['u'][:]
                    #self.u = temp[region_t,:,region_e]
                    #temp = data.variables['v'][:]
                    #self.v = temp[region_t,:,region_e]
                    #temp = data.variables['ua'][:]
                    #self.ua = temp[region_t,region_e]
                    #temp = data.variables['va'][:]
                    #self.va = temp[region_t,region_e]
                    #TR Comment: using scipy netcdf we can use .data
                    #self.w = data.variables['ww'][region_t,:,region_e]
                    #self.u = data.variables['u'][region_t,:,region_e]
                    #self.v = data.variables['v'][region_t,:,region_e]
                    #self.ua = data.variables['ua'][region_t,region_e]
                    #self.va = data.variables['va'][region_t,region_e]
                    self.w = data.variables['ww'].data[region_t,:,region_e]
                    self.u = data.variables['u'].data[region_t,:,region_e]
                    self.v = data.variables['v'].data[region_t,:,region_e]
                    self.ua = data.variables['ua'].data[region_t,region_e]
                    self.va = data.variables['va'].data[region_t,region_e]
                    # invisible variables
                    self._3D = True
                except KeyError:
                    #TR Comment: using scipy netcdf we can use .data
                    #self.ua = data.variables['ua'][region_t,region_e]
                    #self.va = data.variables['va'][region_t,region_e]
                    self.ua = data.variables['ua'].data[region_t,region_e]
                    self.va = data.variables['va'].data[region_t,region_e]
                    self._3D = False          
        else:
            #-Append message to History field
            text = 'Full temporal domain'
            self._History.append(text)
            # get time and adjust it to matlab datenum
            self.julianTime = data.variables['time']
            self.matlabTime = self.julianTime[:] + 678942
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]

            #Check if bounding box has been defined
            if grid._ax==[]:
                if debug:
                    print 'Linking variables...'
                # elev timeseries
                self.el = data.variables['zeta']           
                try:
                    self.w = data.variables['ww']
                    self.u = data.variables['u']
                    self.v = data.variables['v']
                    self.ua = data.variables['ua']
                    self.va = data.variables['va']
                    # invisible variables
                    self._3D = True
                except KeyError:
                    self.ua = data.variables['ua']
                    self.va = data.variables['va']
                    self._3D = False
            else:
                if debug:
                    print 'Caching variables...'
                #TR comment: very slow...gonna need optimisation down the line   
                #Bounding box
                region_e = grid._element_index
                region_n = grid._node_index
                #Redefine variables in bounding box & time period
                # elev timeseries
                #TR: accelerating caching but have a memory cost
                #temp = data.variables['zeta'][:]
                #self.el = temp[:,region_n]
                #TR Comment: using scipy netcdf we can use .data
                #self.el = data.variables['zeta'][:,region_n]
                self.el = data.variables['zeta'].data[:,region_n]
                try:
                    #TR: accelerating caching but have a memory cost
                    #temp = data.variables['ww'][:]
                    #self.w = temp[:,:,region_e]
                    #temp = data.variables['u'][:]
                    #self.u = temp[:,:,region_e]
                    #temp = data.variables['v'][:]
                    #self.v = temp[:,:,region_e]
                    #temp = data.variables['ua'][:]
                    #self.ua = temp[:,region_e]
                    #temp = data.variables['va'][:]
                    #self.va = temp[:,region_e]
                    #TR Comment: using scipy netcdf we can use .data
                    #self.w = data.variables['ww'][:,:,region_e]
                    #self.u = data.variables['u'][:,:,region_e]
                    #self.v = data.variables['v'][:,:,region_e]
                    #self.ua = data.variables['ua'][:,region_e]
                    #self.va = data.variables['va'][:,region_e]
                    #TR comment: looping on time indexes is a trick from Mitchell
                    #            to improve loading time
                    self.w = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                    self.u = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                    self.v = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                    self.ua = np.zeros((grid.ntime, grid.nele))
                    self.va = np.zeros((grid.ntime, grid.nele))
                    for i in range(grid.ntime):
                        self.w[i,:,:] = np.transpose(data.variables['ww'].data[i,:,region_e])
                        self.u[i,:,:] = np.transpose(data.variables['u'].data[i,:,region_e])
                        self.v[i,:,:] = np.transpose(data.variables['v'].data[i,:,region_e])
                        self.ua[i,:] = np.transpose(data.variables['ua'].data[i,region_e])
                        self.va[i,:] = np.transpose(data.variables['va'].data[i,region_e])
                    #TR comment: no idea why I have to transpose here but I do !!!
                    # invisible variables
                    self._3D = True
                except KeyError:
                    #TR: accelerating caching but have a memory cost
                    #temp = data.variables['ua'][:]
                    #self.ua = temp[:,region_e]
                    #temp = data.variables['va'][:]
                    #self.va = temp[:,region_e]  
                    #TR Comment: using scipy netcdf we can use .data                 
                    #self.ua = data.variables['ua'][:,region_e]
                    #self.va = data.variables['va'][:,region_e]
                    self.ua = data.variables['ua'].data[:,region_e]
                    self.va = data.variables['va'].data[:,region_e]
                    # invisible variables
                    self._3D = False
        if debug:
            print '...Passed'

    def _t_region(self, tx, quiet=False, debug=False):
        '''Return time indexes included in time period, aka tx'''
        debug = debug or self._debug      
        if debug:
            print 'Computing region_t...'
        if not tx:
            region_t = range(self.julianTime.shape[0])
        else:
            #Conversion to julian day
            dStart = tx[0].split('.')
            dEnd = tx[1].split('.')
            tS = gcal2jd(dStart[0], dStart[1],dStart[2])[1]
            tE = gcal2jd(dEnd[0], dEnd[1],dEnd[2])[1]
            #finding time period
            region_t = np.argwhere((self.julianTime >= tS) &
                                   (self.julianTime <= tE))          
        if debug:
            print '...Passed'

        # Add metadata entry
        if not quiet:
            text = 'Time period =' + str(tx)
            self._History.append(text)
            print '-Now working in time box-'
        return region_t

class _load_grid:
    '''
'Grid' subset in FVCOM class contains grid related quantities:
-------------------------------------------------------------
  Each variable possesses in-build set of descriptors in Data, ex:
                             _long_name = 'zonal longitude'
                            |_units = 'degrees_east'
    Data.varialbes['lonc']._|_dimensions = 'nele'
                            |_...
  Some others shall be generated as methods are being called, ex:
    triangle = triangulation object for plotting purposes
    ...     
    '''
    def __init__(self, data, ax, History, debug=False):
        debug = debug or self._debug     
        if debug:
            print 'Caching grid...'
        #Pointer to History
        self._History = History
        History = self._History
        #Load grid variables on the entire domain:
        #self.lon = data.variables['lon'][:]
        #self.lat = data.variables['lat'][:]
        #self.lonc = data.variables['lonc'][:]
        #self.latc = data.variables['latc'][:]
        #self.x = data.variables['x'][:]
        #self.y = data.variables['y'][:]
        #self.xc = data.variables['xc'][:]
        #self.yc = data.variables['yc'][:]
        #self.a1u = data.variables['a1u'][:]
        #self.a2u = data.variables['a2u'][:]
        #self.aw0 = data.variables['aw0'][:]
        #self.awx = data.variables['awx'][:]
        #self.awy = data.variables['awy'][:]
        #self.trinodes = data.variables['nv'][:].T - 1
        #self.triele = data.variables['nbe'][:].T
        #TR Comment: using scipy netcdf we can use .data
        self.lon = data.variables['lon'].data
        self.lat = data.variables['lat'].data
        self.lonc = data.variables['lonc'].data
        self.latc = data.variables['latc'].data
        self.x = data.variables['x'].data
        self.y = data.variables['y'].data
        self.xc = data.variables['xc'].data
        self.yc = data.variables['yc'].data
        self.a1u = data.variables['a1u'].data
        self.a2u = data.variables['a2u'].data
        self.aw0 = data.variables['aw0'].data
        self.awx = data.variables['awx'].data
        self.awy = data.variables['awy'].data
        self.trinodes = np.transpose(data.variables['nv'].data) - 1
        self.triele = np.transpose(data.variables['nbe'].data)
        if ax==[]:
            #Append message to History field
            text = 'Full spatial domain'
            self._History.append(text)
            #Define the rest of the grid variables
            self.h = data.variables['h'][:]
            try:
                self.siglay = data.variables['siglay'][:]
                self.siglev = data.variables['siglev'][:]
                self.nlevel = self.siglay.shape[0]
            except KeyError:
                pass
            self.nele = data.dimensions['nele']
            self.node = data.dimensions['node']
            #Define bounding box
            self._ax = ax
        else:
            print 'Re-indexing may take some time...'   
            Data = regioner(self, ax, debug=debug)   
            self.lon = Data['lon'][:]
            self.lat = Data['lat'][:]
            self.lonc = Data['lonc'][:]
            self.latc = Data['latc'][:]
            self.x = Data['x'][:]
            self.y = Data['y'][:]
            self.xc = Data['xc'][:]
            self.yc = Data['yc'][:]
            self.a1u = Data['a1u'][:]
            self.a2u = Data['a2u'][:]
            self.aw0 = Data['aw0'][:]
            self.awx = Data['awx'][:]
            self.awy = Data['awy'][:]
            self.trinodes = Data['nv'][:]
            self.triele = Data['nbe'][:]
            self.triangle = Data['triangle']
            #Only load the element within the box
            self._node_index = Data['node_index']
            self._element_index = Data['element_index']
            #TR Comment: using scipy netcdf we can use .data
            #self.h = data.variables['h'][self._node_index]
            self.h = data.variables['h'].data[self._node_index]
            try:
                #self.siglay = data.variables['siglay'][:,self._node_index]
                #self.siglev = data.variables['siglev'][:,self._node_index]
                #TR Comment: using scipy netcdf we can use .data
                self.siglay = data.variables['siglay'].data[:,self._node_index]
                self.siglev = data.variables['siglev'].data[:,self._node_index]
                self.nlevel = self.siglay.shape[0]
            except KeyError:
                pass
            self.nele = Data['element_index'].shape[0]
            self.node = Data['node_index'].shape[0]
            del Data
            #Define bounding box
            self._ax = ax
            # Add metadata entry
            text = 'Bounding box =' + str(ax)
            self._History.append(text)
            print '-Now working in bounding box-'
    
        if debug:
            print '...Passed'

    #TR comment: probably not needed anymore
    #def _ele_region(self, ax, debug=False):
    #    '''Return element indexes included in bounding box, aka ax'''       
    #    if debug:
    #        print 'Computing region_e...'

    #    region_e = np.argwhere((self.lonc >= ax[0]) &
    #                           (self.lonc <= ax[1]) &
    #                           (self.latc >= ax[2]) &
    #                           (self.latc <= ax[3]))          
    #    if debug or self._debug:
    #        print '...Passed'

    #    return region_e

    #TR comment: probably not needed anymore
    #def _node_region(self, ax, debug=False):
    #    '''Return node indexes included in bounding box, aka ax'''
    #    if debug:
    #        print 'Computing region_n...'

    #    region_n = np.argwhere((self.lon >= ax[0]) &
    #                           (self.lon <= ax[1]) &
    #                           (self.lat >= ax[2]) &
    #                           (self.lat <= ax[3]))
    #    if debug or self._debug:
    #        print '...Passed'

    #    return region_n

    #TR comment: probably not needed anymore
    #def _bounding_box(self, ax, quiet=False, debug=False):
    #    """
    #    Define bounding box and reset the box by default.
    #    Input ex:
    #    -------- 
    #      _bounding_box(ax=[min lon, max lon, min lat, max lat])
    #    """
    #    if not ax:
    #        region_e = range(self.nele)
    #        region_n = range(self.node)
    #    else:
    #        region_e = self._ele_region(ax, debug=debug)
    #        region_n = self._node_region(ax, debug=debug)

    #    # Add metadata entry
    #    if not quiet:
    #        text = 'Bounding box =' + str(ax)
    #        self._History.append(text)
    #        print '-Now working in bounding box-'
    #    return region_e, region_n

