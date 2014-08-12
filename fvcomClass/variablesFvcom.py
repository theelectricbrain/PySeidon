#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
#from jdcal import gcal2jd
import numpy as np
import matplotlib.tri as Tri
from itertools import groupby
from operator import itemgetter
#Local import
from regioner import *
from miscallaneous import time_to_index

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
        if not tx==[]:
            #Time period           
            region_t = self._t_region(tx, debug=debug)
            #Quick reshape
            region_t = region_t.T[0,:]
            self._region_time = region_t
            # define time bounds
            ts = region_t[0]
            te = region_t[-1]
            # get time and adjust it to matlab datenum
            self.julianTime = data.variables['time'][ts:te]
            self.matlabTime = self.julianTime + 678942
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]

            #Check if bounding box has been defined
            if grid._ax==[]:
                if debug:
                    print 'Loading variables...'
                # elev timeseries
                self.el = data.variables['zeta'].data[ts:te,:]            
                try:
                    #different loading technique if using OpenDap server
                    if type(data.variables).__name__=='DatasetType':
                        self.w = data.variables['ww'].data[ts:te,:,:]
                        self.u = data.variables['u'].data[ts:te,:,:]
                        self.v = data.variables['v'].data[ts:te,:,:]
                        self.ua = data.variables['ua'].data[ts:te,:]
                        self.va = data.variables['va'].data[ts:te,:]
                    else:
                        #TR comment: looping on time indexes is a trick from Mitchell
                        #            to improve loading time
                        #TR comment: no idea why I have to transpose here but I do !!
                        self.w = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.u = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.v = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.ua = np.zeros((grid.ntime, grid.nele))
                        self.va = np.zeros((grid.ntime, grid.nele))
                        for i in region_t:
                            self.w[i,:,:] = np.transpose(
                                            data.variables['ww'].data[i,:,:])
                            self.u[i,:,:] = np.transpose(
                                            data.variables['u'].data[i,:,:])
                            self.v[i,:,:] = np.transpose(
                                            data.variables['v'].data[i,:,:])
                            self.ua[i,:] = np.transpose(
                                           data.variables['ua'].data[i,:])
                            self.va[i,:] = np.transpose(
                                           data.variables['va'].data[i,:])
                    # invisible variables
                    self._3D = True
                except KeyError:
                    #different loading technique if using OpenDap server
                    if type(data.variables).__name__=='DatasetType':
                        self.ua = data.variables['ua'].data[ts:te,:]
                        self.va = data.variables['va'].data[ts:te,:]
                    else:
                        #TR comment: looping on time indexes is a trick from Mitchell
                        #            to improve loading time
                        #TR comment: no idea why I have to transpose here but I do !!
                        self.ua = np.zeros((grid.ntime, grid.nele))
                        self.va = np.zeros((grid.ntime, grid.nele))
                        for i in region_t:
                            self.ua[i,:] = np.transpose(
                                           data.variables['ua'].data[i,:])
                            self.va[i,:] = np.transpose(
                                           data.variables['va'].data[i,:])

                    # invisible variables
                    self._3D = False
            else:
                if debug:
                    print 'Loading variables...'
                #TR comment: very slow...gonna need optimisation down the line  
                #Bounding box
                region_e = grid._element_index
                region_n = grid._node_index
                #Redefine variables in bounding box & time period
                try:
                    #different loading technique if using OpenDap server
                    if type(data.variables).__name__=='DatasetType':
                        #Split into consecutive integers to optimise loading
                        #TR comment: data.variables['ww'].data[:,:,region_n] doesn't
                        #            work with non consecutive indices
                        H=0
                        for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.w = data.variables['ww'].data[ts:te,:,ID[0]:(ID[-1]+1)]
                                self.u = data.variables['u'].data[ts:te,:,ID[0]:(ID[-1]+1)]
                                self.v = data.variables['v'].data[ts:te,:,ID[0]:(ID[-1]+1)]
                                self.ua = data.variables['ua'].data[ts:te,ID[0]:(ID[-1]+1)]
                                self.va = data.variables['va'].data[ts:te,ID[0]:(ID[-1]+1)]
                            else:
                                self.w = np.dstack((self.w,
                                data.variables['ww'].data[ts:te,:,ID[0]:(ID[-1]+1)]))
                                self.u = np.dstack((self.u,
                                data.variables['u'].data[ts:te,:,ID[0]:(ID[-1]+1)]))
                                self.v = np.dstack((self.v,
                                data.variables['v'].data[ts:te,:,ID[0]:(ID[-1]+1)]))
                                self.ua = np.hstack((self.ua,
                                data.variables['ua'].data[ts:te,ID[0]:(ID[-1]+1)]))
                                self.va = np.hstack((self.va,
                                data.variables['va'].data[ts:te,ID[0]:(ID[-1]+1)]))
                            H=1

                        # elev timeseries
                        H=0
                        for k, g in groupby(enumerate(region_n), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.el = data.variables['zeta']\
                                         .data[ts:te,:,ID[0]:(ID[-1]+1)]
                            else:
                                self.el = np.hstack((self.el,
                                data.variables['zeta'].data[ts:te,ID[0]:(ID[-1]+1)]))
                            H=1
                    else:
                        # elev timeseries
                        self.el = data.variables['zeta'].data[ts:te,region_n] 
                        #TR comment: looping on time indexes is a trick from Mitchell
                        #            to improve loading time
                        #TR comment: no idea why I have to transpose here but I do !!
                        self.w = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.u = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.v = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.ua = np.zeros((grid.ntime, grid.nele))
                        self.va = np.zeros((grid.ntime, grid.nele))
                        for i in region_t:
                            self.w[i,:,:] = np.transpose(
                                            data.variables['ww'].data[i,:,region_e])
                            self.u[i,:,:] = np.transpose(
                                            data.variables['u'].data[i,:,region_e])
                            self.v[i,:,:] = np.transpose(
                                            data.variables['v'].data[i,:,region_e])
                            self.ua[i,:] = np.transpose(
                                           data.variables['ua'].data[i,region_e])
                            self.va[i,:] = np.transpose(
                                           data.variables['va'].data[i,region_e])
                    # invisible variables
                    self._3D = True

                except KeyError:
                    #different loading technique if using OpenDap server
                    if type(data.variables).__name__=='DatasetType':
                        #Split into consecutive integers to optimise loading
                        #TR comment: data.variables['ww'].data[:,:,region_n] doesn't
                        #            work with non consecutive indices
                        H=0
                        for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.ua = data.variables['ua'].data[ts:te,ID[0]:(ID[-1]+1)]
                                self.va = data.variables['va'].data[ts:te,ID[0]:(ID[-1]+1)]
                            else:
                                self.ua = np.hstack((self.ua,
                                data.variables['ua'].data[ts:te,ID[0]:(ID[-1]+1)]))
                                self.va = np.hstack((self.va,
                                data.variables['va'].data[ts:te,ID[0]:(ID[-1]+1)]))
                            H=1

                        # elev timeseries
                        H=0
                        for k, g in groupby(enumerate(region_n), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.el = data.variables['zeta']\
                                         .data[:,:,ID[0]:(ID[-1]+1)]
                            else:
                                self.el = np.hstack((self.el,
                                data.variables['zeta'].data[:,ID[0]:(ID[-1]+1)]))
                            H=1
                    else:
                        # elev timeseries
                        self.el = data.variables['zeta'].data[ts:te,region_n] 
                        #TR comment: looping on time indexes is a trick from Mitchell
                        #            to improve loading time
                        #TR comment: no idea why I have to transpose here but I do !!
                        self.ua = np.zeros((grid.ntime, grid.nele))
                        self.va = np.zeros((grid.ntime, grid.nele))
                        for i in region_t:
                            self.ua[i,:] = np.transpose(
                            data.variables['ua'].data[i,region_e])
                            self.va[i,:] = np.transpose(
                            data.variables['va'].data[i,region_e])
                    # invisible variables
                    self._3D = False          
        else:
            #-Append message to History field
            text = 'Full temporal domain'
            self._History.append(text)
            # get time and adjust it to matlab datenum
            self.julianTime = data.variables['time'].data
            self.matlabTime = self.julianTime[:] + 678942
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]

            #Check if bounding box has been defined
            if grid._ax==[]:
                if debug:
                    print 'Linking variables...'
                # elev timeseries
                self.el = data.variables['zeta'].data           
                try:
                    self.w = data.variables['ww'].data
                    self.u = data.variables['u'].data
                    self.v = data.variables['v'].data
                    self.ua = data.variables['ua'].data
                    self.va = data.variables['va'].data
                    # invisible variables
                    self._3D = True
                except KeyError:
                    self.ua = data.variables['ua'].data
                    self.va = data.variables['va'].data
                    self._3D = False
            else:
                if debug:
                    print 'Loading variables...'
                #TR comment: very slow...gonna need optimisation down the line   
                #Bounding box
                region_e = grid._element_index
                region_n = grid._node_index
                #Redefine variables in bounding box
                try:
                    #different loading technique if using OpenDap server
                    if type(data.variables).__name__=='DatasetType':
                        #Split into consecutive integers to optimise loading
                        #TR comment: data.variables['ww'].data[:,:,region_n] doesn't
                        #            work with non consecutive indices
                        H=0
                        for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.w = data.variables['ww'].data[:,:,ID[0]:(ID[-1]+1)]
                                self.u = data.variables['u'].data[:,:,ID[0]:(ID[-1]+1)]
                                self.v = data.variables['v'].data[:,:,ID[0]:(ID[-1]+1)]
                                self.ua = data.variables['ua'].data[:,ID[0]:(ID[-1]+1)]
                                self.va = data.variables['va'].data[:,ID[0]:(ID[-1]+1)]
                            else:
                                self.w = np.dstack((self.w,
                                         data.variables['ww'].data[:,:,ID[0]:(ID[-1]+1)]))
                                self.u = np.dstack((self.u,
                                         data.variables['u'].data[:,:,ID[0]:(ID[-1]+1)]))
                                self.v = np.dstack((self.v,
                                         data.variables['v'].data[:,:,ID[0]:(ID[-1]+1)]))
                                self.ua = np.hstack((self.ua,
                                          data.variables['ua'].data[:,ID[0]:(ID[-1]+1)]))
                                self.va = np.hstack((self.va,
                                          data.variables['va'].data[:,ID[0]:(ID[-1]+1)]))
                            H=1

                        # elev timeseries
                        H=0
                        for k, g in groupby(enumerate(region_n), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.el = data.variables['zeta']\
                                         .data[:,:,ID[0]:(ID[-1]+1)]
                            else:
                                self.el = np.hstack((self.el,
                                data.variables['zeta'].data[:,ID[0]:(ID[-1]+1)]))
                            H=1
                    else:
                        # elev timeseries
                        self.el = data.variables['zeta'].data[:,region_n]
                        #TR comment: looping on time indexes is a trick from Mitchell
                        #            to improve loading time
                        self.w = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.u = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.v = np.zeros((grid.ntime, grid.nlevel, grid.nele))
                        self.ua = np.zeros((grid.ntime, grid.nele))
                        self.va = np.zeros((grid.ntime, grid.nele))
                        for i in range(grid.ntime):
                            #TR comment: no idea why I have to transpose here but I do !!!
                            self.w[i,:,:] = np.transpose(
                                            data.variables['ww'].data[i,:,region_e])
                            self.u[i,:,:] = np.transpose(
                                            data.variables['u'].data[i,:,region_e])
                            self.v[i,:,:] = np.transpose(
                                            data.variables['v'].data[i,:,region_e])
                            self.ua[i,:] = np.transpose(
                                           data.variables['ua'].data[i,region_e])
                            self.va[i,:] = np.transpose(
                                           data.variables['va'].data[i,region_e])
                    # invisible variables
                    self._3D = True
                except KeyError:
                    #different loading technique if using OpenDap server
                    if type(data.variables).__name__=='DatasetType':
                        H=0
                        for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.ua = data.variables['ua'].data[:,ID[0]:(ID[-1]+1)]
                                self.va = data.variables['va'].data[:,ID[0]:(ID[-1]+1)]
                            else:
                                self.ua = np.hstack((self.ua,
                                          data.variables['ua'].data[:,ID[0]:(ID[-1]+1)]))
                                self.va = np.hstack((self.va,
                                          data.variables['va'].data[:,ID[0]:(ID[-1]+1)]))
                            H=1

                        # elev timeseries
                        H=0
                        for k, g in groupby(enumerate(region_n), lambda (i,x):i-x):
                            ID = map(itemgetter(1), g)
                            if debug: print 'Index bound: ' +\
                                      str(ID[0]) + '-' + str(ID[-1]+1)
                            if H==0:
                                self.el = data.variables['zeta']\
                                         .data[:,:,ID[0]:(ID[-1]+1)]
                            else:
                                self.el = np.hstack((self.el,
                                data.variables['zeta'].data[:,ID[0]:(ID[-1]+1)]))
                            H=1
                    else:
                        # elev timeseries
                        self.el = data.variables['zeta'].data[:,region_n]
                        #TR comment: looping on time indexes is a trick from Mitchell
                        #            to improve loading time
                        #TR comment: no idea why I have to transpose here but I do !!!
                        self.ua = np.zeros((grid.ntime, grid.nele))
                        self.va = np.zeros((grid.ntime, grid.nele))
                        for i in range(grid.ntime):
                            self.ua[i,:] = np.transpose(
                                           data.variables['ua'].data[i,region_e])
                            self.va[i,:] = np.transpose(
                                           data.variables['va'].data[i,region_e])
                    # invisible variables
                    self._3D = False
        if debug:
            print '...Passed'

    def _t_region(self, tx, quiet=False, debug=False):
        '''Return time indexes included in time period, aka tx'''
        debug = debug or self._debug      
        if debug:
            print 'Computing region_t...'
        #if not tx:
        #    region_t = range(self.julianTime.shape[0])
        #else:
        #    #Conversion to julian day
        #    dStart = tx[0].split('.')
        #    dEnd = tx[1].split('.')
        #    tS = gcal2jd(dStart[0], dStart[1],dStart[2])[1]
        #    tE = gcal2jd(dEnd[0], dEnd[1],dEnd[2])[1]
        #    #finding time period
        #    region_t = np.argwhere((self.julianTime >= tS) &
        #                           (self.julianTime <= tE)) 
        region_t = time_to_index(tx[0], t[1], self.julian + 678942, debug=debug)
         
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
            print 'Loading grid...'
        #Pointer to History
        self._History = History
        History = self._History
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
            self.siglay = data.variables['siglay'][:]
            self.siglev = data.variables['siglev'][:]
            self.nlevel = self.siglay.shape[0]
            try:
                self.nele = data.dimensions['nele']
                self.node = data.dimensions['node']
            except AttributeError:
                #TR: bug due to difference in Pydap's data sturcture
                self.nele = self.lonc.shape[0]
                self.node = data.lon.shape[0]
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

            #different loading technique if using OpenDap server
            if type(data.variables).__name__=='DatasetType':
                #Split into consecutive integers to optimise loading
                #TR comment: data.variables['ww'].data[:,:,region_n] doesn't
                #            work with non consecutive indices
                H=0
                for k, g in groupby(enumerate(self._node_index), lambda (i,x):i-x):
                    ID = map(itemgetter(1), g)
                    if debug: print 'Index bound: ' +\
                        str(ID[0]) + '-' + str(ID[-1]+1)
                    if H==0:
                        self.h = data.variables['h'].data[ID[0]:(ID[-1]+1)]
                        self.siglay = data.variables['siglay'].data[:,ID[0]:(ID[-1]+1)]
                        self.siglev = data.variables['siglev'].data[:,ID[0]:(ID[-1]+1)]
                    else:
                        self.h = np.hstack((self.h,
                        data.variables['h'].data[ID[0]:(ID[-1]+1)]))
                        self.siglay = np.hstack((self.siglay,
                        data.variables['siglay'].data[:,ID[0]:(ID[-1]+1)]))
                        self.siglev = np.hstack((self.siglev,
                        data.variables['siglev'].data[:,ID[0]:(ID[-1]+1)]))
                    H=1
            else:
                self.h = data.variables['h'].data[self._node_index]
                self.siglay = data.variables['siglay'].data[:,self._node_index]
                self.siglev = data.variables['siglev'].data[:,self._node_index]
            #Dimensions
            self.nlevel = self.siglay.shape[0]
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
