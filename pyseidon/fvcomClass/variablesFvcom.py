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
from miscellaneous import time_to_index
from miscellaneous import mattime_to_datetime

class _load_var:
    """
'Variables' subset in FVCOM class contains the numpy arrays:
-----------------------------------------------------------
Some variables are directly passed on from FVCOM output:
                   _el = elevation (m), 2D array (ntime, nnode)
                  |_julianTime = julian date, 1D array (ntime)
                  |_matlabTime = matlab time, 1D array (ntime)
                  |_ua = depth averaged u velocity component (m/s),
                  |      2D array (ntime, nele)
                  |_va = depth averaged v velocity component (m/s),
 FVCOM.Variables._|      2D array (ntime, nele)
                  |_u = u velocity component (m/s),
                  |     3D array (ntime, nlevel, nele)
                  |_v = v velocity component (m/s),
                  |     3D array (ntime, nlevel, nele)
                  |_w = w velocity component (m/s),
                        3D array (ntime, nlevel, nele)

Some others shall be generated as methods are being called, ex:
                  ...
                  |_hori_velo_norm = horizontal velocity norm (m/s),
                  |                  2D array (ntime, nele)
                  |_velo_norm = velocity norm (m/s),
                  |             3D array (ntime, nlevel, nele)
                  |_verti_shear = vertical shear (1/s),
                  |               3D array (ntime, nlevel, nele)
                  |_vorticity...            
    """
    def __init__(self, data, grid, tx, History, debug=False):
        self._debug = debug
        #Pointer to History
        self._History = History
        History = self._History

        #List of keywords
        kwl2D = ['ua', 'va', 'zeta','depth_av_flow_dir', 'hori_velo_norm',
                 'depth_av_vorticity', 'depth_av_power_density',
                 'depth_av_power_assessment']
        kwl3D = ['ww', 'u', 'v', 'gls', 'tke', 'flow_dir', 'velo_norm',
                 'verti_shear', 'vorticity', 'power_density']
        #List of aliaSes
        al2D = ['ua', 'va', 'el','depth_av_flow_dir', 'hori_velo_norm',
               'depth_av_vorticity', 'depth_av_power_density',
               'depth_av_power_assessment']
        al3D = ['w', 'u', 'v', 'gls', 'tke', 'flow_dir', 'velo_norm',
                 'verti_shear', 'vorticity', 'power_density'] 

        #Check if time period defined
        try:
            self.julianTime = data.variables['time']
        except KeyError: #exeception due to Save_as(netcdf)
            self.julianTime = data.variables['julianTime']      
        if not tx==[]:
            #Time period           
            region_t = self._t_region(tx, debug=debug)
            #Quick reshape
            #region_t = region_t.T[0,:]
            self._region_time = region_t
            # define time bounds
            ts = region_t[0]
            te = region_t[-1] + 1
            # get time and adjust it to matlab datenum
            self.julianTime = data.variables['time'][ts:te]
            self.matlabTime = self.julianTime + 678942.0
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]
            if debug: print "ntime: ", grid.ntime
            if debug: print "region_t shape: ", region_t.shape

            #Check if bounding box has been defined
            if grid._ax==[]:
                if debug:
                    print 'Loading variables...'
                #Check if OpenDap variables or not
                if type(data.variables).__name__=='DatasetType':
                    #loading hori data
                    keyCount = 0
                    for key, aliaS in zip(kwl2D, al2D):
                        try:
                            try:
                                setattr(self, aliaS, data.variables[key].data[ts:te,:])
                            except AttributeError: #exeception due nc.Dataset
                                setattr(self, aliaS, data.variables[key][ts:te,:]) 
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Horizontal variables are missing---"
                    self._3D = False 
                    
                    #loading verti data
                    keyCount = 0
                    for key, aliaS in zip(kwl3D, al3D):
                        try:
                            try:
                                setattr(self, aliaS, data.variables[key].data[ts:te,:,:])
                            except AttributeError: #exeception due nc.Dataset
                                setattr(self, aliaS, data.variables[key][ts:te,:,:])
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Vertical variables are missing---"
                    else:
                        self._3D = True

                #Not OpenDap
                else:
                    #loading hori data
                    keyCount = 0
                    for key, aliaS in zip(kwl2D, al2D):
                        try:
                            if key=='zeta':
                                setattr(self, aliaS, np.zeros((grid.ntime, grid.nnode)))
                            else:
                                setattr(self, aliaS, np.zeros((grid.ntime, grid.nele)))
                            I = 0
                            for i in region_t:
                                #TR comment: looping on time indices is a trick from
                                #            Mitchell to improve loading time
                                #TR comment: no idea why I have to transpose here but
                                #            I do !!
                                try:
                                    getattr(self, aliaS)[I,:] =\
                                    np.transpose(data.variables[key].data[i,:])
                                except AttributeError: #exeception due nc.Dataset
                                    getattr(self, aliaS)[I,:] = data.variables[key][i,:]
                                I += 1
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Horizontal variables are missing---"
                    self._3D = False 
                    
                    #loading verti data
                    keyCount = 0
                    for key, aliaS in zip(kwl3D, al3D):
                        try:
                            testKey = data.variables[key]
                            del testKey
                            if aliaS == 'verti_shear':#exception for vertical shear
                                setattr(self, aliaS,
                                np.zeros((grid.ntime,grid.nlevel-1,grid.nele)))
                            else:
                                setattr(self, aliaS,
                                np.zeros((grid.ntime,grid.nlevel, grid.nele)))

                            I = 0
                            for i in region_t:
                                #TR comment: looping on time indices is a trick from
                                #            Mitchell to improve loading time
                                #TR comment: no idea why I have to transpose here but
                                #            I do !!
                                try:
                                    try:
                                        getattr(self, aliaS)[I,:,:] =\
                                        np.transpose(data.variables[key].data[i,:,:])
                                    except AttributeError: #exeception due nc.Dataset
                                        getattr(self, aliaS)[I,:,:] =\
                                        data.variables[key][i,:,:]
                                except ValueError:
                                    try:
                                        getattr(self, aliaS)[I,:,:] =\
                                        data.variables[key].data[i,:,:]
                                    except AttributeError: #exeception due nc.Dataset
                                        getattr(self, aliaS)[I,:,:] =\
                                        data.variables[key][i,:,:]
                                I += 1
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Vertical variables are missing---"
                    else:
                        self._3D = True 

            #Time period defined and region defined           
            else:
                if debug:
                    print 'Loading variables...'
                #TR comment: very slow...gonna need optimisation down the line  
                #Bounding box
                region_e = grid._element_index
                region_n = grid._node_index
                #Redefine variables in bounding box & time period
                #Check if OpenDap variables or not
                if type(data.variables).__name__=='DatasetType':
                    #Special loading for zeta
                    H = 0 #local counter
                    key = kwl2D.pop(2)
                    aliaS = al2D.pop(2)
                    for k, g in groupby(enumerate(region_n), lambda (i,x):i-x):
                        ID = map(itemgetter(1), g)
                        if debug: print 'Index bound: ' +\
                            str(ID[0]) + '-' + str(ID[-1]+1)
                        if H==0:
                            try:
                                setattr(self, aliaS,
                                data.variables[key].data[ts:te,ID[0]:(ID[-1]+1)])
                            except AttributeError: #exeception due nc.Dataset
                                setattr(self, aliaS,
                                data.variables[key][ts:te,ID[0]:(ID[-1]+1)])
                            H=1
                        else:
                            try:
                                setattr(self, aliaS,
                                np.hstack((getattr(self, aliaS),
                                data.variables[key].data[ts:te,ID[0]:(ID[-1]+1)])))
                            except AttributeError: #exeception due nc.Dataset
                                setattr(self, aliaS,
                                np.hstack((getattr(self, aliaS),
                                data.variables[key][ts:te,ID[0]:(ID[-1]+1)])))
                    #loading hori data
                    keyCount = 0
                    for key, aliaS in zip(kwl2D, al2D):
                        try:
                            H = 0 #local counter
                            for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                                ID = map(itemgetter(1), g)
                                if debug: print 'Index bound: ' +\
                                          str(ID[0]) + '-' + str(ID[-1]+1)
                                if H==0:
                                    try:
                                        setattr(self, aliaS,
                                        data.variables[key].data[ts:te,ID[0]:(ID[-1]+1)])
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS,
                                        data.variables[key][ts:te,ID[0]:(ID[-1]+1)])
                                    H=1
                                else:
                                    try:
                                        setattr(self, aliaS,
                                        np.hstack((getattr(self, aliaS),
                                        data.variables[key].data[ts:te,ID[0]:(ID[-1]+1)])))
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS,
                                        np.hstack((getattr(self, aliaS),
                                        data.variables[key][ts:te,ID[0]:(ID[-1]+1)])))
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Horizontal variables are missing---"
                    self._3D = False 
                    
                    #loading verti data
                    keyCount = 0
                    for key, aliaS in zip(kwl3D, al3D):
                        try:
                            H = 0 #local counter
                            for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                                ID = map(itemgetter(1), g)
                                if debug: print 'Index bound: ' +\
                                          str(ID[0]) + '-' + str(ID[-1]+1)
                                if H==0:
                                    try:
                                        setattr(self, aliaS, data.variables[key].data\
                                        [ts:te,:,ID[0]:(ID[-1]+1)])
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS, data.variables[key]\
                                        [ts:te,:,ID[0]:(ID[-1]+1)])
                                    H=1
                                else:
                                    try:
                                        setattr(self, aliaS,
                                        np.dstack((getattr(self, aliaS),
                                        data.variables[key].data[ts:te,:,ID[0]:(ID[-1]+1)])))
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS,
                                        np.dstack((getattr(self, aliaS),
                                        data.variables[key][ts:te,:,ID[0]:(ID[-1]+1)])))
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Vertical variables are missing---"
                    else:
                         self._3D = True

                #Not OpenDap
                else:
                    #loading hori data
                    keyCount = 0
                    for key, aliaS in zip(kwl2D, al2D):
                        try:
                            if key=='zeta':
                                setattr(self, aliaS, np.zeros((grid.ntime, grid.nnode)))
                                I = 0
                                for i in region_t:
                                    #TR comment: looping on time indices is a trick from
                                    #            Mitchell to improve loading time
                                    #TR comment: no idea why I have to transpose here but
                                    #            I do !!
                                    try:
                                        getattr(self, aliaS)[I,:] =\
                                        np.transpose(data.variables[key].data[i,region_n])
                                    except AttributeError: #exeception due nc.Dataset
                                        getattr(self, aliaS)[I,:] =\
                                        data.variables[key][i,region_n]
                                    I += 1
                            else:
                                setattr(self, aliaS, np.zeros((grid.ntime, grid.nele)))
                                I = 0
                                for i in region_t:
                                    #TR comment: looping on time indices is a trick from
                                    #            Mitchell to improve loading time
                                    #TR comment: no idea why I have to transpose here but
                                    #            I do !!
                                    try:
                                        getattr(self, aliaS)[I,:] =\
                                        np.transpose(data.variables[key].data[i,region_e])
                                    except AttributeError: #exeception due nc.Dataset
                                        getattr(self, aliaS)[I,:] =\
                                        data.variables[key][i,region_e]
                                    I += 1
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Horizontal variables are missing---"
                    self._3D = False 
                    
                    #loading verti data
                    keyCount = 0
                    for key in zip(kwl3D, al3D):
                        try:
                            testKey = data.variables[key]
                            del testKey
                            if aliaS == 'verti_shear':#exception for vertical shear
                                setattr(self, aliaS,
                                np.zeros((grid.ntime,grid.nlevel-1,grid.nele)))
                            else:
                                setattr(self, aliaS,
                                np.zeros((grid.ntime,grid.nlevel, grid.nele)))
                            I=0
                            for i in region_t:
                                #TR comment: looping on time indices is a trick from
                                #            Mitchell to improve loading time
                                #TR comment: no idea why I have to transpose here but
                                #            I do !!
                                try:
                                    getattr(self, aliaS)[I,:,:] =\
                                    np.transpose(data.variables[key].data[i,:,region_e])
                                except AttributeError: #exeception due nc.Dataset
                                    getattr(self, aliaS)[I,:,:] =\
                                    data.variables[key][i,:,region_e]
                                I += 1
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Vertical variables are missing---"
                    else:
                        self._3D = True 
  
        #No time period define    
        else:
            # get time and adjust it to matlab datenum
            try:
                self.julianTime = data.variables['time'].data
            except KeyError: #exeception due to Save_as(netcdf)
                self.julianTime = data.variables['julianTime'] 
            self.matlabTime = self.julianTime[:] + 678942.0
            #-Append message to History field
            start = mattime_to_datetime(self.matlabTime[0])
            end = mattime_to_datetime(self.matlabTime[-1])
            text = 'Full temporal domain from ' + str(start) +\
                   ' to ' + str(end)
            self._History.append(text)
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]

            #Check if bounding box has been defined
            if grid._ax==[]:
                if debug:
                    print 'Linking variables...'

                #loading hori data
                keyCount = 0
                for key, aliaS in zip(kwl2D, al2D):
                    try:
                        try:
                            setattr(self, aliaS, data.variables[key].data)
                        except AttributeError: #exeception due nc.Dataset
                            setattr(self, aliaS, data.variables[key])#[:])
                        keyCount +=1
                    except KeyError:
                        if debug: print key, " is missing !"
                        continue
                if keyCount==0:
                    print "---Horizontal variables are missing---"
                self._3D = False 
                    
                #loading verti data
                keyCount = 0
                for key, aliaS in zip(kwl3D, al3D):
                    try:
                        try:
                            setattr(self, aliaS, data.variables[key].data)
                        except AttributeError: #exeception due nc.Dataset
                            setattr(self, aliaS, data.variables[key])#[:])
                        keyCount +=1
                    except KeyError:
                        if debug: print key, " is missing !"
                        continue
                if keyCount==0:
                    print "---Vertical variables are missing---"
                else:
                    self._3D = True

            #No time period defined but region defined
            else:
                if debug:
                    print 'Loading variables...'
                #TR comment: very slow...gonna need optimisation down the line   
                #Bounding box
                region_e = grid._element_index
                region_n = grid._node_index
                #Redefine variables in bounding box
                #Check if OpenDap variables or not
                if type(data.variables).__name__=='DatasetType':
                    #loading hori data
                    keyCount = 0
                    for key, aliaS in zip(kwl2D, al2D):
                        #Special loading for zeta
                        H = 0 #local counter
                        if key == 'zeta':
                            for k, g in groupby(enumerate(region_n), lambda (i,x):i-x):
                                ID = map(itemgetter(1), g)
                                if debug: print 'Index bound: ' +\
                                    str(ID[0]) + '-' + str(ID[-1]+1)
                                if H==0:
                                    try:
                                        setattr(self, aliaS,
                                        data.variables[key].data[:,ID[0]:(ID[-1]+1)])
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS,
                                        data.variables[key][:,ID[0]:(ID[-1]+1)])
                                    H=1
                                else:
                                    try:
                                        setattr(self, aliaS,
                                        np.hstack((getattr(self, aliaS),
                                        data.variables[key].data[:,ID[0]:(ID[-1]+1)])))
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS,
                                        np.hstack((getattr(self, aliaS),
                                        data.variables[key][:,ID[0]:(ID[-1]+1)])))
                        else:
                            try:                        
                                for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                                    ID = map(itemgetter(1), g)
                                    if debug: print 'Index bound: ' +\
                                              str(ID[0]) + '-' + str(ID[-1]+1)
                                    if H==0:
                                        try:
                                            setattr(self, aliaS,
                                            data.variables[key].data[:,ID[0]:(ID[-1]+1)])
                                        except AttributeError: #exeception due nc.Dataset
                                            setattr(self, aliaS,
                                            data.variables[key][:,ID[0]:(ID[-1]+1)])
                                        H=1
                                    else:
                                        try:
                                            setattr(self, aliaS,
                                            np.hstack((getattr(self, aliaS),
                                            data.variables[key].data[:,ID[0]:(ID[-1]+1)])))
                                        except AttributeError: #exeception due nc.Dataset
                                            setattr(self, aliaS,
                                            np.hstack((getattr(self, aliaS),
                                            data.variables[key][:,ID[0]:(ID[-1]+1)])))
                                    keyCount +=1
                            except KeyError:
                                if debug: print key, " is missing !"
                                continue
                    if keyCount==0:
                        print "---Horizontal variables are missing---"
                    self._3D = False 
                    
                    #loading verti data
                    keyCount = 0
                    for key, aliaS in zip(kwl3D, al3D):
                        try:
                            H = 0 #local counter
                            for k, g in groupby(enumerate(region_e), lambda (i,x):i-x):
                                ID = map(itemgetter(1), g)
                                if debug: print 'Index bound: ' +\
                                          str(ID[0]) + '-' + str(ID[-1]+1)
                                if H==0:
                                    try:
                                        setattr(self, aliaS,
                                        data.variables[key].data[:,:,ID[0]:(ID[-1]+1)])
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS,
                                        data.variables[key][:,:,ID[0]:(ID[-1]+1)])
                                    H=1
                                else:
                                    try:
                                        setattr(self, aliaS,
                                        np.dstack((getattr(self, aliaS),
                                        data.variables[key].data[:,:,ID[0]:(ID[-1]+1)])))
                                    except AttributeError: #exeception due nc.Dataset
                                        setattr(self, aliaS,
                                        np.dstack((getattr(self, aliaS),
                                        data.variables[key][:,:,ID[0]:(ID[-1]+1)])))
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Vertical variables are missing---"
                    else:
                        self._3D = True

                #Not OpenDap
                else:
                    #loading hori data
                    keyCount = 0
                    for key, aliaS in zip(kwl2D, al2D):
                        try:
                            if key=='zeta':
                                setattr(self, aliaS, np.zeros((grid.ntime, grid.nnode)))
                                for i in range(grid.ntime):
                                    #TR comment: looping on time indices is a trick from
                                    #            Mitchell to improve loading time
                                    #TR comment: no idea why I have to transpose here but
                                    #            I do !!
                                    try:
                                        getattr(self, aliaS)[i,:] =\
                                        np.transpose(data.variables[key].data[i,region_n])
                                    except AttributeError: #exeception due nc.Dataset
                                        getattr(self, aliaS)[i,:] =\
                                        data.variables[key][i,region_n]
                            else:
                                setattr(self, aliaS, np.zeros((grid.ntime, grid.nele)))
                                for i in range(grid.ntime):
                                    #TR comment: looping on time indices is a trick from
                                    #            Mitchell to improve loading time
                                    #TR comment: no idea why I have to transpose here but
                                    #            I do !!
                                    try:
                                        getattr(self, aliaS)[i,:] =\
                                        np.transpose(data.variables[key].data[i,region_e])
                                    except AttributeError: #exeception due nc.Dataset
                                        getattr(self, aliaS)[i,:] =\
                                        data.variables[key][i,region_e]
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Horizontal variables are missing---"
                    self._3D = False 
                    
                    #loading verti data
                    keyCount = 0
                    for key, aliaS in zip(kwl3D, al3D):
                        try:
                            testKey = data.variables[key]
                            del testKey
                            if aliaS == 'verti_shear':#exception for vertical shear
                                setattr(self, aliaS,
                                np.zeros((grid.ntime,grid.nlevel-1,grid.nele)))
                            else:
                                setattr(self, aliaS,
                                np.zeros((grid.ntime,grid.nlevel, grid.nele)))
                            for i in range(grid.ntime):
                                #TR comment: looping on time indices is a trick from
                                #            Mitchell to improve loading time
                                #TR comment: no idea why I have to transpose here but
                                #            I do !!
                                try:
                                    getattr(self, aliaS)[i,:,:] =\
                                    np.transpose(data.variables[key].data[i,:,region_e])
                                except AttributeError: #exeception due nc.Dataset
                                    getattr(self, aliaS)[i,:,:] =\
                                    data.variables[key][i,:,region_e]
                            keyCount +=1
                        except KeyError:
                            if debug: print key, " is missing !"
                            continue
                    if keyCount==0:
                        print "---Vertical variables are missing---"
                    else:
                        self._3D = True 
        if debug:
            print '...Passed'

    def _t_region(self, tx, debug=False):
        '''Return time indices included in time period, aka tx'''
        debug = debug or self._debug      
        if debug:
            print 'Computing region_t...'
        region_t = time_to_index(tx[0], tx[1],
                   (self.julianTime[:] + 678942.0), debug=debug)       
        if debug:
            print '...Passed'
        # Add metadata entry
        text = 'Time period =' + str(tx)
        self._History.append(text)
        print '-Now working in time box-'

        return region_t

class _load_grid:
    '''
'Grid' subset in FVCOM class contains grid related quantities:
-------------------------------------------------------------
Some grid data are directly passed on from FVCOM output:
              _lon = longitudes at nodes (deg.), 2D array (ntime, nnode)
             |_lonc = longitudes at elements (deg.), 2D array (ntime, nele)
             |_lat = latitudes at nodes (deg.), 2D array (ntime, nnode)
             |_latc = latitudes at elements (deg.), 2D array (ntime, nele)   
 FVCOM.Grid._|_x = x coordinates at nodes (m), 2D array (ntime, nnode)
             |_xc = x coordinates at elements (m), 2D array (ntime, nele)
             |_y = y coordinates at nodes (m), 2D array (ntime, nnode)
             |_yc = y coordinates at nodes (m), 2D array (ntime, nele)
             |_h = bathymetry (m), 2D array (ntime, nnode)
             |_nele = element dimension, integer
             |_nnode = node dimension, integer
             |_nlevel = vertical level dimension, integer
             |_ntime = time dimension, integer
             |_trinodes = surrounding node indices, 2D array (3, nele)
             |_trinodes = surrounding element indices, 2D array (3, nele)
             |_siglay = sigma layers, 2D array (nlevel, nnode)
             |_siglay = sigma levels, 2D array (nlevel+1, nnode)
             |_and a all bunch of grid parameters...
             | i.e. a1u, a2u, aw0, awx, awy


Some others shall be generated as methods are being called, ex:
             ...
             |_triangle = triangulation object for plotting purposes    
    '''
    def __init__(self, data, ax, History, debug=False):
        self._debug = debug   
        if debug:
            print 'Loading grid...'
        #Pointer to History
        self._History = History
        History = self._History
        #list of grid variable
        gridvar = ['lon','lat','lonc','latc','x','y','xc','yc',
                   'a1u','a2u','aw0','awx','awy']
        for key in gridvar:
            try:
                setattr(self, key, data.variables[key].data)
            except AttributeError: #exception for nc.dataset type data
                setattr(self, key, data.variables[key][:])

        #special treatment for triele & trinodes due to Save_as(netcdf)
        datavar = data.variables.keys()
        if "trinodes" in datavar:
            try:
                setattr(self, 'trinodes', data.variables['trinodes'].data)
            except AttributeError: #exception for nc.dataset type data
                setattr(self, 'trinodes', data.variables['trinodes'][:])
        else:
            self.trinodes = np.transpose(data.variables['nv'].data) - 1
        if "triele" in datavar:
            try:
                setattr(self, 'triele', data.variables['triele'].data)
            except AttributeError: #exception for nc.dataset type data
                setattr(self, 'triele', data.variables['triele'][:])
        else:
            self.triele = np.transpose(data.variables['nbe'].data) - 1
        #special treatment for depth2D & depth due to Save_as(netcdf)
        if "depth2D" in datavar:
            setattr(self, "depth2D", data.variables["depth2D"][:])
        if "depth" in datavar:
            setattr(self, "depth", data.variables["depth"][:])

        if ax==[]:
            #Define bounding box
            self._ax = []
            #Append message to History field
            text = 'Full spatial domain'
            self._History.append(text)
            #Define the rest of the grid variables
            self.h = data.variables['h'][:]
            self.siglay = data.variables['siglay'][:]
            self.siglev = data.variables['siglev'][:]
            self.nlevel = self.siglay.shape[0]
            self.nele = self.lonc.shape[0]
            self.nnode = self.lon.shape[0]
            #try:
            #    self.nele = data.dimensions['nele']
            #    self.nnode = data.dimensions['node']
            #except AttributeError:
            #    #TR: bug due to difference in Pydap's data sturcture
            #    self.nele = self.lonc.shape[0]
            #    self.nnode = self.lon.shape[0]
        else:
            #Checking for pre-defined regions
            if ax=='GP':
                ax=[-66.36, -66.31, 44.24, 44.3]
            elif ax=='PP':
                ax=[-66.23, -66.19, 44.37, 44.41]
            elif ax=='DG':
                ax=[-65.84, -65.73, 44.64, 44.72]
            elif ax=='MP':
                ax=[-65.5, -63.3, 45.0, 46.0]
           
            print 'Re-indexing may take some time...'   
            Data = regioner(self, ax, debug=debug)
            #list of grid variable
            gridvar = ['lon','lat','lonc','latc','x','y','xc','yc',
                       'a1u','a2u','aw0','awx','awy','nv','nbe']
            for key in gridvar:
                setattr(self, key, Data[key][:])
            #Special treatment here
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
                        try:
                            self.h = data.variables['h'].data[ID[0]:(ID[-1]+1)]
                            self.siglay = data.variables['siglay'].data[:,ID[0]:(ID[-1]+1)]
                            self.siglev = data.variables['siglev'].data[:,ID[0]:(ID[-1]+1)]
                        except AttributeError: #exception for nc.dataset type data
                            self.h = data.variables['h'][ID[0]:(ID[-1]+1)]
                            self.siglay = data.variables['siglay'][:,ID[0]:(ID[-1]+1)]
                            self.siglev = data.variables['siglev'][:,ID[0]:(ID[-1]+1)]
                    else:
                        try:
                            self.h = np.hstack((self.h,
                            data.variables['h'].data[ID[0]:(ID[-1]+1)]))
                            self.siglay = np.hstack((self.siglay,
                            data.variables['siglay'].data[:,ID[0]:(ID[-1]+1)]))
                            self.siglev = np.hstack((self.siglev,
                            data.variables['siglev'].data[:,ID[0]:(ID[-1]+1)]))
                        except AttributeError: #exception for nc.dataset type data
                            self.h = np.hstack((self.h,
                            data.variables['h'][ID[0]:(ID[-1]+1)]))
                            self.siglay = np.hstack((self.siglay,
                            data.variables['siglay'][:,ID[0]:(ID[-1]+1)]))
                            self.siglev = np.hstack((self.siglev,
                            data.variables['siglev'][:,ID[0]:(ID[-1]+1)]))
                    H=1
            else:
                try:  
                    self.h = data.variables['h'].data[self._node_index]
                    self.siglay = data.variables['siglay'].data[:,self._node_index]
                    self.siglev = data.variables['siglev'].data[:,self._node_index]
                except AttributeError: #exception for nc.dataset type data
                    self.h = data.variables['h'][self._node_index]
                    self.siglay = data.variables['siglay'][:,self._node_index]
                    self.siglev = data.variables['siglev'][:,self._node_index]
            #Dimensions
            self.nlevel = self.siglay.shape[0]
            self.nele = Data['element_index'].shape[0]
            self.nnode = Data['node_index'].shape[0]

            del Data
            #Define bounding box
            self._ax = ax
            # Add metadata entry
            text = 'Bounding box =' + str(ax)
            self._History.append(text)
            print '-Now working in bounding box-'
    
        if debug:
            print '...Passed'
