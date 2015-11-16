#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
from itertools import groupby
from operator import itemgetter
import datetime
import gc
# Parallel computing
#import multiprocessing as mp
#Local import
from pyseidon_dvt.utilities.regioner import *
from pyseidon_dvt.utilities.miscellaneous import time_to_index
from pyseidon_dvt.utilities.miscellaneous import mattime_to_datetime

class _load_grid:
    """
    **'Grid' subset in FVCOM class contains grid related quantities**

    Some grid data are directly passed on from FVCOM output: ::
                  _lon = longitudes at nodes (deg.), 2D array (ntime, nnode)
                 |_lonc = longitudes at elements (deg.), 2D array (ntime, nele)
                 |_lat = latitudes at nodes (deg.), 2D array (ntime, nnode)
                 |_latc = latitudes at elements (deg.), 2D array (ntime, nele)
                 |_x = x coordinates at nodes (m), 2D array (ntime, nnode)
                 |_xc = x coordinates at elements (m), 2D array (ntime, nele)
                 |_y = y coordinates at nodes (m), 2D array (ntime, nnode)
                 |_yc = y coordinates at nodes (m), 2D array (ntime, nele)
     FVCOM.Grid._|_h = bathymetry (m), 2D array (ntime, nnode)
                 |_nele = element dimension, integer
                 |_nnode = node dimension, integer
                 |_nlevel = vertical level dimension, integer
                 |_ntime = time dimension, integer
                 |_trinodes = surrounding node indices, 2D array (3, nele)
                 |_triele = surrounding element indices, 2D array (3, nele)
                 |_siglay = sigma layers, 2D array (nlevel, nnode)
                 |_siglay = sigma levels, 2D array (nlevel+1, nnode)
                 |_and a all bunch of grid parameters...
                 | i.e. a1u, a2u, aw0, awx, awy

    Some others shall be generated as methods are being called, ex: ::
                 ...
                 |_triangle = triangulation object for plotting purposes

    """
    def __init__(self, data, ax, History, debug=False):
        self._debug = debug
        if debug:
            print 'Loading grid...'
        #Pointer to History
        setattr(self, '_History', History)

        #list of required grid variable
        gridvar = ['lon','lat','lonc','latc','x','y','xc','yc',
                   'a1u','a2u','aw0','awx','awy']

        # Figure out which quantity to treat
        self._gridvar = []
        for key in gridvar:
            if key in data.variables.keys():
                self._gridvar.append(key)
            else:
                if key in ["a1u", "a2u", "aw0", "awx", "awy"]:
                    print "--- "+key+" is missing. Some interpolation functions will not work ---"
                if key in ["lonc", "latc", "xc", "yc"]:
                    print "--- "+key+" is missing. Some element based functions will not work ---"

        for key in self._gridvar:
            try:
                setattr(self, key, data.variables[key].data)
            except AttributeError: #exception for nc.dataset type data
                setattr(self, key, data.variables[key])#[:])

        #special treatment for triele & trinodes due to Save_as(netcdf)
        datavar = data.variables.keys()
        if "trinodes" in datavar:
            try:
                setattr(self, 'trinodes', data.variables['trinodes'].data)
            except AttributeError: #exception for nc.dataset type data
                setattr(self, 'trinodes', data.variables['trinodes'])#[:])
        elif "nv" in datavar:
            try:
                self.trinodes = np.transpose(data.variables['nv'].data) - 1
            except AttributeError: #exception for nc.dataset type data
                self.trinodes = np.transpose(data.variables['nv'][:]) - 1
        else:
            print "--- surrounding node indices (nv) missing. Some functions will not work ---"
        if "triele" in datavar:
            try:
                setattr(self, 'triele', data.variables['triele'].data)
            except AttributeError: #exception for nc.dataset type data
                setattr(self, 'triele', data.variables['triele'])#[:])
        elif "nbe" in datavar:
            try:
                self.triele = np.transpose(data.variables['nbe'].data) - 1
            except AttributeError: #exception for nc.dataset type data
                self.triele = np.transpose(data.variables['nbe'][:]) - 1
        else:
            print "--- surrounding element indices (nbe) missing. Some functions will not work ---"

        #special treatment for depth2D & depth due to Save_as(netcdf)
        if "depth2D" in datavar:
            setattr(self, "depth2D", data.variables["depth2D"])#[:])
        if "depth" in datavar:
            setattr(self, "depth", data.variables["depth"])#[:])

        if ax==[]:
            #Define bounding box
            self._ax = []
            #Append message to History field
            text = 'Full spatial domain'
            self._History.append(text)
            #Define the rest of the grid variables
            try: self.h = data.variables['h'][:]
            except KeyError: pass
            try: self.siglay = data.variables['siglay'][:]
            except KeyError: pass
            try: self.siglev = data.variables['siglev'][:]
            except KeyError: pass
            try: self.nlevel = self.siglay.shape[0]
            except AttributeError: pass
            try: self.nele = self.lonc.shape[0]
            except AttributeError: pass
            try: self.nnode = self.lon.shape[0]
            except: pass
        else:
            #Checking for pre-defined regions
            if ax=='GP': ax=[-66.36, -66.31, 44.24, 44.3]
            elif ax=='PP': ax=[-66.23, -66.19, 44.37, 44.41]
            elif ax=='DG': ax=[-65.84, -65.73, 44.64, 44.72]
            elif ax=='MP': ax=[-65.5, -63.3, 45.0, 46.0]

            print 'Re-indexing may take some time...'
            Data = regioner(self, ax, debug=debug)
            #list of grid variable
            gridvar = ['lon','lat','lonc','latc','x','y','xc','yc',
                       'a1u','a2u','aw0','awx','awy','nv','nbe']

            # Figure out which quantity to treat
            self._gridvar = []
            for key in gridvar:
                if key in data.variables.keys():
                    self._gridvar.append(key)
                else:
                    if debug: print "Grid related field "+key+" is missing !"

            for key in self._gridvar:
                setattr(self, key, Data[key][:])
            # Special treatment here
            self.trinodes = Data['nv'][:]
            self.triele = Data['nbe'][:]
            self.triangle = Data['triangle']
            # Only load the element within the box
            self._node_index = Data['node_index']
            self._element_index = Data['element_index']

            # different loading technique if using OpenDap server
            if type(data.variables).__name__ == 'DatasetType':
                # Split into consecutive integers to optimise loading
                # TR comment: data.variables['ww'].data[:,:,region_n] doesn't
                #            work with non consecutive indices
                H=0
                for k, g in groupby(enumerate(self._node_index), lambda (i,x):i-x):
                    ID = map(itemgetter(1), g)
                    # if debug: print 'Index bound: ' + str(ID[0]) + '-' + str(ID[-1]+1)
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

        return

class _load_var:
    """
    **'Variables' subset in FVCOM class contains the numpy arrays**

    Some variables are directly passed on from FVCOM output: ::

                       _el = elevation (m), 2D array (ntime, nnode)
                      |_julianTime = julian date, 1D array (ntime)
                      |_matlabTime = matlab time, 1D array (ntime)
                      |_tauc = bottom shear stress (m2/s2),
                      |      2D array (ntime, nele)
                      |_ua = depth averaged u velocity component (m/s),
                      |      2D array (ntime, nele)
                      |_va = depth averaged v velocity component (m/s),
     FVCOM.Variables._|      2D array (ntime, nele)
                      |_u = u velocity component (m/s),
                      |     3D array (ntime, nlevel, nele)
                      |_v = v velocity component (m/s),
                      |     3D array (ntime, nlevel, nele)
                      |_w = w velocity component (m/s),
                      |     3D array (ntime, nlevel, nele)

    Some others shall be generated as methods are being called, ex: ::
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
        self._3D = False
        self._opendap = type(data.variables).__name__=='DatasetType'

        # Pointer to History
        setattr(self, '_History', History)

        # Parallel computing attributs
        #self._cpus = mp.cpu_count()

        #List of keywords
        kwl2D = ['ua', 'va', 'zeta','depth_av_flow_dir', 'hori_velo_norm',
                 'depth_av_vorticity', 'depth_av_power_density',
                 'depth_av_power_assessment', 'tauc']
        kwl3D = ['ww', 'u', 'v', 'gls', 'tke', 'flow_dir', 'velo_norm',
                 'verti_shear', 'vorticity', 'power_density']
        #List of aliaSes
        al2D = ['ua', 'va', 'el','depth_av_flow_dir', 'hori_velo_norm',
               'depth_av_vorticity', 'depth_av_power_density',
               'depth_av_power_assessment', 'tauc']
        al3D = ['w', 'u', 'v', 'gls', 'tke', 'flow_dir', 'velo_norm',
                 'verti_shear', 'vorticity', 'power_density']

        # Figure out which quantity to treat
        self._kwl2D = []
        self._al2D = []
        for key, aliaS in zip(kwl2D, al2D):
            if key in data.variables.keys():
                self._kwl2D.append(key)
                self._al2D.append(aliaS)
            else:
                if debug: print key, " is missing !"


        self._kwl3D = []
        self._al3D = []
        for key, aliaS in zip(kwl3D, al3D):
            if key in data.variables.keys():
                self._kwl3D.append(key)
                self._al3D.append(aliaS)
            else:
                if debug: print key, " is missing !"

        if not len(self._kwl3D)==0: self._3D = True

        #Loading time stamps
        try:
            julianTime = data.variables['julianTime']
        except KeyError:
            # exception due to Save_as(netcdf)
            julianTime=data.variables['time']
            #  Work out if time is in julian time
            timeFlag = 0.0
            try:
                for key in julianTime.attributes:
                    if "julian" in julianTime.attributes[key].lower():
                        timeFlag += 1.0
            except AttributeError:  # pydap lib error
                if "julian" in julianTime.format.lower():
                    timeFlag += 1.0
            # if not julian time, convert in days in needed
            if timeFlag == 0.0:
                try:
                    for key in julianTime.attributes:
                        if "second" in julianTime.attributes[key].lower():
                            timeFlag += 1.0
                except AttributeError:  # pydap lib error
                    if "second" in julianTime.units.lower():
                        timeFlag += 1.0
                if not timeFlag == 0.0:  # TR: this conversion needs to be improved by introducing the right epoch
                    dayTime = julianTime[:] / (24*60*60) # convert in days
                    julianTime = dayTime

        if tx==[]:
            # get time and adjust it to matlab datenum
            self.julianTime = julianTime[:]
            self.matlabTime = self.julianTime[:] + 678942.0
            #-Append message to History field
            start = mattime_to_datetime(self.matlabTime[0])
            end = mattime_to_datetime(self.matlabTime[-1])
            text = 'Full temporal domain from ' + str(start) +\
                   ' to ' + str(end)
            self._History.append(text)
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]
            if debug: print 'Full temporal domain'
        else:
            #Time period           
            region_t = self._t_region(tx, julianTime, debug=debug)
            self._region_time = region_t
            ts = self._region_time[0]
            te = self._region_time[-1] + 1
            # get time and adjust it to matlab datenum
            self.julianTime = julianTime[ts:te]
            self.matlabTime = self.julianTime + 678942.0
            #-Append message to History field
            start = mattime_to_datetime(self.matlabTime[0])
            end = mattime_to_datetime(self.matlabTime[-1])
            text = 'Temporal domain from ' + str(start) +\
                   ' to ' + str(end)
            self._History.append(text)
            #Add time dimension to grid variables
            grid.ntime = self.julianTime.shape[0]
            if debug: print "ntime: ", grid.ntime
            if debug: print "region_t shape: ", region_t.shape

        # Define which loading function to use
        if grid._ax==[] and tx==[]:
            loadVar = self._load_full_time_full_region
            for key, aliaS in zip(self._kwl2D, self._al2D):
                loadVar(data, key, aliaS, debug=debug)
            for key, aliaS in zip(self._kwl3D, self._al3D):
                loadVar(data, key, aliaS, debug=debug)
        else:

            if grid._ax!=[] and tx!=[]:
                loadVar = self._load_partial_time_partial_region
            elif grid._ax==[] and tx!=[]:
                loadVar = self._load_partial_time_full_region
            else:
                loadVar = self._load_full_time_partial_region

            # Loading 2D variables
            for key, aliaS in zip(self._kwl2D, self._al2D):
                gc.collect()  # force garbage collector in an attempt to free up some RAM
                loadVar(data, grid, key, aliaS, debug=debug)

            # Loading 3D variables
            for key, aliaS in zip(self._kwl3D, self._al3D):
                gc.collect()  # force garbage collector in an attempt to free up some RAM
                loadVar(data, grid, key, aliaS, debug=debug)

            ##-------Parallelized loading block-------
            # if debug: startT = time.time()
            #
            # divisor = len(self._kwl2D)//self._cpus
            # remainder = len(self._kwl2D)%self._cpus
            #
            # if debug: print "Parallel loading 2D vars..."
            # if debug: start2D = time.time()
            #
            # for i in range(divisor):
            #     start = self._cpus * i
            #     end = start + (self._cpus-1)
            #     processes = [mp.Process(target=loadVar, args=(data, grid, key, aliaS, debug))\
            #                  for key, aliaS in zip(self._kwl2D[start:end], self._al2D[start:end])]
            #     # Run processes
            #     for p in processes:
            #         p.start()
            #     # Exit the completed processes
            #     for p in processes:
            #         p.join()
            #
            # # Remaining vars
            # if remainder != 0:
            #     start = int(-1 * remainder)
            #     processes = [mp.Process(target=loadVar, args=(data, grid, key, aliaS, debug))\
            #                  for key, aliaS in zip(self._kwl2D[start:], self._al2D[start:])]
            #     # Run processes
            #     for p in processes:
            #         p.start()
            #     # Exit the completed processes
            #     for p in processes:
            #         p.join()
            #
            # if debug: end2D = time.time()
            # if debug: print "...processing time: ", (end2D - start2D)
            #
            # if debug: print "Parallel loading 3D vars..."
            # if debug: start3D = time.time()
            #
            # for i in range(divisor):
            #     start = self._cpus * i
            #     end = start + (self._cpus-1)
            #     processes = [mp.Process(target=loadVar, args=(data, grid,key, aliaS, debug))\
            #                  for key, aliaS in zip(self._kwl3D[start:end], self._al3D[start:end])]
            #     # Run processes
            #     for p in processes:
            #         p.start()
            #     # Exit the completed processes
            #     for p in processes:
            #         p.join()
            #
            # # Remaining vars
            # if remainder != 0:
            #     start = int(-1 * remainder)
            #     processes = [mp.Process(target=loadVar, args=(data, grid, key, aliaS, debug))\
            #                  for key, aliaS in zip(self._kwl3D[start:], self._al3D[start:])]
            #     # Run processes
            #     for p in processes:
            #         p.start()
            #     # Exit the completed processes
            #     for p in processes:
            #         p.join()
            #
            # if debug: end3D = time.time()
            # if debug: endT = time.time()
            # if debug: print "...-loading 3D- processing time: ", (end3D - start3D)
            # if debug: print "...-loading 2D & 3D- processing time: ", (endT - startT)
            # #-------end-------

        if debug: print '...Passed'

        # Define method loadvar for use in import
        self._loadVar = loadVar

        return

    def _load_full_time_full_region(self, data, key, aliaS, debug=False):
        """
        loading variables for full time and space domains

        Inputs:
          - key = FVCOM variable name, str
          - aliaS = PySeidon variable alias, str

        Options:
          - debug = debug flag, boolean
        """
        if debug: print "loading " + str(aliaS) +"..."
        try:
            setattr(self, aliaS, data.variables[key].data)
        except AttributeError: #exeception due nc.Dataset
            setattr(self, aliaS, data.variables[key])

    def _load_partial_time_partial_region(self, data, grid, key, aliaS, debug=False):
        """
        loading variables for partial time and space domains

        Inputs:
          - key = FVCOM variable name, str
          - aliaS = PySeidon variable alias, str

        Options:
          - debug = debug flag, boolean
        """
        if debug: print "loading " + str(aliaS) +"..."

        # define time bounds
        ts = self._region_time[0]
        te = self._region_time[-1] + 1

        if key == 'zeta':
            region = grid._node_index
            horiDim = grid.nnode
        else:
            region = grid._element_index
            horiDim = grid.nele

        if key == 'verti_shear':
            vertiDim = grid.nlevel-1
        else:
            vertiDim = grid.nlevel

        # Find out if using netCDF4 or scipy
        try:
            Test = data.variables[key].data
            self._scipynetcdf = True
        except AttributeError: # exeception due nc.Dataset
            self._scipynetcdf = False

        if self._opendap:
        # loop over contiguous indexes for opendap
            H = 0 #local counter
            for k, g in groupby(enumerate(region), lambda (i,x):i-x):
                ID = map(itemgetter(1), g)
                #if debug: print 'Index bound: ' + str(ID[0]) + '-' + str(ID[-1]+1)
                if key in self._kwl2D:
                    if self._scipynetcdf:
                        #TR : Don't I need to transpose here?
                        var = data.variables[key].data[ts:te,ID[0]:(ID[-1]+1)]
                    else:
                        var = data.variables[key][ts:te,ID[0]:(ID[-1]+1)]
                    if H == 0:
                        setattr(self, aliaS,var)
                        H = 1
                    else:
                        setattr(self, aliaS, np.hstack((getattr(self, aliaS), var)))
                else:
                    if self._scipynetcdf:
                        #TR : Don't I need to transpose here?
                        var = data.variables[key].data[ts:te,:,ID[0]:(ID[-1]+1)]
                    else:
                        var = data.variables[key][ts:te,:,ID[0]:(ID[-1]+1)]
                    if H == 0:
                        setattr(self, aliaS,var)
                        H = 1
                    else:
                        setattr(self, aliaS, np.dstack((getattr(self, aliaS), var)))
        # TR comment: looping on time indices is a trick from Mitchell O'Flaherty-Sproul to improve loading time
        else:
            I = 0
            if key in self._kwl2D:
                setattr(self, aliaS, np.zeros((grid.ntime, horiDim)))
                for i in self._region_time:
                    if self._scipynetcdf:
                        getattr(self, aliaS)[I,:] = np.transpose(data.variables[key].data[i, region])
                    else:
                        getattr(self, aliaS)[I,:] = (data.variables[key][i, region])
                    I += 1
            else:
                setattr(self, aliaS, np.zeros((grid.ntime, vertiDim, horiDim)))
                for i in self._region_time:
                    if self._scipynetcdf:
                        getattr(self, aliaS)[I,:,:] = np.transpose(data.variables[key].data[i, :, region])
                    else:
                        getattr(self, aliaS)[I,:,:] = (data.variables[key][i, :, region])
                    I += 1


    def _load_full_time_partial_region(self, data, grid, key, aliaS, debug=False):
        """
        loading variables for full time domain and partial space domain

        Inputs:
          - key = FVCOM variable name, str
          - aliaS = PySeidon variable alias, str

        Options:
          - debug = debug flag, boolean
        """
        if debug: print "loading " + str(aliaS) +"..."
        if key == 'zeta':
            region = grid._node_index
            horiDim = grid.nnode
        else:
            region = grid._element_index
            horiDim = grid.nele

        if key == 'verti_shear':
            vertiDim = grid.nlevel-1
        else:
            vertiDim = grid.nlevel

        # Find out if using netCDF4 or scipy
        try:
            Test = data.variables[key].data
            self._scipynetcdf = True
        except AttributeError: # exeception due nc.Dataset
            self._scipynetcdf = False

        if self._opendap:
        # loop over contiguous indexes for opendap
            H = 0 #local counter
            for k, g in groupby(enumerate(region), lambda (i,x):i-x):
                ID = map(itemgetter(1), g)
                #if debug: print 'Index bound: ' + str(ID[0]) + '-' + str(ID[-1]+1)
                if key in self._kwl2D:
                    if self._scipynetcdf:
                        #TR : Don't I need to transpose here?
                        var = data.variables[key].data[:,ID[0]:(ID[-1]+1)]
                    else:
                        var = data.variables[key][:,ID[0]:(ID[-1]+1)]
                    if H == 0:
                        setattr(self, aliaS,var)
                        H = 1
                    else:
                        setattr(self, aliaS, np.hstack((getattr(self, aliaS), var)))
                else:
                    if self._scipynetcdf:
                        #TR : Don't I need to transpose here?
                        var = data.variables[key].data[:,:,ID[0]:(ID[-1]+1)]
                    else:
                        var = data.variables[key][:,:,ID[0]:(ID[-1]+1)]
                    if H == 0:
                        setattr(self, aliaS,var)
                        H = 1
                    else:
                        setattr(self, aliaS, np.dstack((getattr(self, aliaS), var)))
        else:
            # TR comment: looping on time indices is a trick from Mitchell O'Flaherty-Sproul to improve loading time
            if key in self._kwl2D:
                setattr(self, aliaS, np.zeros((grid.ntime, horiDim)))
                for i in range(grid.ntime):
                    if self._scipynetcdf:
                        getattr(self, aliaS)[i,:] = np.transpose(data.variables[key].data[i, region])
                    else:
                        getattr(self, aliaS)[i,:] = (data.variables[key][i, region])
            else:
                setattr(self, aliaS, np.zeros((grid.ntime, vertiDim, horiDim)))
                for i in range(grid.ntime):
                    if self._scipynetcdf:
                        getattr(self, aliaS)[i,:,:] = np.transpose(data.variables[key].data[i, :, region])
                    else:
                        getattr(self, aliaS)[i,:,:] = (data.variables[key][i, :, region])

    def _load_partial_time_full_region(self, data, grid, key, aliaS, debug=False):
        """
        loading variables for partial time domain and full space domain

        Inputs:
          - key = FVCOM variable name, str
          - aliaS = PySeidon variable alias, str

        Options:
          - debug = debug flag, boolean
        """
        if debug: print "loading " + str(aliaS) +"..."

        # define time bounds
        ts = self._region_time[0]
        te = self._region_time[-1] + 1

        if key == 'zeta':
            horiDim = grid.nnode
        else:
            horiDim = grid.nele

        if key == 'verti_shear':
            vertiDim = grid.nlevel-1
        else:
            vertiDim = grid.nlevel

        # Find out if using netCDF4 or scipy
        try:
            Test = data.variables[key].data
            self._scipynetcdf = True
        except AttributeError: # exeception due nc.Dataset
            self._scipynetcdf = False

        if self._opendap:
            if key in self._kwl2D:
                if self._scipynetcdf:
                    var = data.variables[key].data[ts:te,:]
                else:
                    var = data.variables[key][ts:te,:]
            else:
                if self._scipynetcdf:
                    var = data.variables[key].data[ts:te,:,:]
                else:
                    var = data.variables[key][ts:te,:,:]
            setattr(self, aliaS,var)
        else:
            I = 0
            # TR comment: looping on time indices is a trick from Mitchell O'Flaherty-Sproul to improve loading time
            if key in self._kwl2D:
                setattr(self, aliaS, np.zeros((grid.ntime, horiDim)))
                for i in self._region_time:
                    if self._scipynetcdf:
                        getattr(self, aliaS)[I,:] = np.transpose(data.variables[key].data[i, :])
                    else:
                        getattr(self, aliaS)[I,:] = (data.variables[key][i, :])
                    I += 1
            else:
                setattr(self, aliaS, np.zeros((grid.ntime, vertiDim, horiDim)))
                for i in self._region_time:
                    if self._scipynetcdf:
                        getattr(self, aliaS)[I,:,:] = np.transpose(data.variables[key].data[i, :, :])
                    else:
                        getattr(self, aliaS)[I,:,:] = (data.variables[key][i, :, :])
                    I += 1

    def _t_region(self, tx, julianTime, debug=False):
        """Return time indices included in time period, aka tx"""
        debug = debug or self._debug      
        if debug: print 'Computing region_t...'
        start = datetime.datetime.strptime(tx[0], '%Y-%m-%d %H:%M:%S')
        end = datetime.datetime.strptime(tx[1], '%Y-%m-%d %H:%M:%S')
        region_t = time_to_index(start, end, julianTime[:], debug=debug)
        if debug: print '...Passed'
        print '-Now working in time box-'
        return region_t
