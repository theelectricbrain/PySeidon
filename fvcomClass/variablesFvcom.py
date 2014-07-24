#!/usr/bin/python2.7
# encoding: utf-8

from jdcal import gcal2jd
import numpy as np

class _load_var:
    """
'Variables' subset in FVCOM class contains the numpy arrays:
-----------------------------------------------------------
  Some variables are directly passed on from FVCOM output
  (i.e. el, julianTime, ww, u, v, ua, va) and possess in-build
  set of descriptors, ex:
         _long_name = 'Vertically Averaged x-velocity'
        |_units = 'meters s-1'
    ua._|_grid = 'fvcom_grid'
        |_...
  Some others shall be generated as methods are being called, ex:
    hori_velo_norm = horizontal velocity norm
    velo_norm = velocity norm
    verti_shear = vertical shear
    ...            
    """
    def __init__(self, data, grid, tx, QC, debug=False):
        #Pointer to QC
        self._QC = QC
        QC = self._QC

        if debug:
            print 'Loading variables...'
        #Check if time period defined
        self.julianTime = data.variables['time'][:]      
        if tx:
            #Time period           
            region_t = self._t_region(tx, debug=debug)
            #Quick reshape
            region_t = region_t.T[0,:]
            self._region_time = region_t
            # get time and adjust it to matlab datenum
            self.julianTime = data.variables['time'][region_t]
            self.matlabTime = self.julianTime + 678942

            #Check if bounding box has been defined
            if not hasattr(grid, '_region_e'):
                # elev timeseries
                self.el = data.variables['zeta'][region_t,:]           
                try:
                    self.ww = data.variables['ww'][region_t,:,:]
                    self.u = data.variables['u'][region_t,:,:]
                    self.v = data.variables['v'][region_t,:,:]
                    self.ua = data.variables['ua'][region_t,:]
                    self.va = data.variables['va'][region_t,:]
                    # invisible variables
                    self._3D = True
                except KeyError:
                    self.ua = data.variables['ua'][region_t,:]
                    self.va = data.variables['va'][region_t,:]
                    self._3D = False
            else:       
                #Bounding box
                region_e = grid._region_e
                region_n = grid._region_n
                #Redefine variables in bounding box & time period
                # elev timeseries
                self.el = data.variables['zeta'][region_t,region_n] 
                try:
                    self.ww = data.variables['ww'][region_t,:,region_e]
                    self.u = data.variables['u'][region_t,:,region_e]
                    self.v = data.variables['v'][region_t,:,region_e]
                    self.ua = data.variables['ua'][region_t,region_e]
                    self.va = data.variables['va'][region_t,region_e]
                    # invisible variables
                    self._3D = True
                except KeyError:
                    self.ua = data.variables['ua'][region_t,region_e]
                    self.va = data.variables['va'][region_t,region_e]
                    self._3D = False          
        else:
            #-Append message to QC field
            text = 'Full temporal domain'
            self._QC.append(text)
            # get time and adjust it to matlab datenum
            self.julianTime = data.variables['time'][:]
            self.matlabTime = self.julianTime + 678942

            #Check if bounding box has been defined
            if not hasattr(grid, '_region_e'):
                # elev timeseries
                self.el = data.variables['zeta'][:,:]           
                try:
                    self.ww = data.variables['ww'][:,:,:]
                    self.u = data.variables['u'][:,:,:]
                    self.v = data.variables['v'][:,:,:]
                    self.ua = data.variables['ua'][:,:]
                    self.va = data.variables['va'][:,:]
                    # invisible variables
                    self._3D = True
                except KeyError:
                    self.ua = data.variables['ua'][:,:]
                    self.va = data.variables['va'][:,:]
                    self._3D = False
            else:       
                #Bounding box
                region_e = grid._region_e
                region_n = grid._region_n
                #Redefine variables in bounding box & time period
                # elev timeseries
                self.el = data.variables['zeta'][:,region_n] 
                try:
                    self.ww = data.variables['ww'][:,:,region_e]
                    self.u = data.variables['u'][:,:,region_e]
                    self.v = data.variables['v'][:,:,region_e]
                    self.ua = data.variables['ua'][:,region_e]
                    self.va = data.variables['va'][:,region_e]
                    # invisible variables
                    self._3D = True
                except KeyError:
                    self.ua = data.variables['ua'][:,region_e]
                    self.va = data.variables['va'][:,region_e]
                    self._3D = False
        if debug:
            print '...Passed'

    def _t_region(self, tx, quiet=False, debug=False):
        '''Return time indexes included in time period, aka tx'''       
        if debug or self._debug:
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
        if debug or self._debug:
            print '...Passed'

        # Add metadata entry
        if not quiet:
            text = 'Time period =' + str(tx)
            self._QC.append(text)
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
    region_e = ???
    region_n = ???
    ax = bounding 
    ...     
    '''
    def __init__(self, data, ax, QC, debug=False):
        if debug:
            print 'Loading grid...'
        #Pointer to QC
        self._QC = QC
        QC = self._QC
        if not ax:
            #Define grid variables on the entire domain:
            self.lon = data.variables['lon'][:]
            self.lat = data.variables['lat'][:]
            self.lonc = data.variables['lonc'][:]
            self.latc = data.variables['latc'][:]
            self.x = data.variables['x'][:]
            self.y = data.variables['y'][:]
            self.xc = data.variables['xc'][:]
            self.yc = data.variables['yc'][:]
            self.h = data.variables['h'][:]
            self.siglay = data.variables['siglay'][:,:]
            self.siglev = data.variables['siglev'][:,:]
            #TR_comments: what the hell are the a* parameters???
            self.a1u = data.variables['a1u'][:,:]
            self.a2u = data.variables['a2u'][:,:]
            self.aw0 = data.variables['aw0'][:,:]
            self.awx = data.variables['awx'][:,:]
            self.awy = data.variables['awy'][:,:]
            self.nv = data.variables['nv'][:,:]
            self.nbe = data.variables['nbe'][:,:]
            #-Need to use len to get size of dimensions
            self.nele = len(data.dimensions['nele'])
            self.node = len(data.dimensions['node'])
            #-Append message to QC field
            text = 'Full spatial domain'
            self._QC.append(text)
        else:           
            #Define bounding box
            self.lon = data.variables['lon'][:]
            self.lat = data.variables['lat'][:]
            self.lonc = data.variables['lonc'][:]
            self.latc = data.variables['latc'][:]
            self.nv = data.variables['nv'][:,:]
            self.nbe = data.variables['nbe'][:,:]
            region_e, region_n = self._bounding_box(ax, debug=debug)
            #Quick reshape
            region_e = region_e.T[0,:]
            region_n = region_n.T[0,:]
            self._region_e = region_e
            self._region_n = region_n
            #TR: add check neighbouring nodes too and append to region_n
            #Re-indexing nv and nbe, i+1 because of python indexing
            #if debug:
            #    print "re-indexing"
            #j = 0
            #NV = np.empty(self.nbe.shape, dtype=int)
            #NBE = np.empty(self.nbe.shape, dtype=int) 
            #for i in region_n:
            #    ne = np.where(self.nv==i+1)
            #    NV[ne] = j
            #    j += 1
            #j =0
            #for i in region_n:
            #    no = np.where(self.nbe==i+1)
            #    NBE[no] = j
            #    j += 1
            #self.trinodes = NV.T
            #self.triele = NBE.T

            #Redefine variables in bounding box
            self.nele = len(region_e)
            self.node = len(region_n)
            self.lon = data.variables['lon'][region_n]
            self.lat = data.variables['lat'][region_n]
            self.lonc = data.variables['lonc'][region_e]
            self.latc = data.variables['latc'][region_e]
            self.x = data.variables['x'][region_n]
            self.y = data.variables['y'][region_n]
            self.xc = data.variables['xc'][region_e]
            self.yc = data.variables['yc'][region_e]
            self.h = data.variables['h'][region_n]
            self.siglay = data.variables['siglay'][:,region_n]
            self.siglev = data.variables['siglev'][:,region_n]
            #TR_comments: what the hell are the a* parameters???
            self.a1u = data.variables['a1u'][:,region_e]
            self.a2u = data.variables['a2u'][:,region_e]
            self.aw0 = data.variables['aw0'][:,region_e]
            self.awx = data.variables['awx'][:,region_e]
            self.awy = data.variables['awy'][:,region_e]
            self.nv = data.variables['nv'][:,region_e]
            self.nbe = data.variables['nbe'][:,region_e]
            #-Need to use len to get size of dimensions
            self.nele = self.h.shape[0]
            self.node = self.xc.shape[0]

            if debug:
                print '...Passed'

    def _ele_region(self, ax, debug=False):
        '''Return element indexes included in bounding box, aka ax'''       
        if debug:
            print 'Computing region_e...'

        region_e = np.argwhere((self.lonc >= ax[0]) &
                               (self.lonc <= ax[1]) &
                               (self.latc >= ax[2]) &
                               (self.latc <= ax[3]))          
        if debug or self._debug:
            print '...Passed'

        return region_e

    def _node_region(self, ax, debug=False):
        '''Return node indexes included in bounding box, aka ax'''
        if debug:
            print 'Computing region_n...'

        region_n = np.argwhere((self.lon >= ax[0]) &
                               (self.lon <= ax[1]) &
                               (self.lat >= ax[2]) &
                               (self.lat <= ax[3]))
        if debug or self._debug:
            print '...Passed'

        return region_n

    def _bounding_box(self, ax, quiet=False, debug=False):
        """
        Define bounding box and reset the box by default.
        Input ex:
        -------- 
          _bounding_box(ax=[min lon, max lon, min lat, max lat])
        """
        if not ax:
            region_e = range(self.nele)
            region_n = range(self.node)
        else:
            region_e = self._ele_region(ax, debug=debug)
            region_n = self._node_region(ax, debug=debug)

        # Add metadata entry
        if not quiet:
            text = 'Bounding box =' + str(ax)
            self._QC.append(text)
            print '-Now working in bounding box-'
        return region_e, region_n

