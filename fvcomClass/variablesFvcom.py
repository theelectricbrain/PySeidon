#!/usr/bin/python2.7
# encoding: utf-8

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
    norma = depth averaged velocity norm
    norm = velocity norm
    ...            
    """
    def __init__(self, data, debug=False):
        if debug:
            print 'Loading variables...'
        # elev timeseries
        self.el = data.variables['zeta']
        # get time and adjust it to matlab datenum
        self.julianTime = data.variables['time'][:]
        self.matlabTime = self.julianTime + 678942
        try:
            self.ww = data.variables['ww']
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
        if debug:
            print '...Passed'

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
    def __init__(self, data, debug=False):
        if debug:
            print 'Loading grid...'
        self.x = data.variables['x'][:]
        self.y = data.variables['y'][:]
        self.xc = data.variables['xc'][:]
        self.yc = data.variables['yc'][:]
        self.lon = data.variables['lon'][:]
        self.lat = data.variables['lat'][:]
        self.lonc = data.variables['lonc'][:]
        self.latc = data.variables['latc'][:]
        self.nv = data.variables['nv'][:]
        self.h = data.variables['h'][:]
        self.nbe = data.variables['nbe'][:]
        self.trinodes = data.variables['nv'][:]
        self.siglay = data.variables['siglay'][:]
        self.siglev = data.variables['siglev'][:]
        #TR_comments: what the hell are the a* parameters???
        self.a1u = data.variables['a1u'][:]
        self.a2u = data.variables['a2u'][:]
        self.aw0 = data.variables['aw0'][:]
        self.awx = data.variables['awx'][:]
        self.awy = data.variables['awy'][:]
        # Need to use len to get size of dimensions
        self.nele = len(data.dimensions['nele'])
        self.node = len(data.dimensions['node'])
        # Custom grid variables
        if debug:
            print '...Passed'


