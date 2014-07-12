#!/usr/bin/python2.7
# encoding: utf-8

class _load_var:
    """
'Variables' subset in FVCOM class contains the following numpy arrays:
    x = 
    y = 
    xc = 
    yc = 
    lon = 
    lat =
    lonc = 
    latc =
    nv = 
    h = bathymetry
    nbe = 
    a1u = 
    a2u = 
    aw0 = 
    awx = 
    awy =  
The other variables possess in-build set of descriptors, ex:
         _long_name = 'Vertically Averaged x-velocity'
        |_units = 'meters s-1'
    ua._|_grid = 'fvcom_grid'
        |_...
             
    """
    def __init__(self, data, debug=False):
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
        self.a1u = data.variables['a1u'][:]
        self.a2u = data.variables['a2u'][:]
        self.aw0 = data.variables['aw0'][:]
        self.awx = data.variables['awx'][:]
        self.awy = data.variables['awy'][:]

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
            #invisible variables
            self._D3 = True

        except KeyError:
            self.ua = data.variables['ua']
            self.va = data.variables['va']
            self._D3 = False

class _load_grid:
    '''
'Grid' subset in FVCOM class contains the following numpy arrays:
    trinodes =
    siglay =
    siglev = 
    nele = 
    node = 
    '''
    def __init__(self, data, debug=False):
        self.trinodes = data.variables['nv'][:]
        self.siglay = data.variables['siglay'][:]
        self.siglev = data.variables['siglev'][:]

        # Need to use len to get size of dimensions
        self.nele = data.dimensions['nele']
        self.node = data.dimensions['node']

