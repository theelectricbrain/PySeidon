import netCDF4 as nc


def createNC(data):

    ncFile = nc.Dataset('test.nc', 'w', format='NETCDF4')

    #ncgrp = ncFile.createGroup('regioned')
    #ncgrp.createDimension('dim', None)
    #ncgrp = ncFile.createGroup('regioned')
    #ncFile.createDimension('dimTest', None)

    ncFile.createDimension('dim', None)

    time = ncFile.createVariable('time', 'f8', ('dim',))
    time[:] = data['time']

    x = ncFile.createVariable('x', 'f8', ('dim',))
    x[:] = data['x']

    y = ncFile.createVariable('y', 'f8', ('dim',))
    y[:] = data['y']

    xc = ncFile.createVariable('xc', 'f8', ('dim',))
    xc[:] = data['xc']

    yc = ncFile.createVariable('yc', 'f8', ('dim',))
    yc[:] = data['yc']

    h = ncFile.createVariable('h', 'f8', ('dim',))
    h[:] = data['h']

    lon = ncFile.createVariable('lon', 'f8', ('dim',))
    lon[:] = data['lon']

    lat = ncFile.createVariable('lat', 'f8', ('dim',))
    lat[:] = data['lat']

    lonc = ncFile.createVariable('lonc', 'f8', ('dim',))
    lonc[:] = data['lonc']

    latc = ncFile.createVariable('latc', 'f8', ('dim',))
    latc[:] = data['latc']

    elev = ncFile.createVariable('elev', 'f8', ('dim', 'dim'))
    elev[:] = data['elev']

    ua = ncFile.createVariable('ua', 'f8', ('dim', 'dim'))
    ua[:] = data['ua']

    va = ncFile.createVariable('va', 'f8', ('dim', 'dim'))
    va[:] = data['va']

    node_index = ncFile.createVariable('node_index', 'f8', ('dim',))
    node_index[:] = data['node_index']

    element_index = ncFile.createVariable('element_index', 'f8', ('dim',))
    element_index[:] = data['element_index']

    nbe = ncFile.createVariable('nbe', 'f8', ('dim', 'dim'))
    nbe[:] = data['nbe']

    nv = ncFile.createVariable('nv', 'f8', ('dim', 'dim'))
    nv[:] = data['nv']

    a1u = ncFile.createVariable('a1u', 'f8', ('dim', 'dim'))
    a1u[:] = data['a1u']

    a2u = ncFile.createVariable('a2u', 'f8', ('dim', 'dim'))
    a2u[:] = data['a2u']

    aw0 = ncFile.createVariable('aw0', 'f8', ('dim', 'dim'))
    aw0[:] = data['aw0']

    awx = ncFile.createVariable('awx', 'f8', ('dim', 'dim'))
    awx[:] = data['awx']

    awy = ncFile.createVariable('awy', 'f8', ('dim', 'dim'))
    awy[:] = data['awy']

    siglay = ncFile.createVariable('siglay', 'f8', ('dim', 'dim'))
    siglay[:] = data['siglay']

    siglev = ncFile.createVariable('siglev', 'f8', ('dim', 'dim'))
    siglev[:] = data['siglev']

    ncFile.close()
