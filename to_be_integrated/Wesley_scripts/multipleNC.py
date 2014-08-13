import fnmatch
import os
import netCDF4 as nc
import numpy as np
import cPickle as pickle


def findFiles(filename, name):
    '''
    Wesley comment: the name needs to be a linux expression to find files you
    want. For multiple station files, this would work
    name = '*station*.nc'

    For just dngrid_0001 and no restart files:
    name = 'dngrid_0*.nc'
    will work
    '''

    matches = []
    for root, dirnames, filenames in os.walk(filename):
        for filename in fnmatch.filter(filenames, name):
            matches.append(os.path.join(root, filename))

    return matches


if __name__ == '__main__':
    filename = '/array/data1/073208o/workspace_matlab/runs/2013_run'

    name = 'dngrid_0*.nc'
    matches = findFiles(filename, name)

    time = np.array([])
    ua = np.array([])
    for i,v in enumerate(matches):
        data = nc.Dataset(v, 'r')
        t = data.variables['time'][:]
        u = data.variables['ua'][:, 0]

        time = np.hstack((time, t))
        ua = np.hstack((ua, u))

    test = {'time':time, 'ua':ua}
    pickle.dump(test, open('test.p', 'wb'))
