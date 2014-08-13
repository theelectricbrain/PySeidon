import numpy as np
import scipy.io as sio
import h5py

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

class ADCP:

    def __init__(self, filename, percent_of_depth=1):

        self.percent_of_depth = percent_of_depth
        self.load(filename)
        self.depthAverage()

    def __str__(self):
        '''
        What shows when
        print datais run
        '''
        return 'ADCP Class'

    def __repr__(self):
        return 'ADCP Class'

    def load(self, filename):
        try:
            self.mat = sio.loadmat(filename,
                                struct_as_record=False, squeeze_me=True)

            self.lat = self.mat['lat']
            self.lon = self.mat['lon']

            self.data = self.mat['data']
            self.bins = self.data.bins[:].flatten()
            self.north_vel = self.data.north_vel[:]
            self.east_vel = self.data.east_vel[:]
            self.vert_vel = self.data.vert_vel[:]
            self.dir_vel = self.data.dir_vel[:]
            self.mag_signed_vel = self.data.mag_signed_vel[:]
            self.ucross = self.data.ucross[:]
            self.ualong = self.data.ualong[:]

            self.pressure = self.mat['pres']
            self.surf = self.pressure.surf[:]

            self.time = self.mat['time']
            self.mtime = self.time.mtime[:]

        except NotImplementedError:
            self.mat = h5py.File(filename)

            self.lat = self.mat['lat'][:]
            self.lon = self.mat['lon'][:]

            #self.data = self.mat['data']
            self.data = Struct(**self.mat['data'])
            self.bins = self.data.bins[:].flatten()
            self.north_vel = self.data.north_vel[:].T
            self.east_vel = self.data.east_vel[:].T
            self.vert_vel = self.data.vert_vel[:].T
            self.dir_vel = self.data.dir_vel[:].T
            self.mag_signed_vel = self.data.mag_signed_vel[:].T
            self.ucross = self.data.Ucross[:].T
            self.ualong = self.data.Ualong[:].T

            self.pressure = Struct(**self.mat['pres'])
            self.surf = self.pressure.surf[:].flatten()

            self.time = Struct(**self.mat['time'])
            self.mtime = self.time.mtime[:].flatten()

    def depthAverage(self):
        '''Find the depth average of a variable based on percent_of_depth
        choosen by the user. Currently only working for east_vel (u) and
        north_vel (v) '''

        ind = np.argwhere(self.bins < self.percent_of_depth * self.surf[:,None])
        index = ind[np.r_[ind[1:,0] != ind[:-1,0], True]]
        data_ma_u = np.ma.array(self.east_vel, mask=np.arange(self.east_vel.shape[1]) > index[:, 1, None])
        data_ma_v = np.ma.array(self.north_vel, mask=np.arange(self.north_vel.shape[1]) > index[:, 1, None])
        self.ua = np.array(data_ma_u.mean(axis=1))
        self.va = np.array(data_ma_v.mean(axis=1))



if __name__ == '__main__':
    filename = 'Flow_GP-100915-BPa.mat'
    filename = '/home/wesley/github/karsten/project/140703-EcoEII_database/data/GP-120726-BPd_10minavg.mat'
    data = ADCP(filename)
