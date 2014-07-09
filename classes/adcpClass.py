import numpy as np
import scipy.io as sio


class ADCP:

    def __init__(self, filename):

        self.load(filename)

    def load(self, filename):
        self.mat = sio.loadmat(filename,
                               struct_as_record=False, squeeze_me=True)

        self.lat = self.mat['lat']
        self.lon = self.mat['lon']

        self.data = self.mat['data']
        self.north_vel = self.data.north_vel
        self.east_vel = self.data.east_vel
        self.vert_vel = self.data.vert_vel
        self.dir_vel = self.data.dir_vel
        self.mag_signed_vel = self.data.mag_signed_vel
        self.ucross = self.data.ucross
        self.ualong = self.data.ualong

        self.pressure = self.mat['pres']
        self.surf = self.pressure.surf

        self.time = self.mat['time']
        self.mtime = self.time.mtime


if __name__ == '__main__':
    filename = 'Flow_GP-100915-BPa.mat'
    data = ADCP(filename)
