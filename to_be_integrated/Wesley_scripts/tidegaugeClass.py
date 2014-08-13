import numpy as np
import scipy.io as sio
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr


class Tidegauge:
    def __init__(self, filename):
        self.load(filename)

    def load(self, filename):
        print 'Loading...'

        self.mat = sio.loadmat(filename,
                               struct_as_record=False, squeeze_me=True)

        self.RBR = self.mat['RBR']
        self.data = self.RBR.data
        self.time = self.RBR.date_num_Z
        self.lat = self.RBR.lat
        self.lon = self.RBR.lon

        self.elev = self.data - np.mean(self.data)

        print 'Done'

    def harmonics(self, **kwarg):

        self.coef = ut_solv(self.time,
                            self.elev, [],
                            self.lat, **kwarg)

    def reconstr(self, time):
        self.ts_recon, _ = ut_reconstr(time, self.coef)


if __name__ == '__main__':
    filename = 'Westport_015892_20140325_1212_Z.mat'
    tide = Tidegauge(filename)

    ut_constits = ['M2','S2','N2','K2','K1','O1','P1','Q1']
    tide.harmonics(method='ols', cnstit=ut_constits, ordercnstit='frq',
                   nodiagn=True)

    tide.reconstr(tide.time)
