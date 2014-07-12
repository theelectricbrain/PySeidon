import numpy as np
#import scipy.io as sio
import h5py


class ADCP:
    ''' 
A class/structure for ADCP data.
  Only takes a file name as input, ex: testAdcp=ADCP('./pat_to_matlab_file/filename')
  Functionality strutured as follows:
               _Data. = raw matlab file data
              |_Variables. = useable adcp variables and quantities
              |_QC. = Quality Control metadata
    testAdcp._|_Utils. = set of useful functions
              |_Plots. = plotting functions
              |_method_1
              | ...      = methods and analysis techniques intrinsic to ADCPs
              |_method_n
    '''
    def __init__(self, filename, debug=False):
        ''' Initialize ADCP class. Notes: only handle raw ADCP matlab data at the mo.'''    
        self._debug = debug
        if debug:
            print '-Debug mode on-' 

        self.QC = ['Raw data']
        self.Data = h5py.File(filename)
        

    def _load(self, filename):
        #self.mat = sio.loadmat(filename,struct_as_record=False, squeeze_me=True)
        self.mat = h5py.File(filename)

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
