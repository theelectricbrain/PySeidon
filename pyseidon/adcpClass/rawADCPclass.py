from __future__ import division
import numpy as np
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr
#from shortest_element_path import shortest_element_path
#import matplotlib.pyplot as plt
#import matplotlib.tri as Tri
#import matplotlib.ticker as ticker
#import seaborn
import scipy.io as sio
import h5py
from os import path

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


class rawADCP:
    def __init__(self, filename):
        self.QC = ['raw data']
        self.load(filename)
        self.Params_Stn4_SWNSreport(filename)
        self.load_rbrdata()

        ## set options
        self.options = {}
        self.options['showPA'] = 1
        self.options['showRBRavg'] = 1

        ## save a flow file in BPformat
        #save_FlowFile_BPFormat(fileinfo,adcp,rbr,saveparams,options)

    def load(self, filename):

        try:
            self.mat = sio.loadmat(filename,
                                struct_as_record=False, squeeze_me=True)

            self.adcp = self.mat['adcp']

        except NotImplementedError:
            self.mat = h5py.File(filename)
            self.adcp = self.mat['adcp']
            #self.adcp = Struct(**self.mat['adcp'])

    def Params_Stn4_SWNSreport(self, filename):
        fname = filename.split('/')
        filebase = fname[-1].split('_')[0]
        self.fileinfo = {}
        self.fileinfo['datadir'] = path.join(*fname[:-1]) + '/'
        self.fileinfo['ADCP'] = filebase + '_raw'
        self.fileinfo['outdir'] = path.join(*fname[:-1]) + '/'
        self.fileinfo['flowfile'] = filebase + '_Flow'
        self.fileinfo['rbr']= 'station4_grandPassageII_RBRSN_011857.mat'
        self.fileinfo['paramfile']= 'Params_Stn4_SWNSreport'

        #%% ADCP parameters
        self.saveparams = {}
        self.saveparams['tmin'] = 209
        self.saveparams['tmax'] = 240
        self.saveparams['zmin'] = 0
        self.saveparams['zmax'] = 20
        self.saveparams['approxdepth'] = 15.5
        self.saveparams['flooddir'] = 0
        self.saveparams['declination'] = -17.25
        self.saveparams['lat'] =  44.2605
        self.saveparams['lon'] = -66.3354
        self.saveparams['dabADCP'] = 0.5
        self.saveparams['dabPS'] = -0.6
        self.saveparams['rbr_hr_offset'] = 3

    def load_rbrdata(self):
        rbrFile = self.fileinfo['datadir'] + self.fileinfo['rbr']

        try:
            rbrMat = sio.loadmat(rbrFile,
                                struct_as_record=False, squeeze_me=True)

        except NotImplementedError:
            rbrMat = h5py.File(rbrFile)

        rbr = rbrMat['rbr']
        rbrout = {}
        rbrout['mtime'] = rbr.yd

        rbrout['temp'] = rbr.temperature
        rbrout['pres'] = rbr.pressure
        rbrout['depth'] = rbr.depth
        rbrout['mtime'] = rbr.yd
        self.rbr = rbrout





if __name__ == '__main__':
    #filename = 'GP-120726-BPd_raw.mat'
    filename = '140703-EcoEII_database/data/GP-120726-BPd_raw.mat'
    data = rawADCP(filename)




#stn = 'GP-120726-BPd';
#%% File information
#fileinfo.datadir = 	'../data/'; 	 %path to raw data files
#fileinfo.ADCP = [stn '_raw']; 	 %name of ADCP file
#fileinfo.outdir = '../data/'; 	 %path to output directory
#fileinfo.flowfile =  [stn,'_Flow']; 	 %name of output file with Flow data
#fileinfo.rbr = ['station4_grandPassageII_RBRSN_011857.mat'];
#fileinfo.paramfile = mfilename;
#
#%% ADCP parameters
#saveparams.tmin = 209; 	 %tmin (year day)
#saveparams.tmax = 240; 	 %tmax (year day)
#saveparams.zmin = 0;     %minimum z to include in saves file
#saveparams.zmax = 20;
#saveparams.approxdepth = 15.5; %Approximate depth
#saveparams.flooddir= 0; 	 %Flood direction (relative to true north, CW is positive)
#saveparams.declination = -17.25;%Declination angle
#saveparams.lat =  44.2605; 	 %latitude
#saveparams.lon = -66.3354; 	 %longitude
#saveparams.dabADCP = 0.5; 	 %depth above bottom of ADCP
#saveparams.dabPS = -0.6; 	 %depth above bottom of pressure sensor
#saveparams.rbr_hr_offset = 3;    % hour offset to convert rbr time to UTC
