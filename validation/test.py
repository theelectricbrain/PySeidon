#This program simply combines all adcp files into one

import os
from pyseidon_dvt import *
import numpy as np
from utide import  *
import scipy.interpolate


#Initialize ADCP data

dir = '/EcoII/acadia_uni/workspace/observed/GP/ADCP/'
list = os.listdir(dir)
list.sort()
listADCP = []
 
for i in list:
        if ".mat" in i:
                listADCP.append(ADCP(dir+i))
                
for index in range(len(listADCP)):
	adcps=listADCP[index]

#Perform Harmonic Analysis of velocities
harmo=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False)

harmo1=adcps.Utils.Harmonic_analysis(elevation=False,velocity=True)

#Reconstruct velocities
velos=adcps.Utils.Harmonic_reconstruction(harmo)

velos1=adcps.Utils.Harmonic_reconstruction(harmo1)


#Plot #1: time el"
#np.polyfit(adcps.Variables.matlabTime,adcps.Variables.el,4)
#scipy.interpolate.interp1d(adcps.Variables.matlabTime,adcps.Variables.el,kind='cubic')


adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.el)
adcps.Plots.Histogram(adcps.Variables.matlabTime,adcps.Variables.el)
#Plot #2: times u vel"
adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.u)

#Plot #3:  times v vel"
adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.v)

#Plot #4:  time surf"
adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.surf)

#Plot #5:  time va vel"
adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.va)

#Plot #6:  times ua vel"
adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.ua)









