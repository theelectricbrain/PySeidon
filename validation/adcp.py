#This program combines all ADCP files into one

import os
from pyseidon_dvt import *
import numpy as np
from utide import  *

#Initialize ADCP data

dir = '/EcoII/acadia_uni/workspace/observed/GP/ADCP/GPcabled/'
pile = ""
list = os.listdir(dir)
list.sort()
listADCP = []
count = 0
for i in list:
	if ".mat" in i:
		listADCP.append(pile+i)
		count = count+1
print "There are",count,"adcp files in this directory"

adcp1=ADCP(dir+listADCP[0])
adcp2=ADCP(dir+listADCP[1])
adcp3=ADCP(dir+listADCP[2])
adcp4=ADCP(dir+listADCP[3])
adcp5=ADCP(dir+listADCP[4])
adcp6=ADCP(dir+listADCP[5])
adcp7=ADCP(dir+listADCP[6])
adcp8=ADCP(dir+listADCP[7])
adcp9=ADCP(dir+listADCP[8])
adcp10=ADCP(dir+listADCP[9])
adcp11=ADCP(dir+listADCP[10])
adcp12=ADCP(dir+listADCP[11])
adcp13=ADCP(dir+listADCP[12])

#This does nothing
for index in range(len(listADCP)):
	adcps=ADCP(dir+listADCP[index])

#This combines all adcp elevations
adcps.Variables.el=np.hstack((adcp1.Variables.el,adcp2.Variables.el,adcp3.Variables.el,adcp4.Variables.el,adcp5.Variables.el,adcp6.Variables.el,adcp7.Variables.el,adcp8.Variables.el,adcp9.Variables.el,adcp10.Variables.el,adcp11.Variables.el,adcp12.Variables.el,adcp13.Variables.el))
print "elevations combined"

#This combines all adcp matlab Times	
adcps.Variables.matlabTime=np.hstack((adcp1.Variables.matlabTime,adcp2.Variables.matlabTime,adcp3.Variables.matlabTime,adcp4.Variables.matlabTime,adcp5.Variables.matlabTime,adcp6.Variables.matlabTime,adcp7.Variables.matlabTime,adcp8.Variables.matlabTime,adcp9.Variables.matlabTime,adcp10.Variables.matlabTime,adcp11.Variables.matlabTime,adcp12.Variables.matlabTime,adcp13.Variables.matlabTime))
print "matlabTime combined"

#This combines all adcp surfs
adcps.Variables.surf=np.hstack((adcp1.Variables.surf,adcp2.Variables.surf,adcp3.Variables.surf,adcp4.Variables.surf,adcp5.Variables.surf,adcp6.Variables.surf,adcp7.Variables.surf,adcp8.Variables.surf,adcp9.Variables.surf,adcp10.Variables.surf,adcp11.Variables.surf,adcp12.Variables.surf,adcp13.Variables.surf))
print "surf combined"

#This combines all adcp pressures
adcps.Variables.pressure=np.hstack((adcp1.Variables.pressure,adcp2.Variables.pressure,adcp3.Variables.pressure,adcp4.Variables.pressure,adcp5.Variables.pressure,adcp6.Variables.pressure,adcp7.Variables.pressure,adcp8.Variables.pressure,adcp9.Variables.pressure,adcp10.Variables.pressure,adcp11.Variables.pressure,adcp12.Variables.pressure,adcp13.Variables.pressure))
print "pressure combined"

#This combines all adcp v velocities
adcps.Variables.v=np.vstack((adcp1.Variables.v,adcp2.Variables.v,adcp3.Variables.v,adcp4.Variables.v,adcp5.Variables.v,adcp6.Variables.v,adcp7.Variables.v,adcp8.Variables.v,adcp9.Variables.v,adcp10.Variables.v,adcp11.Variables.v,adcp12.Variables.v,adcp13.Variables.v))
print "v velocities combined"

#This combines all adcp u velocites
adcps.Variables.u=np.vstack((adcp1.Variables.u,adcp2.Variables.u,adcp3.Variables.u,adcp4.Variables.u,adcp5.Variables.u,adcp6.Variables.u,adcp7.Variables.u,adcp8.Variables.u,adcp9.Variables.u,adcp10.Variables.u,adcp11.Variables.u,adcp12.Variables.u,adcp13.Variables.u))
print "u velocities combined"

#Perform Harmonic Analysis of velocities
harmo=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False)

#Reconstruct velocities 
velos=adcps.Utils.Harmonic_reconstruction(harmo)

#Plot time el
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.el)

#Plot times u vel
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.u)

#Plot times v vel
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.v)

#Plot time surf
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.surf)

#This combines all adcp va velocities
adcps.Variables.va=np.hstack((adcp1.Variables.va,adcp2.Variables.va,adcp3.Variables.va,adcp4.Variables.va,adcp5.Variables.va,adcp6.Variables.va,adcp7.Variables.va,adcp8.Variables.va,adcp9.Variables.va,adcp10.Variables.va,adcp11.Variables.va,adcp12.Variables.va,adcp13.Variables.va))
print "va velocities combined"

#This combines all adcp ua velocites
adcps.Variables.ua=np.hstack((adcp1.Variables.ua,adcp2.Variables.ua,adcp3.Variables.ua,adcp4.Variables.ua,adcp5.Variables.ua,adcp6.Variables.ua,adcp7.Variables.ua,adcp8.Variables.ua,adcp9.Variables.ua,adcp10.Variables.ua,adcp11.Variables.ua,adcp12.Variables.ua,adcp13.Variables.ua))
print "ua velocities combined"

harmo2=solve(adcps.Variables.matlabTime,adcps.Variables.ua,adcps.Variables.va,adcps.Variables.lat)


