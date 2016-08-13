################################################################################################################################################

# This program takes ADCP data with gaps and bad data points and attempts to rid both gaps and bad data points

# Table of Contents:

# 1. Elevation Data

# 2. Ua Velocity

# 3. Va Velocity

################################################################################################################################################

# Combining ADCP files

from pylab import *
import os
from pyseidon_dvt import *
import numpy as np
from utide import  *
from scipy.interpolate import interp1d
import csv
import matplotlib.pyplot as plt 
import copy as cp

# Initialize ADCP data

dir = '/EcoII/acadia_uni/workspace/observed/GP/ADCP/GPcabled/'
list = os.listdir(dir)
list.sort()
listADCP = []

for i in list:
        if ".mat" in i:	
		 listADCP.append(dir+i)

for index in range(len(listADCP)):
	adcp=ADCP(listADCP[index])
	
	if index == 0:
		adcps=cp.deepcopy(adcp)
		adcp12=cp.deepcopy(adcp)
		         
	elif index == 1:
		adcp12.Variables.el=np.hstack((adcp12.Variables.el,adcp.Variables.el))
		adcp12.Variables.ua=np.hstack((adcp12.Variables.ua,adcp.Variables.ua))
		adcp12.Variables.va=np.hstack((adcp12.Variables.va,adcp.Variables.va))
		adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))
	
	if index != 0:
		adcps.Variables.el=np.hstack((adcps.Variables.el,adcp.Variables.el))
		adcps.Variables.ua=np.hstack((adcps.Variables.ua,adcp.Variables.ua))
		adcps.Variables.va=np.hstack((adcps.Variables.va,adcp.Variables.va))
		adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))
		
#newadcp.Variables.el=np.stack((adcp.Variables.el,newadcp.Variables.el))
#newadcp.Variables.ua=np.hstack((adcp.Variables.ua,newadcp.Variables.ua))
#newadcp.Variables.va=np.hstack((adcp.Variables.va,newadcp.Variables.va))
modifiedtime=np.arange(adcps.Variables.matlabTime[0],adcps.Variables.matlabTime[-1],(adcps.Variables.matlabTime[1]-adcps.Variables.matlabTime[0]))
newadcp=cp.deepcopy(adcp)               
newadcp.Variables.matlabTime=modifiedtime

################################################################################################################################################

# 1. Elevation 


# Perform Harmonic Analysis of Elevation
harmo1=adcp12.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)

# Reconstruct with original measured data
velos1=adcps.Utils.Harmonic_reconstruction(harmo1)

# Difference
diff_el = abs(adcps.Variables.el - velos1['h'])

harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
velos2=adcps.Utils.Harmonic_reconstruction(harmo2)

# The following while loop takes ~50 mins to run
while (diff_el > (np.mean(diff_el)+2.5*np.std(diff_el))).any():
	harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
	velos2=adcps.Utils.Harmonic_reconstruction(harmo2)
	diff_el = abs(adcps.Variables.el - velos2['h'])
	x = np.where(diff_el>(np.mean(diff_el)+2.5*np.std(diff_el)))
	diff_el=np.delete(diff_el,x)
        adcps.Variables.el = np.delete(adcps.Variables.el,x)
        adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,x)
	print shape(diff_el)
	print shape(x)

harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
velos2=adcps.Utils.Harmonic_reconstruction(harmo2)
diff_el=abs(adcps.Variables.el-velos2['h'])

#plt.plot(adcps.Variables.matlabTime,diff_el,'.')
#show()

#reconstruct1 = newadcp.Utils.Harmonic_reconstruction(harmo2)

#for i in diff_el:
#	residual[i]=(diff_el-velos2['h'])/(.5*abs(diff_el)+abs(velos2['h']))

#mr=np.mean(residual)
#sdr=np.std(residual)

#plt.plot()
#show()
#plt.plot(modifedtime,velos2['h'])
#show()

#velos=newadcp.Utils.Harmonic_reconstruction(harmo2)
#newadcp.Variables.el=velos['h']
#plt.plot(adcps.Variables.matlabTime,adcps.Variables.el,'.')

#####################################################################################################################################


# 2. Ua Velocity

harmo3=adcp12.Utils.Harmonic_analysis(elevation=False,velocity=True)
velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
diff_ua = abs(adcps.Variables.ua - velos3['u'])
#acdps.Variables.ua has bad data taken out from previous loop

for index in range(len(listADCP)):
        adcp=ADCP(listADCP[index])
        
	if index == 0:
                adcps=cp.deepcopy(adcp)
                adcp12=cp.deepcopy(adcp)
	        
	else index == 1:
                adcp12.Variables.el=np.hstack((adcp12.Variables.el,adcp.Variables.el))
               	adcp12.Variables.ua=np.hstack((adcp12.Variables.ua,adcp.Variables.ua))
               	adcp12.Variables.va=np.hstack((adcp12.Variables.va,adcp.Variables.va))
               	adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))

	if index != 0:
                adcps.Variables.el=np.hstack((adcps.Variables.el,adcp.Variables.el))
                adcps.Variables.ua=np.hstack((adcps.Variables.ua,adcp.Variables.ua))
                adcps.Variables.va=np.hstack((adcps.Variables.va,adcp.Variables.va))
                adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))


while (diff_ua > (np.mean(diff_ua)+2.5*np.std(diff_ua))).any():
        harmo3=(adcps.Utils.Harmonic_analysis(elevation=False,velocity=True))
        velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
        diff_ua = abs(adcps.Variables.ua - velos3['u'])
        x = np.where(diff_ua>(np.mean(diff_ua)+2.5*np.std(diff_ua)))
        diff_ua=np.delete(diff_ua,x)
        adcps.Variables.ua = np.delete(adcps.Variables.ua,x)
        adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,x)
       # print shape(diff_ua)
        print shape(x)

harmo3=adcps.Utils.Harmonic_analysis(elevation=False,velocity=True)
velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
diff_ua=abs(adcps.Variables.ua-velos3['u'])

#plt.plot(adcps.Variables.matlabTime,diff_ua,'.')
#show()

#################################################################################################################################################

# 3. Va Velocity

harmo4=adcp12.Utils.Harmonic_analysis(elevation=False,velocity=True)
velos4=adcps.Utils.Harmonic_reconstruction(harmo4)
diff_va = abs(adcps.Variables.va - velos4['v'])

for index in range(len(listADCP)):
        adcp=ADCP(listADCP[index])

        if index == 0:
                adcps=cp.deepcopy(adcp)
                adcp12=cp.deepcopy(adcp)

        else index == 1:
                adcp12.Variables.el=np.hstack((adcp12.Variables.el,adcp.Variables.el))
                adcp12.Variables.ua=np.hstack((adcp12.Variables.ua,adcp.Variables.ua))
                adcp12.Variables.va=np.hstack((adcp12.Variables.va,adcp.Variables.va))
                adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))

        if index != 0:
                adcps.Variables.el=np.hstack((adcps.Variables.el,adcp.Variables.el))
                adcps.Variables.ua=np.hstack((adcps.Variables.ua,adcp.Variables.ua))
                adcps.Variables.va=np.hstack((adcps.Variables.va,adcp.Variables.va))
                adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))

while (diff_va > (np.mean(diff_va)+2.5*np.std(diff_va))).any():
        harmo4=(adcps.Utils.Harmonic_analysis(elevation=False,velocity=True))
        velos4=adcps.Utils.Harmonic_reconstruction(harmo4)
        diff_va = abs(adcps.Variables.va - velos4['v'])
        x = np.where(diff_va>(np.mean(diff_va)+2.5*np.std(diff_va)))
        diff_va=np.delete(diff_va,x)
        adcps.Variables.va = np.delete(adcps.Variables.va,x)
        adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,x)
        print shape(diff_va)
        print shape(x)

harmo4=(adcps.Utils.Harmonic_analysis(elevation=False,velocity=True))
velos4=adcps.Utils.Harmonic_reconstruction(harmo4)
diff_va=abs(adcps.Variables.va-velos2['v'])

#plt.plot(adcps.Variables.matlabTime,diff_ua,'.')
#show()


################################################################################################################################################

