#This program simply combines all adcp files into one
from pylab import *
import os
from pyseidon_dvt import *
import numpy as np
from utide import  *
import scipy.interpolate
import csv
import matplotlib.pyplot as plt 
import copy as cp
#Initialize ADCP data

dir = '/EcoII/acadia_uni/workspace/observed/GP/ADCP/GPcabled/'
list = os.listdir(dir)
list.sort()
listADCP = []

for i in list:
        if ".mat" in i:	
		 listADCP.append(dir+i)

modifiedtime=np.arange(735840.74652777787,736227.5381944445,0.0034722221316779719)
for index in range(len(listADCP)):
	adcp=ADCP(listADCP[index])
	
	if index == 0:
		adcps=cp.copy(adcp)
		adcp12=cp.copy(adcp)
		         
	if index == 1:
		adcp12.Variables.el=np.hstack((adcp12.Variables.el,adcp.Variables.el))
		adcp12.Variables.ua=np.hstack((adcp12.Variables.ua,adcp.Variables.ua))
		adcp12.Variables.va=np.hstack((adcp12.Variables.va,adcp.Variables.va))
		adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))
	else:
		adcps.Variables.el=np.hstack((adcps.Variables.el,adcp.Variables.el))
		adcps.Variables.ua=np.hstack((adcps.Variables.ua,adcp.Variables.ua))
		adcps.Variables.va=np.hstack((adcps.Variables.va,adcp.Variables.va))
		adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))
		
#newadcp.Variables.el=np.stack((adcp.Variables.el,newadcp.Variables.el))
#newadcp.Variables.ua=np.hstack((adcp.Variables.ua,newadcp.Variables.ua))
#newadcp.Variables.va=np.hstack((adcp.Variables.va,newadcp.Variables.va))
newadcp=cp.copy(adcp)               
newadcp.Variables.matlabTime=modifiedtime



#Perform Harmonic Analysis of velocities and elevation
harmo1=adcp12.Utils.Harmonic_analysis(elevation=True,velocity=False)

harmo2=adcp12.Utils.Harmonic_analysis(elevation=False,velocity=True)


#Need a reconstruction with original measured data

velos1=adcps.Utils.Harmonic_reconstruction(harmo1)

velos2=adcps.Utils.Harmonic_reconstruction(harmo2)

#Reconstruct velocities
#not yet, reconstucted without gaps
#velos3=newadcp.Utils.Harmonic_reconstruction(harmo)

#velos4=newadcp.Utils.Harmonic_reconstruction(harmo1)

#newadcp.Plots.plot_xy(modifiedtime,velos3['h'])
#plt.savefig('Measured_el_data_no_gap_bad_data',format='png')
#newadcp.Plots.plot_xy(modifiedtime,velos4['u'])
#plt.savefig('Measure_uvel_data_no_gap_bad_data',format='png')
#newadcp.Plots.plot_xy(modifiedtime,velos4['v'])
#plt.savefig('Measure_vvel_data_no_gap_bad_data',format='png')

#Differences
diff_el = abs(adcps.Variables.el - velos1['h'])
diff_ua = abs(adcps.Variables.ua - velos2['u'])
diff_va = abs(adcps.Variables.va - velos2['v'])

#Save plots
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos1['h'])
#plt.savefig('Measured_el_data_w_gap__bad_data',format='png')
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos2['u'])
#plt.savefig('Measured_uvel_data_w_gap__bad_data',format='png')
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos2['v'])
#plt.savefig('Measured_vvel_data_w_gap__bad_data',format='png')

#elevation

harmo7=(adcps.Utils.Harmonic_analysis(elevation=True,velocity=False))
velos7=adcps.Utils.Harmonic_reconstruction(harmo7)

while (diff_el > (np.mean(diff_el))3.5*np.std(diff_el))).any():
	harmo7=(adcps.Utils.Harmonic_analysis(elevation=True,velocity=False))
	velos7=adcps.Utils.Harmonic_reconstruction(harmo7)
	diff_el = abs(adcps.Variables.el - velos7['h'])
	x = np.where(diff_el>(np.mean(diff_el)+3.5*np.std(diff_el)))
	diff_el=np.delete(diff_el,x)
        adcps.Variables.el = np.delete(adcps.Variables.el,x)
        adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,x)
	print shape(diff_el)
	print shape(x)
harmo7=(adcps.Utils.Harmonic_analysis(elevation=True,velocity=False))
velos7=adcps.Utils.Harmonic_reconstruction(harmo7)
diff_el=abs(adcps.Variables.el-velos7['h'])

plt.plot(adcps.Variables.matlabTime,diff_el,'.')
show()


#velos8 = newadcp.Utils.Harmonic_reconstruction(harmo8)
#plt.plot(modifedtime,velos8['h'])
#show()




#Perform a new Harmonic Analysis of velocities and elevation
#harmo3=adcp12.Utils.Harmonic_analysis(elevation=True,velocity=False)

#harmo4=adcp12.Utils.Harmonic_analysis(elevation=False,velocity=True)


#Need another reconstruction with looped through data (no bad data)

#velos5=adcps.Utils.Harmonic_reconstruction(harmo3)

#velos6=adcps.Utils.Harmonic_reconstruction(harmo4)



#ua vel
#m = 0
#x = []

#for i in diff_ua:
#        if i > 1.5*(np.mean(adcps.Variables.ua)):
#                x.append(m)
#                m+=1

#b=np.delete(diff_ua,x)
#bb=np.delete(adcps.Variables.matlabTime,x)
#plt.plot(bb,b)
#show()


#va vel
#m = 0
#x = []

#for i in diff_va:
#        if i > 1.5*(np.mean(adcps.Variables.va)):
#                x.append(m)
#                m+=1
#
#c=np.delete(diff_va,x)
#cc=np.delete(adcps.Variables.matlabTime,x)
#plt.plot(cc,c)
#show()

#Need to perform a new harmonic analysis now


#Old plots

#Plot #1: time el
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.el)
#adcps.Plots.plot_xy(modifiedtime,adcps.Variables.el)
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos['h'])
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos1['u'])
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos1['v'])
	#save in el, ua, and va

#Plot #2: times u vel
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.u)

#Plot #3:  times v vel
##adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.v)

#Plot #4:  time surf
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.surf)

#Plot #5:  time va vel
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.va)

#Plot #6:  times ua vel
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.ua)






#Project for Tom

#with open('el1.csv','w') as fp:
#	a = csv.writer(fp,delimiter=',')
#	data=np.vstack((adcps.Variables.matlabTime,adcps.Variables.el,velos['h'],(velos['h']-adcps.Variables.el)))
#	a.writerows(data.T)
	#plt.savefig(adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.el),format='png')
	#plt.savefig(adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos['h']),format='png')
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos['h']-adcps.Variables.el)

#with open('ua1.csv','w') as fq:
#	a = csv.writer(fq,delimiter=',')
#	data=np.vstack((adcps.Variables.matlabTime,adcps.Variables.ua,velos1['u'],(adcps.Variables.ua-velos1['u'])))
	#plt.savefig(adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.ua))
#	plt.savefig(adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos1['u']))
#	a.writerows(data.T)

#with open('va1.csv','w') as fy:
#	a = csv.writer(fy,delimiter=',')
#	data=np.vstack((adcps.Variables.matlabTime,adcps.Variables.va,velos1['v'],(adcps.Variables.va-velos1['v'])))
#	a.writerows(data.T)
#	plt.savefig(adcps.Plots.plot_xy(adcps.Variables.matlabTime,adcps.Variables.va))
#	plt.savefig(adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos1['v']))
	


