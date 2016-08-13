

#Combining all ADCP files into one file


#Script for validating one specific ADCP observed data with one specific 2D FVCOM simulated data
import os
from pyseidon_dvt import *

#Initialize ADCP data
adcp = ADCP('/EcoII/acadia_uni/workspace/observed/DG/ADCP/DG-140115-BPb_avg10_NaNsurf_fixed.mat')

#Choose and Initialize 2D FVCOM data
fvcom2Dpath = "/EcoII/acadia_uni/workspace/simulated/FVCOM/dncoarse/calibration/bottom_roughness/2D/"
fvcom2Drun = "/output/dn_coarse_0001.nc"
brconst = raw_input("Which bottom roughness constant to test? (0.0010, 0.0015, 0.0020, 0.002486, 0.0025, 0.0030, 0.0035): ")
fvcom = FVCOM(fvcom2Dpath+brconst+fvcom2Drun)

#Change directory and validate data
os.chdir("/array/home/137463l/Validation/Results2D")
val = Validation(adcp, fvcom, flow='daf')
plot_choice = raw_input("View Plots? (y/n)")
if plot_choice == "y":
	val.validate_data(filename='val2D_'+brconst, save_csv=True, plot=True)
else:
	val.validate_data(filename='val2D_'+brconst, save_csv=True, plot=False)

