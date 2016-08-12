# Rough Work


























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


#Save plots
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos1['h'])
#plt.savefig('Measured_el_data_w_gap__bad_data',format='png')
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos2['u'])
#plt.savefig('Measured_uvel_data_w_gap__bad_data',format='png')
#adcps.Plots.plot_xy(adcps.Variables.matlabTime,velos2['v'])
#plt.savefig('Measured_vvel_data_w_gap__bad_data',format='png')


#Perform a new Harmonic Analysis of velocities and elevation
#harmo3=adcp12.Utils.Harmonic_analysis(elevation=True,velocity=False)

#harmo4=adcp12.Utils.Harmonic_analysis(elevation=False,velocity=True)


#Need another reconstruction with looped through data (no bad data)

#velos5=adcps.Utils.Harmonic_reconstruction(harmo3)

#velos6=adcps.Utils.Harmonic_reconstruction(harmo4)

