#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import sys
from datetime import datetime, timedelta

#Local import
from pyseidon.utilities.interpolation_utils import *

# Custom error
from pyseidon.utilities.pyseidon_error import PyseidonError

class _load_validation:
    """
    **'Variables' subset in Validation class**

    It contains the following items: ::

                             _obs. = measurement/observational variables
      Validation.Variables._|_sim. = simulated variables
                            |_struct. = dictionnary structure for validation purposes
    """
    def __init__(self, observed, simulated, flow='sf', debug=False, debug_plot=False):
        if debug: print "..variables.."
        self.obs = observed.Variables
        self.sim = simulated.Variables
        self.struct = np.array([])
        self._3D = simulated.Variables._3D

        #Check if times coincide
        obsMax = self.obs.matlabTime[~np.isnan(self.obs.matlabTime)].max()
        obsMin = self.obs.matlabTime[~np.isnan(self.obs.matlabTime)].min()
        simMax = self.sim.matlabTime.max()
        simMin = self.sim.matlabTime.min()
        absMin = max(obsMin, simMin)
        absMax = min(obsMax, simMax)
        A = set(np.where(self.sim.matlabTime[:] >= absMin)[0].tolist())
        B = set(np.where(self.sim.matlabTime[:] <= absMax)[0].tolist())
        C = list(A.intersection(B))
        #-Correction by J.Culina 2014-
        C = sorted(C)
        #-end-
        self._C = C

        a = set(np.where(self.obs.matlabTime[:] >= absMin)[0].tolist())
        b = set(np.where(self.obs.matlabTime[:] <= absMax)[0].tolist())
        c = list(a.intersection(b))
        #-Correction by J.Culina 2014-
        c = sorted(c)
        #-end-
        self._c = c

        if len(C) == 0:
            raise PyseidonError("---Time between simulation and measurement does not match up---")

        #Check what kind of simulated data it is
        if simulated.__module__=='pyseidon.stationClass.stationClass':
            self._simtype = 'station'
            #Find closest point to ADCP
            ind = closest_points([self.obs.lon], [self.obs.lat],
                                simulated.Grid.lon[:],
                                simulated.Grid.lat[:])
            nameSite = ''.join(simulated.Grid.name[ind,:][0,:])
            print "Station site: " + nameSite
            self.sim.lat = simulated.Grid.lat[ind]
            el = self.sim.el[:, ind].flatten()
            ua = self.sim.ua[:, ind].flatten()
            va = self.sim.va[:, ind].flatten()
            if self._3D:
                u = np.squeeze(self.sim.u[:, :,ind])
                v = np.squeeze(self.sim.v[:, :,ind])
                sig = np.squeeze(simulated.Grid.siglay[:, ind])

        #Alternative simulation type
        elif simulated.__module__=='pyseidon.fvcomClass.fvcomClass':
            self._simtype = 'fvcom'
            #Different treatment measurements come from drifter
            if not observed.__module__=='pyseidon.drifterClass.drifterClass':
                if debug: print "...Interpolation at measurement location..."
                el=simulated.Util2D.interpolation_at_point(self.sim.el,
                                                           self.obs.lon, self.obs.lat)
                ua=simulated.Util2D.interpolation_at_point(self.sim.ua,
                                                           self.obs.lon, self.obs.lat)
                va=simulated.Util2D.interpolation_at_point(self.sim.va,
                                                           self.obs.lon, self.obs.lat)
                if self._3D:
                   u=simulated.Util3D.interpolation_at_point(self.sim.u,
                                                             self.obs.lon, self.obs.lat)
                   v=simulated.Util3D.interpolation_at_point(self.sim.v,
                                                             self.obs.lon, self.obs.lat)
                   sig=simulated.Util3D.interpolation_at_point(simulated.Grid.siglay,
                                                               self.obs.lon, self.obs.lat)
            else: #Interpolation for drifter
                if debug: print "...Interpolation at measurement locations & times..."
                if self._3D:
                    lock=True
                    while lock:
                        userInp = flow
                        if userInp == 'daf':
                            if debug:
                                print 'flow comparison is depth-averaged'
                            uSim = np.squeeze(self.sim.ua[self._C,:])
                            vSim = np.squeeze(self.sim.va[self._C,:])
                            self._3D = False
                            lock=False
                        elif userInp == 'sf':
                            #Import only the surface velocities
                            #TR_comment: is surface vertical indice -1 or 0?
                            #KC : 0 is definitely the surface...
                            uSim = np.squeeze(self.sim.u[self._C,0,:])
                            vSim = np.squeeze(self.sim.v[self._C,0,:])
                            self._3D = False
                            lock=False
                            if debug:
                                print 'flow comarison at surface'
                        elif type(userInp) == float:
                            if debug:
                                print 'flow comparison at depth level ', float
                            if userInp > 0.0: userInp = userInp*-1.0
                            uInterp = simulated.Util3D.interp_at_depth(self.sim.u[:], userInp, debug=debug)
                            uSim = np.squeeze(uInterp[self._C,:])
                            vInterp = simulated.Util3D.interp_at_depth(self.sim.v[:], userInp, debug=debug)
                            vSim = np.squeeze(vInterp[self._C,:])
                            self._3D = False
                            lock=False
                        else:
                            print "compare flow by 'daf', 'sf' or a float number only!!!"
                else:
                    uSim = np.squeeze(self.sim.ua[self._C,:])
                    vSim = np.squeeze(self.sim.va[self._C,:])

                #Finding the closest Drifter time to simulated data assuming measurement
                #time step way faster than model one
                indClosest = []
                for i in self._C:
                    ind = np.abs(self.obs.matlabTime[:]-self.sim.matlabTime[i]).argmin()
                    indClosest.append(ind)
                #Keep only unique values to avoid sampling in measutement gaps
                #TR: this doesn't not guarantee that the unique value kept is indeed
                #    the closest one among the values relative to the same indice!!!
                uniqCloInd, uniqInd = np.unique(indClosest, return_index=True)

                #KC: Convert of matlabTime to datetime, smooth drifter temporally!
                self.datetimes = [datetime.fromordinal(int(x))+timedelta(days=x%1)\
                   - timedelta(days = 366) for x in self.obs.matlabTime[self._c]]

                #uObs, vObs, t_s, dt_start = smooth(self.obs.u[self._c], \
                #        self.datetimes, self.obs.v[self._c], \
                #        self.datetimes, delta_t=1, debug=True)

                uObs = self.obs.u[uniqCloInd]
                vObs = self.obs.v[uniqCloInd]
                uSim = np.squeeze(uSim[uniqInd[:],:])
                vSim = np.squeeze(vSim[uniqInd[:],:])

                #print 'uObs: \n', uObs, '\nvObs: \n', vObs
                #print 'uSim: \n', vSim.shape, '\nuSim: \n', vSim.shape

                #Interpolation of timeseries at drifter's trajectory points
                uSimInterp = np.zeros(len(uniqCloInd))
                vSimInterp = np.zeros(len(uniqCloInd))
                for i in range(len(uniqCloInd)):
                    uSimInterp[i]=simulated.Util2D.interpolation_at_point(uSim[i,:],
                                                self.obs.lon[uniqCloInd[i]],
                                                self.obs.lat[uniqCloInd[i]])
                    vSimInterp[i]=simulated.Util2D.interpolation_at_point(vSim[i,:],
                                                self.obs.lon[uniqCloInd[i]],
                                                self.obs.lat[uniqCloInd[i]])
                #print 'vSimInterp: \n', vSimInterp, '\nuSimInterp: \n', uSimInterp

        else:
            raise PyseidonError("-This type of simulations is not supported yet-")

        #Store in dict structure for compatibility purposes (except for drifters)
        if not observed.__module__=='pyseidon.drifterClass.drifterClass':
            if not self._3D:
                sim_mod={'ua':ua[C],'va':va[C],'elev':el[C]}
            else:
                sim_mod={'ua':ua[C],'va':va[C],'elev':el[C],'u':u[C,:],
                         'v':v[C,:],'siglay':sig[:]}


            #Check what kind of observed data it is
            if observed.__module__=='pyseidon.adcpClass.adcpClass' or observed.__module__ == 'adcpClass':
                self._obstype = 'adcp'
                obstype='ADCP'
                obs_mod={'ua':self.obs.ua[c],'va':self.obs.va[c],'elev':self.obs.surf[c],
                         'u':self.obs.east_vel[c,:],'v':self.obs.north_vel[c,:],
                         'bins':self.obs.bins[:]}

            #Alternative measurement type
            elif observed.__module__=='pyseidon.tidegaugeClass.tidegaugeClass':
                self._obstype = 'tidegauge'
                obstype='TideGauge'
                obs_mod = {'data':self.obs.RBR.data, 'elev':self.obs.el[c]}

            else:
                raise PyseidonError("-This type of measurements is not supported yet-")
        else:
            self._obstype = 'drifter'
            obstype='Drifter'

        #Store in dict structure for compatibility purposes
        #Common block for 'struct'
        if not observed.__module__=='pyseidon.drifterClass.drifterClass':
            self.struct = {'name': observed.History[0].split(' ')[-1],
                           'type':obstype,
                           'lat':self.obs.lat,
                           'lon':self.obs.lon,
                           'obs_timeseries':obs_mod,
                           'mod_timeseries':sim_mod,
                           'obs_time':self.obs.matlabTime[c],
                           'mod_time':self.sim.matlabTime[C]}
        else:#Drifter's case
            self.struct = {'name': observed.History[0].split(' ')[-1],
                           'type':obstype,
                           'lat':self.obs.lat[uniqCloInd],
                           'lon':self.obs.lon[uniqCloInd],
                           'obs_timeseries':{'u': uObs, 'v': vObs},
                           'mod_timeseries':{'u': uSimInterp, 'v': vSimInterp},
                           'obs_time':self.obs.matlabTime[uniqCloInd],
                           'mod_time':self.sim.matlabTime[C]}

        if debug: print "..done"
