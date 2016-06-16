#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from datetime import datetime, timedelta
from os import makedirs
from os.path import exists

# Local import
from pyseidon_dvt.utilities.interpolation_utils import *

# Custom error
from pyseidon_dvt.utilities.pyseidon_error import PyseidonError

class _load_validation:
    """
    **'Variables' subset in Validation class**

    It contains the following items: ::

                             _obs. = measurement/observational variables
      Validation.Variables._|_sim. = simulated variables
                            |_struct. = dictionnary structure for validation purposes

    """
    def __init__(self, outpath, observed, simulated, flow='sf', nn=True, debug=False, debug_plot=False):
        if debug: print "..variables.."
        self.obs = observed.Variables
        self.sim = simulated.Variables
        self._nn =  nn
        # Compatibility test
        if (observed.__module__.split('.')[-1] == 'drifterClass' and
           simulated.__module__.split('.')[-1] == 'stationClass'):
            raise PyseidonError("---Station and Drifter are incompatible objects---")

        self.struct = np.array([])
        if flow == 'daf':
            self._3D = False
        else:
            self._3D = simulated.Variables._3D
            
        try:
            # Check if times coincide
            obsMax = self.obs.matlabTime[~np.isnan(self.obs.matlabTime)].max()
            obsMin = self.obs.matlabTime[~np.isnan(self.obs.matlabTime)].min()
            simMax = self.sim.matlabTime.max()
            simMin = self.sim.matlabTime.min()
            absMin = max(obsMin, simMin)
            absMax = min(obsMax, simMax)
            A = set(np.where(self.sim.matlabTime[:] >= absMin)[0].tolist())
            B = set(np.where(self.sim.matlabTime[:] <= absMax)[0].tolist())
            C = list(A.intersection(B))
            # -Correction by J.Culina 2014-
            C = sorted(C)
            # -end-
            self._C = np.asarray(C)

            a = set(np.where(self.obs.matlabTime[:] >= absMin)[0].tolist())
            b = set(np.where(self.obs.matlabTime[:] <= absMax)[0].tolist())
            c = list(a.intersection(b))
            # -Correction by J.Culina 2014-
            c = sorted(c)
            # -end-
            self._c = np.asarray(c)

            if len(C) == 0:
                raise PyseidonError("---Time between simulation and measurement does not match up---")
        except AttributeError:
            raise PyseidonError("---Observations missing matlabTime comparison is impossible---")

        # Check which variables are available in the observations for comparison
        self._obs_vars=dir(self.obs)
        if ('lon' and 'lat') not in self._obs_vars:
            raise PyseidonError("---Observations missing lon/lat comparison is impossible---")  

        # Check what kind of simulated data it is
        if simulated.__module__ .split('.')[-1] == 'stationClass':
            self._simtype = 'station'
            # Find closest point to ADCP
            ind = closest_points([self.obs.lon], [self.obs.lat],
                                simulated.Grid.lon[:],
                                simulated.Grid.lat[:])
            try:
                nameSite = ''.join(simulated.Grid.name[ind,:][0,:])
            except IndexError:  # TR: quick fix
                nameSite = ''.join(simulated.Grid.name[ind])
            print "Station site: " + nameSite
            self.sim.lat = simulated.Grid.lat[ind]
            el = self.sim.el[:, ind].ravel()
            ua = self.sim.ua[:, ind].ravel()
            va = self.sim.va[:, ind].ravel()
            if self._3D:
                u = np.squeeze(self.sim.u[:, :,ind])
                v = np.squeeze(self.sim.v[:, :,ind])
                sig = np.squeeze(simulated.Grid.siglay[:, ind])
                h = simulated.Grid.h[:]

        # Alternative simulation type
        elif simulated.__module__.split('.')[-1] == 'fvcomClass':
            self._simtype = 'fvcom'
            # Different treatment measurements come from drifter
            if not observed.__module__.split('.')[-1] == 'drifterClass':
                if debug: print "...Interpolation at measurement location..."
                el = simulated.Util2D.interpolation_at_point(self.sim.el,
                                                             self.obs.lon, self.obs.lat, nn=self._nn)
                ua = simulated.Util2D.interpolation_at_point(self.sim.ua,
                                                             self.obs.lon, self.obs.lat, nn=self._nn)
                va = simulated.Util2D.interpolation_at_point(self.sim.va,
                                                             self.obs.lon, self.obs.lat, nn=self._nn)
                if self._3D:
                    u = simulated.Util3D.interpolation_at_point(self.sim.u,
                                                                self.obs.lon, self.obs.lat, nn=self._nn)
                    v = simulated.Util3D.interpolation_at_point(self.sim.v,
                                                                self.obs.lon, self.obs.lat, nn=self._nn)
                    sig = simulated.Util3D.interpolation_at_point(simulated.Grid.siglay,
                                                                  self.obs.lon, self.obs.lat, nn=self._nn)
                    h = simulated.Util3D.interpolation_at_point(simulated.Grid.h[:],
                                                                self.obs.lon, self.obs.lat, nn=self._nn)
            else:  # Interpolation for drifter
                if debug: print "...Interpolation at measurement locations & times..."
                if self._3D:
                    lock=True
                    userInp = flow
                    while lock:
                        if userInp == 'daf':
                            if debug:
                                print 'flow comparison is depth-averaged'
                            # TR: temporary fix for proxy access
                            if self.sim._opendap:
                                uSim = np.zeros((self._C.shape[0], self.sim.ua.shape[1]))
                                vSim = np.zeros((self._C.shape[0], self.sim.va.shape[1]))
                                for i, j in enumerate(self._C):
                                    uSim[i,:] = self.sim.ua[j,:]
                                    vSim[i,:] = self.sim.va[j,:]
                            else:
                                uSim = np.squeeze(self.sim.ua[self._C,:])
                                vSim = np.squeeze(self.sim.va[self._C,:])
                            self._3D = False
                            lock = False
                        elif userInp == 'sf':
                            # Import only the surface velocities
                            # TR_comment: is surface vertical indice -1 or 0?
                            # KC : 0 is definitely the surface...
                            # TR: temporary fix for proxy access
                            if self.sim._opendap:
                                uSim = np.zeros((self._C.shape[0], self.sim.u.shape[2]))
                                vSim = np.zeros((self._C.shape[0], self.sim.v.shape[2]))
                                for i, j in enumerate(self._C):
                                    uSim[i,:] = self.sim.u[j,0,:]
                                    vSim[i,:] = self.sim.v[j,0,:]
                            else:
                                uSim = np.squeeze(self.sim.u[self._C,0,:])
                                vSim = np.squeeze(self.sim.v[self._C,0,:])
                            self._3D = False
                            lock = False
                            if debug:
                                print 'flow comparison at surface'
                        elif type(userInp) == float:
                            if debug:
                                print 'flow comparison at depth level ', float
                            # if userInp > 0.0: userInp = userInp*-1.0
                            uInterp = simulated.Util3D.interp_at_depth(self.sim.u[:], userInp, debug=debug)
                            vInterp = simulated.Util3D.interp_at_depth(self.sim.v[:], userInp, debug=debug)
                            # TR: temporary fix for proxy access
                            if self.sim._opendap:
                                uSim = np.zeros((self._C.shape[0], uInterp.shape[1]))
                                vSim = np.zeros((self._C.shape[0], vInterp.shape[1]))
                                for i, j in enumerate(self._C):
                                    uSim[i,:] = uInterp[j,:]
                                    vSim[i,:] = vInterp[j,:]
                            else:
                                uSim = np.squeeze(uInterp[self._C,:])
                                vSim = np.squeeze(vInterp[self._C,:])
                            self._3D = False
                            lock = False
                        else:
                            userInp = input("compare flow by 'daf', 'sf' or a float number only!!!")
                else:
                    # TR: temporary fix for proxy access
                    if self.sim._opendap:
                        uSim = np.zeros((self._C.shape[0], self.sim.ua.shape[1]))
                        vSim = np.zeros((self._C.shape[0], self.sim.va.shape[1]))
                        for i, j in enumerate(self._C):
                            uSim[i,:] = self.sim.ua[j,:]
                            vSim[i,:] = self.sim.va[j,:]
                    else:
                        uSim = np.squeeze(self.sim.ua[self._C,:])
                        vSim = np.squeeze(self.sim.va[self._C,:])

                # Finding the closest Drifter time to simulated data assuming measurement
                # time step way faster than model one
                indClosest = []
                for i in self._C:
                    ind = np.abs(self.obs.matlabTime[:]-self.sim.matlabTime[i]).argmin()
                    indClosest.append(ind)
                # Keep only unique values to avoid sampling in measurement gaps
                # TR: this doesn't not guarantee that the unique value kept is indeed
                #    the closest one among the values relative to the same indice!!!
                uniqCloInd, uniqInd = np.unique(indClosest, return_index=True)

                # KC: Convert of matlabTime to datetime, smooth drifter temporally!
                self.datetimes = [datetime.fromordinal(int(x))+timedelta(days=x%1) -
                                  timedelta(days=366) for x in self.obs.matlabTime[self._c]]

                uObs = self.obs.u[uniqCloInd]
                vObs = self.obs.v[uniqCloInd]
                uSim = np.squeeze(uSim[uniqInd[:],:])
                vSim = np.squeeze(vSim[uniqInd[:],:])

                # print 'uObs: \n', uObs, '\nvObs: \n', vObs
                # print 'uSim: \n', vSim.shape, '\nuSim: \n', vSim.shape

                # Interpolation of timeseries at drifter's trajectory points
                uSimInterp = np.zeros(len(uniqCloInd))
                vSimInterp = np.zeros(len(uniqCloInd))
                for i in range(len(uniqCloInd)):
                    uSimInterp[i]=simulated.Util2D.interpolation_at_point(uSim[i,:],
                                                self.obs.lon[uniqCloInd[i]],
                                                self.obs.lat[uniqCloInd[i]])
                    vSimInterp[i]=simulated.Util2D.interpolation_at_point(vSim[i,:],
                                                self.obs.lon[uniqCloInd[i]],
                                                self.obs.lat[uniqCloInd[i]])
                # print 'vSimInterp: \n', vSimInterp, '\nuSimInterp: \n', uSimInterp

        else:
            raise PyseidonError("-This type of simulations is not supported yet-")

        # Store in dict structure for compatibility purposes (except for drifters)
        if not observed.__module__.split('.')[-1] == 'drifterClass':
            if not self._3D:
                sim_mod = {'ua': ua[C], 'va': va[C], 'el': el[C]}
            else:
                sim_mod = {'ua': ua[C], 'va': va[C], 'el': el[C], 'u': u[C,:],
                           'v': v[C,:], 'siglay': sig[:], 'h': h}

            # Check what kind of observed data it is
            if observed.__module__.split('.')[-1] == 'adcpClass' or observed.__module__ == 'adcpClass':
                self._obstype = 'adcp'
                obstype = 'ADCP'
            # Alternative measurement type
            elif observed.__module__.split('.')[-1] == 'tidegaugeClass':
                self._obstype = 'tidegauge'
                obstype = 'TideGauge'
            else:
                raise PyseidonError("---This type of measurements is not supported yet---")

            # This had to be split into two parts so that a time subset could be used.
            # Specify list of data variables from observations that should be used.    
            dictlist=['ua','va','u','v','el'] 
            self._commonlist_data = [var for var in self._obs_vars if var in dictlist]
            if debug:
                print 'Data variables being used'
                print self._commonlist_data         
            obs_mod={}    
            for key in self._commonlist_data:
                obs_mod[key] = getattr(self.obs,key)
                obs_mod[key] = obs_mod[key][c,]
                
            # Specify list of nondata variables that should be used.    
            dictlist = ['bins','data']
            self._commonlist_nondata = [var for var in self._obs_vars if var in dictlist]
            if debug:
                print 'Non data variables being used'
                print self._commonlist_nondata    
            for key in self._commonlist_nondata:
                obs_mod[key] = getattr(self.obs,key)
                
        else:
            self._obstype = 'drifter'
            obstype = 'Drifter'

        #Store in dict structure for compatibility purposes
        #Common block for 'struct'
        if not observed.__module__.split('.')[-1] == 'drifterClass':
            self.struct = {'name': observed.History[0].split(' ')[-1],
                           'type': obstype,
                           'obs_lat': self.obs.lat,
                           'obs_lon': self.obs.lon,
                           'mod_lat': self.obs.lat,
                           'mod_lon': self.obs.lon,
                           'obs_timeseries': obs_mod,
                           'mod_timeseries': sim_mod,
                           'obs_time': self.obs.matlabTime[c],
                           'mod_time': self.sim.matlabTime[C],
                           '_commonlist_data': self._commonlist_data,
                           '_commonlist_nondata': self._commonlist_nondata}
        else: # Drifter's case
            self.struct = {'name': observed.History[0].split(' ')[-1],
                           'type': obstype,
                           'lat': self.obs.lat[uniqCloInd],
                           'lon': self.obs.lon[uniqCloInd],
                           'obs_timeseries': {'u': uObs, 'v': vObs},
                           'mod_timeseries': {'u': uSimInterp, 'v': vSimInterp},
                           'obs_time': self.obs.matlabTime[uniqCloInd],
                           'mod_time': self.sim.matlabTime[C],
                           '_commonlist_data': ['u', 'v']}

        if debug: print "..done"
        
        # find save_path
        # self._save_path = ''
        # for s in observed.History:
        #     if 'Create save_path ' not in s:
        #         continue
        #     else:
        #         self._save_path=s.split()[2]
        #
        # if '' is self._save_path:
        #     name = self.struct['name']
        #     self._save_path = outpath+name.split('/')[-1].split('.')[0]+'/'
        #     while exists(self._save_path):
        #         self._save_path = self._save_path[:-1] + '_bis/'
        #     makedirs(self._save_path)
        #     observed.History.append('Create save_path {}'.format(self._save_path))
        name = self.struct['name']
        self._save_path = outpath+name.split('/')[-1].split('.')[0]+'/'
        while exists(self._save_path):
            self._save_path = self._save_path[:-1] + '_bis/'
        makedirs(self._save_path)

        return
