#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from utide import ut_solv, ut_reconstr
from miscellaneous import mattime_to_datetime 

class FunctionsTidegauge:
    """'Utils' subset of TideGauge class gathers useful functions"""
    def __init__(self, variable, plot, History, debug=False):
        self._var = variable
        self._plot = plot
        self._History = History
        #Create pointer to FVCOM class
        variable = self._var
        History = self._History

    def harmonics(self, time_ind=slice(None), **kwarg):
        """
        Description:
        This function performs a harmonic analysis on the sea surface elevation
        time series or the velocity components timeseries.

        Outputs:
          - harmo = harmonic coefficients, dictionary

        Keywords:
          - time_ind = time indices to work in, list of integers

        Options:
        Options are the same as for ut_solv, which are shown below with
        their default values:
            conf_int=True; cnstit='auto'; notrend=0; prefilt=[]; nodsatlint=0;
            nodsatnone=0; gwchlint=0; gwchnone=0; infer=[]; inferaprx=0;
            rmin=1; method='cauchy'; tunrdn=1; linci=0; white=0; nrlzn=200;
            lsfrqosmp=1; nodiagn=0; diagnplots=0; diagnminsnr=2;
            ordercnstit=[]; runtimedisp='yyy'

        Notes:
        For more detailed information about ut_solv, please see
        https://github.com/wesleybowman/UTide
        """
        harmo = ut_solv(self._var.matlabTime[time_ind],
                       self._var.el, [],
                       self._var.lat, **kwarg)
        return harmo

    def reconstr(self, harmo, time_ind=slice(None), **kwarg):
        """
        Description:
        This function reconstructs the velocity components or the surface elevation
        from harmonic coefficients.
        Harmonic_reconstruction calls ut_reconstr. This function assumes harmonics
        (ut_solv) has already been executed.

        Inputs:
          - Harmo = harmonic coefficient from harmo_analysis

        Output:
          - Reconstruct = reconstructed signal, dictionary

        Keywords:
          - time_ind = time indices to process, list of integers

        Options:
        Options are the same as for ut_reconstr, which are shown below with
        their default values:
            cnstit = [], minsnr = 2, minpe = 0

        Notes:
        For more detailed information about ut_reconstr, please see
        https://github.com/wesleybowman/UTide
        """
        time = self._var.matlabTime[time_ind]
        ts_recon, _ = ut_reconstr(time, harmo, **kwarg)
        return ts_recon

    def mattime2datetime(self, mattime, debug=False):
        """
        Description:
        Output the time (yyyy-mm-dd, hh:mm:ss) corresponding to
        a given matlab time

        Inputs:
          - mattime = matlab time (floats)
        """  
        time = mattime_to_datetime(mattime, debug=debug)   
        print time[0]

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
