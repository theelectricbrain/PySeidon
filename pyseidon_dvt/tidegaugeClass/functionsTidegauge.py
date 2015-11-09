#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

from utide import solve, reconstruct
from pyseidon_dvt.utilities.miscellaneous import mattime_to_datetime

class FunctionsTidegauge:
    """
    **'Utils' subset of TideGauge class gathers useful functions**
    """
    def __init__(self, variable, plot, History, debug=False):
        self._plot = plot
        #Create pointer to FVCOM class
        setattr(self, '_var', variable)
        setattr(self, '_History', History)

    def harmonics(self, time_ind=slice(None), **kwarg):
        """
        This function performs a harmonic analysis on the sea surface elevation
        time series or the velocity components timeseries.

        Outputs:
          - harmo = harmonic coefficients, dictionary

        Options:
          - time_ind = time indices to work in, list of integers

        Utide's options:
        Options are the same as for 'solve', which are shown below with
        their default values:
            conf_int=True; cnstit='auto'; notrend=0; prefilt=[]; nodsatlint=0;
            nodsatnone=0; gwchlint=0; gwchnone=0; infer=[]; inferaprx=0;
            rmin=1; method='cauchy'; tunrdn=1; linci=0; white=0; nrlzn=200;
            lsfrqosmp=1; nodiagn=0; diagnplots=0; diagnminsnr=2;
            ordercnstit=[]; runtimedisp='yyy'

        *Notes*
        For more detailed information about 'solve', please see
        https://github.com/wesleybowman/UTide
        """
        harmo = solve(self._var.matlabTime[time_ind],
                      self._var.el, None,
                      self._var.lat, **kwarg)
        return harmo

    def reconstr(self, harmo, time_ind=slice(None), **kwarg):
        """
        This function reconstructs the velocity components or the surface elevation
        from harmonic coefficients.
        Harmonic_reconstruction calls 'reconstruct'. This function assumes harmonics
        ('solve') has already been executed.

        Inputs:
          - Harmo = harmonic coefficient from harmo_analysis

        Output:
          - Reconstruct = reconstructed signal, dictionary

        Keywords:
          - time_ind = time indices to process, list of integers

        Utide's options:
        Options are the same as for 'reconstruct', which are shown below with
        their default values:
            cnstit = [], minsnr = 2, minpe = 0

        *Notes*
        For more detailed information about 'reconstruct', please see
        https://github.com/wesleybowman/UTide
        """
        time = self._var.matlabTime[time_ind]
        ts_recon = reconstruct(time, harmo, **kwarg)
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
        return time

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
