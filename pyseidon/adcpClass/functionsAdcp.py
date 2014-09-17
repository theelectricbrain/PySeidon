#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np

class FunctionsAdcp:
    ''''Utils' subset of FVCOM class gathers useful functions""" '''
    def __init__(self, variable, plot, History, debug=False):
        self._var = variable
        self._plot = plot
        self._History = History
        #Create pointer to FVCOM class
        variable = self._var
        History = self._History

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
