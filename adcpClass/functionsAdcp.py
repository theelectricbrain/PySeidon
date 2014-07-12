#!/usr/bin/python2.7
# encoding: utf-8

import numpy as np

class FunctionsAdcp:
    ''' '''
    def __init__(self,cls):
        self._var = cls.Variables
        self._debug = cls._debug

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
