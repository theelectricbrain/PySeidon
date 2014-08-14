#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import scipy.io as sio
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr


class Drifter:
    def __init__(self, filename):
        self.load(filename)


