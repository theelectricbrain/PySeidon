#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import seaborn
import pandas as pd

class PlotsTidegauge:
    """
    **'Plots' subset of Tidegauge class gathers plotting functions**
    """
    def __init__(self, variable, debug=False):
        self._debug = debug
        setattr(self, '_var', variable)

    def _def_fig(self):
        """Defines figure window"""
        self._fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')

    def plot_xy(self, x, y, xerror=[], yerror=[],
                title=' ', xLabel=' ', yLabel=' ', dump=False, **kwargs):
        """
        Simple X vs Y plot

        Inputs:
          - x = 1D array
          - y = 1D array

        Options:
          - xerror = error on 'x', 1D array
          - yerror = error on 'y', 1D array
          - title = plot title, string
          - xLabel = title of the x-axis, string
          - yLabel = title of the y-axis, string
          - dump = boolean, dump profile data in csv file
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        #fig = plt.figure(figsize=(18,10))
        #plt.rc('font',size='22')
        self._def_fig()
        self._ax = self._fig.add_subplot(111)         
        self._ax.plot(x, y, label=title)
        scale = 1
        self._ax.set_ylabel(yLabel)
        self._ax.set_xlabel(xLabel)
        self._ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        self._ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        self._ax.grid(b=True, which='major', color='w', linewidth=1.5)
        self._ax.grid(b=True, which='minor', color='w', linewidth=0.5)
        if not yerror==[]:
            self._ax.fill_between(x, y-yerror, y+yerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)        
        if not xerror==[]:
            self._ax.fill_betweenx(y, x-xerror, x+xerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)
        if (not xerror==[]) or (not yerror==[]): 
            blue_patch = mpatches.Patch(color='#089FFF',
                         label='Standard deviation',alpha=0.2)
            plt.legend(handles=[blue_patch],loc=1, fontsize=12)
            #plt.legend([blue_patch],loc=1, fontsize=12)

        self._fig.show()
        if dump: self._dump_profile_data_as_csv(x, y,xerror=xerror, yerror=yerror,
                                                title=title, xLabel=xLabel,
                                                yLabel=yLabel, **kwargs)     

    def _dump_profile_data_as_csv(self, x, y, xerror=[], yerror=[],
                                 title=' ', xLabel=' ', yLabel=' ', **kwargs):
        """
        Dumps profile data in csv file

        Inputs:
          - x = 1D array
          - y = 1D array

        Options:
          - xerror = error on 'x', 1D array
          - yerror = error on 'y', 1D array
          - title = file name, string
          - xLabel = name of the x-data, string
          - yLabel = name of the y-data, string
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        if title == ' ': title = 'dump_profile_data'
        filename=title + '.csv'
        if xLabel == ' ': xLabel = 'X'
        if yLabel == ' ': yLabel = 'Y'
        if not xerror == []:
            df = pd.DataFrame({xLabel:x[:], yLabel:y[:], 'error': xerror[:]})
        elif not yerror == []:
            df = pd.DataFrame({xLabel:x[:], yLabel:y[:], 'error': yerror[:]})
        else:
            df = pd.DataFrame({xLabel:x[:], yLabel:y[:]})
        df.to_csv(filename, encoding='utf-8', **kwargs)    

#TR_comments: templates
#    def whatever(self, debug=False):
#        if debug or self._debug:
#            print 'Start whatever...'
#
#        if debug or self._debug:
#            print '...Passed'
