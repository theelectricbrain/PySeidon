#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division

import os.path as os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf

def plotRegression(tidalStatClass, lr, savepath='', fname='', debug=False):
    """
    Plots a visualization of the output from linear regression,
    including confidence intervals for predictands and slope.

    If a savepath and filename is defined, exports the plot as an image file
    to that location. Filenames should include the image file name and extension.

    Returns:
      fig, ax: maplotlib objects
    """
    if debug : print "Plotting linear regression"

    #define figure frame
    fig = plt.figure(figsize=(18,10))
    plt.rc('font',size='22')
    ax = fig.add_subplot(111)

    ax.scatter(tidalStatClass.model, tidalStatClass.observed, c='b', marker='+', alpha=0.5)

    ## plot regression line
    mod_max = np.amax(tidalStatClass.model)
    mod_min = np.amin(tidalStatClass.model)
    upper_intercept = lr['intercept'] + lr['pred_CI_width']
    lower_intercept = lr['intercept'] - lr['pred_CI_width']
    ax.plot([mod_min, mod_max], [mod_min * lr['slope'] + lr['intercept'],
            mod_max * lr['slope'] + lr['intercept']],
            color='k', linestyle='-', linewidth=2, label='Linear fit')

    ## plot CI's for slope
    ax.plot([mod_min, mod_max], [mod_min * lr['slope_CI'][0] + lr['intercept_CI'][0],
                                 mod_max * lr['slope_CI'][0] + lr['intercept_CI'][0]],
             color='r', linestyle='--', linewidth=2)
    ax.plot([mod_min, mod_max], [mod_min * lr['slope_CI'][1] + lr['intercept_CI'][1],
                                 mod_max * lr['slope_CI'][1] + lr['intercept_CI'][1]],
             color='r', linestyle='--', linewidth=2, label='Slope CI')

    ## plot CI's for predictands
    ax.plot([mod_min, mod_max], [mod_min * lr['slope'] + upper_intercept,
                                 mod_max * lr['slope'] + upper_intercept],
             color='g', linestyle='--', linewidth=2)
    ax.plot([mod_min, mod_max], [mod_min * lr['slope'] + lower_intercept,
                                 mod_max * lr['slope'] + lower_intercept],
             color='g', linestyle='--', linewidth=2, label='Predictand CI')

    ax.set_xlabel('Modeled Data')
    ax.set_ylabel('Observed Data')
    fig.suptitle('Modeled vs. Observed {}: Linear Fit'.format(tidalStatClass.kind))
    plt.legend(loc='lower right', shadow=True)

    r_string = 'R Squared: {}'.format(np.around(lr['r_2'], decimals=3))
    plt.title(r_string)

    # Pretty plot
    # df = DataFrame(data={'model': tidalStatClass.model.ravel(),
    #                      'observed':tidalStatClass.observed.ravel()})
    # seaborn.set(style="darkgrid")
    # color = seaborn.color_palette()[2]
    # g = seaborn.jointplot("model", "observed", data=df, kind="reg",
    #                       xlim=(df.model.min(), df.model.max()),
    #                       ylim=(df.observed.min(), df.observed.max()),
    #                       color=color, size=7)
    # plt.suptitle('Modeled vs. Observed {}: Linear Fit'.format(tidalStatClass.kind))
    # KC: Changed save parameter to be a savepath - making a huge assumption here
    # that people are able to enter in the savepath correctly / exists etc.
    if savepath.strip() and fname.strip():
        if os.exists(savepath):
            fig.savefig(savepath+fname)
            fig.clear()
            plt.close(fig)
    else:
        fig.show()
        plt.show()

    return fig, ax

def plotData(tidalStatClass, graph='time', savepath='', fname='', debug=False):
    """
    Provides a visualization of the data.

    Takes an option which determines the kind of graph to be made.
    time: plots the model data against the observed data over time
    scatter : plots the model data vs. observed data

    If a savepath and filename is defined, exports the plot as an image file
    to that location. Filenames should include the image file name and extension.
    """
    if debug: print "Plotting time-series..."
    #define figure frame
    fig = plt.figure(figsize=(18,10))
    plt.rc('font',size='22')
    ax = fig.add_subplot(111)

    if (graph == 'time'):
        ax.plot(tidalStatClass.times, tidalStatClass.model, label='Model Predictions')
        ax.plot(tidalStatClass.times, tidalStatClass.observed, color='r',
                 label='Observed Data')
        ax.set_xlabel('Time')
        if tidalStatClass.kind == 'elevation':
            ax.set_ylabel('Elevation (m)')
        if tidalStatClass.kind == 'speed':
            ax.set_ylabel('Flow speed (m/s)')
        if tidalStatClass.kind == 'direction':
            ax.set_ylabel('Flow direction (deg.)')
        if tidalStatClass.kind == 'u velocity':
            ax.set_ylabel('U velocity (m/s)')
        if tidalStatClass.kind == 'v velocity':
            ax.set_ylabel('V velocity (m/s)')
        if tidalStatClass.kind == 'velocity':
            ax.set_ylabel('Signed flow speed (m/s)')
        if tidalStatClass.kind == 'cubic speed':
            ax.set_ylabel('Cubic speed (m3/s3)')

        fig.suptitle('Predicted and Observed {}'.format(tidalStatClass.kind))
        ax.legend(shadow=True)

    if (graph == 'scatter'):
        ax.scatter(tidalStatClass.model, tidalStatClass.observed, c='b', alpha=0.5)
        ax.set_xlabel('Predicted Height')
        ax.set_ylabel('Observed Height')
        fig.suptitle('Predicted vs. Observed {}'.format(tidalStatClass.kind))

    if savepath.strip() and fname.strip():
        if os.exists(savepath):
            fig.savefig(savepath+fname)
            fig.clear()
            plt.close(fig)
    else:
        fig.show()
        plt.show()

def taylorDiagram(benchmarks, savepath='', fname='', labels=True, debug=False):
    """
    References:
      Taylor, K.E.:  Summarizing multiple aspects of model performance in a single diagram.
      J. Geophys. Res., 106, 7183-7192, 2001
      (also see PCMDI Report 55, http://www-pcmdi.llnl.gov/publications/ab55.html)

      IPCC, 2001: Climate Change 2001: The Scientific Basis,
      Contribution of Working Group I to the Third Assessment Report of the Intergovernmental Panel on Climate Change
      Houghton, J.T., Y. Ding, D.J. Griggs, M. Noguer, P.J. van der Linden, X. Dai, K. Maskell, and C.A. Johnson (eds.)
      Cambridge University Press, Cambridge, United Kingdom and New York, NY, USA, 881 pp.
      (see http://www.grida.no/climate/ipcc_tar/wg1/317.htm#fig84)

    Code inspired by Yannick Copin's code (see https://github.com/ycopin)
    """
    if debug: print "Plotting time-series..."

    # Setting up graph
    tr = PolarAxes.PolarTransform()
    # Correlation labels
    rlocs = np.concatenate((np.arange(10)/10.,[0.95,0.99]))
    tlocs = np.arccos(rlocs)        # Conversion to polar angles
    gl1 = gf.FixedLocator(tlocs)    # Positions
    tf1 = gf.DictFormatter(dict(zip(tlocs, map(str,rlocs))))
    # Standard deviation axis extent
    smin = 0
    smax = 1.5
    ghelper = fa.GridHelperCurveLinear(tr, extremes=(0,np.pi/2, smin, smax), grid_locator1=gl1, tick_formatter1=tf1)
    fig = plt.figure(figsize=(18,10))
    rect=111
    ax = fa.FloatingSubplot(fig, rect, grid_helper=ghelper)
    fig.add_subplot(ax)
    # Adjust axes
    ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("Correlation")
    ax.axis["left"].set_axis_direction("bottom") # "X axis"
    ax.axis["left"].label.set_text("Standard deviation")
    ax.axis["right"].set_axis_direction("top")   # "Y axis"
    ax.axis["right"].toggle(ticklabels=True)
    ax.axis["right"].major_ticklabels.set_axis_direction("left")
    ax.axis["bottom"].set_visible(False)         # Useless
    # Contours along standard deviations
    ax.grid(False)
    _ax = ax                   # Graphical axes
    ax = ax.get_aux_axes(tr)   # Polar coordinates
    # Reference lines
    x95 = [0.05, 13.9] # For Prcp, this is for 95th level (r = 0.195)
    y95 = [0.0, 71.0]
    x99 = [0.05, 19.0] # For Prcp, this is for 99th level (r = 0.254)
    y99 = [0.0, 70.0]
    ax.plot(x95,y95,color='k')
    ax.plot(x99,y99,color='k')
    # Add reference point and stddev contour
    l, = ax.plot(0.0, 1.0, 'k*', ls='', ms=10, label='Reference')
    t = np.linspace(0, np.pi/2)
    r = np.zeros_like(t) + 1.0
    ax.plot(t,r, 'k--', label='_')
    samplePoints = [l]

    # Plot points
    sampleLenght = benchmarks['Type'].shape[0]
    colors = plt.matplotlib.cm.jet(np.linspace(0,1,sampleLenght))
    for i in range(sampleLenght):
        if labels:
            l, = ax.plot(benchmarks['NRMSE'][i]/100.0, benchmarks['r2'][i],
                        marker='$%d$' % (i+1), ms=10, ls='', mfc=colors[i], mec=colors[i],
                        label= benchmarks['Type'][i] + " " + benchmarks['gear'][i])
        else:
            l, = ax.plot(benchmarks['NRMSE'][i]/100.0, benchmarks['r2'][i],
                        marker='$%d$' % (i+1), ms=10, ls='', mfc=colors[i], mec=colors[i])
        samplePoints.append(l)
    t = np.linspace(0, np.pi/2)
    r = np.zeros_like(t) + 1.0
    ax.plot(t,r, 'k--', label='_')

    # Add NRMS contours, and label them
    rs, ts = np.meshgrid(np.linspace(smin, smax), np.linspace(0,np.pi/2))
    # Compute centered RMS difference
    rms = np.sqrt(1.0 + rs**2 - 2*rs*np.cos(ts))
    contours = ax.contour(ts, rs, rms, 5, colors='0.5')
    ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend
    if labels:
        lgd = fig.legend(samplePoints,[p.get_label() for p in samplePoints],
                         numpoints=1, prop=dict(size='small'), loc='upper right')
    fig.tight_layout()

    if savepath.strip() and fname.strip():
        if os.exists(savepath):
            if labels:
                fig.savefig(savepath+fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
            else:
                fig.savefig(savepath+fname, bbox_inches='tight')
            fig.clear()
            plt.close(fig)
    else:
        fig.show()
        plt.show()

    return fig, ax

def benchmarksMap(benchmarks, adcps, fvcom, savepath='', fname='', debug=False):
    """
    Plots bathymetric map & model validation benchmarks

    Inputs:
      - benchmarks = benchmark attribute from Validation class
      - adcps = list or tuple of ADCP objects
      - fvcom = FVCOM object
    Options:
      - savepath = folder path for saving plot, string
      - fname = filename for saving plot, string
    """
    #if debug: print "Computing flow speed"
    #fvcom.Util2D.hori_velo_norm()
    #speed = np.mean(fvcom.Variables.hori_velo_norm[:],0)
    #if debug: print '...passed'

    # collecting names and locations of adcps
    adcpLoc={}
    try:
        for adcp in adcps:
            try:
                key = adcp.History[0].split(' ')[-1].split('/')[-1].split('.')[0]
                val = benchmarks.loc[key]
                indCS = np.where(val['Type'].values == 'cubic_speed')[0][0]
                adcpLoc[key] = {'location': [adcp.Variables.lon, adcp.Variables.lat],
                                'r2': val['r2'][indCS],
                                'NRMSE': val['NRMSE'][indCS],
                                'bias': val['bias'][indCS]}
            except KeyError:  # in case something different than adcp in the list
                pass
    except TypeError:
        key = adcps.History[0].split(' ')[-1].split('/')[-1].split('.')[0]
        adcpLoc[key] = [adcps.Variables.lon, adcps.Variables.lat]
        val = benchmarks.loc[key]
        indCS = np.where(val['Type'].values == 'cubic_speed')[0][0]
        adcpLoc[key] = {'location': [adcps.Variables.lon, adcps.Variables.lat],
                        'r2': val['r2'][indCS],
                        'NRMSE': val['NRMSE'][indCS],
                        'bias': val['bias'][indCS]}

    # Plot size and color function of R2 and RMSE
    # #background
    #cmap=plt.cm.jet
    #fvcom.Plots.colormap_var(speed, title='Averaged flow speed', mesh=False, cmap=cmap)
    fvcom.Plots.colormap_var(fvcom.Grid.h, title='Bathymetric Map & Model Validation Benchmarks', mesh=False)

    for key in adcpLoc.keys():
        print '...plotting ' + key + '...'
        r2 = adcpLoc[key]['r2']  # r2 for cubic velocity
        nrmse = adcpLoc[key]['NRMSE']  # nrmse for cubic speed
        mk = adcpLoc[key]['bias']  # over or under estimated
        if np.sign(mk) == -1.0 :
            mk = '_'
        else:
            mk = '+'
        fvcom.Plots._ax.scatter(adcpLoc[key]['location'][0],adcpLoc[key]['location'][1],
                                marker=mk, lw=2, s=100, color='red')
        fvcom.Plots._ax.annotate('r2: '+str(round(r2,2))+' |',
                                 xy=(adcpLoc[key]['location'][0],adcpLoc[key]['location'][1]),
                                 xycoords='data', xytext=(-55, -15),
                                 textcoords='offset points', ha='left',
                                 color='white', fontsize=12)
        fvcom.Plots._ax.annotate('nrmse: '+str(round(nrmse,2)),
                                 xy=(adcpLoc[key]['location'][0],adcpLoc[key]['location'][1]),
                                 xycoords='data', xytext=(5, -15),
                                 textcoords='offset points', ha='left',
                                 color='white', fontsize=12)
    plt.figtext(.02, .02,
                "Notes:\n" +
                "'+' represents over-estimation whereas '-' represents under-estimation.\n" +
                "r2 and nrmse are respectively based on cubic signed speed and cubic speed. ",
                size='x-small')

    if savepath.strip() and fname.strip():
        if os.exists(savepath):
            fvcom.Plots._fig.savefig(savepath+fname, bbox_inches='tight')
            fvcom.Plots._fig.clear()
            plt.close(fvcom.Plots._fig)
    else:
        fvcom.Plots._fig.show()
        plt.show()