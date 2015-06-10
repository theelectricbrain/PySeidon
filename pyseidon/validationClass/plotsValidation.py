#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
#import seaborn
import matplotlib.pyplot as plt
#from pandas import DataFrame
import numpy as np
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf

def plotRegression(tidalStatClass, lr, save=False, out_f='', debug=False):
    """
    Plots a visualization of the output from linear regression,
    including confidence intervals for predictands and slope.

    If save is set to True, exports the plot as an image file to out_f.
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
    # df = DataFrame(data={'model': tidalStatClass.model.flatten(),
    #                      'observed':tidalStatClass.observed.flatten()})
    # seaborn.set(style="darkgrid")
    # color = seaborn.color_palette()[2]
    # g = seaborn.jointplot("model", "observed", data=df, kind="reg",
    #                       xlim=(df.model.min(), df.model.max()),
    #                       ylim=(df.observed.min(), df.observed.max()),
    #                       color=color, size=7)
    # plt.suptitle('Modeled vs. Observed {}: Linear Fit'.format(tidalStatClass.kind))

    if save:
        fig.savefig(out_f)
    else:
        fig.show()
        plt.show()

def plotData(tidalStatClass, graph='time', save=False, out_f='', debug=False):
    """
    Provides a visualization of the data.

    Takes an option which determines the kind of graph to be made.
    time: plots the model data against the observed data over time
    scatter : plots the model data vs. observed data

    If save is set to True, saves the image file in out_f.
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

    if save:
        fig.savefig(out_f)
    else:
        fig.show()
        plt.show()

def taylorDiagram(benchmarks, save=False, out_f='', debug=False):
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
        l, = ax.plot(benchmarks['NRMSE'][i]/100.0, benchmarks['r2'][i],
                    marker='$%d$' % (i+1), ms=10, ls='',
                    mfc=colors[i], mec=colors[i], label= benchmarks['Type'][i])
        samplePoints.append(l)
    t = np.linspace(0, np.pi/2)
    r = np.zeros_like(t) + 1.0
    ax.plot(t,r, 'k--', label='_')

    # Add NRMS contours, and label them
    rs,ts = np.meshgrid(np.linspace(smin, smax), np.linspace(0,np.pi/2))
    # Compute centered RMS difference
    rms = np.sqrt(1.0 + rs**2 - 2*rs*np.cos(ts))
    contours = ax.contour(ts, rs, rms, 5, colors='0.5')
    ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend
    fig.legend(samplePoints,[p.get_label() for p in samplePoints],
               numpoints=1, prop=dict(size='small'), loc='upper right')
    fig.tight_layout()

    if save:
        fig.savefig(out_f)
    else:
        fig.show()
        plt.show()