#!/usr/bin/python2.7
# encoding: utf-8

import seaborn
import matplotlib.pyplot as plt
from pandas import DataFrame
import numpy as np

def plotRegression(tidalStatClass, lr, save=False, out_f='', debug=False):
    """
    Plots a visualization of the output from linear regression,
    including confidence intervals for predictands and slope.

    If save is set to True, exports the plot as an image file to out_f.
    """
    if debug : print "Plotting linear regression"
    df = DataFrame(data={'model': tidalStatClass.model.flatten(),
                            'observed':tidalStatClass.observed.flatten()})
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
    seaborn.set(style="darkgrid")
    color = seaborn.color_palette()[2]
    g = seaborn.jointplot("model", "observed", data=df, kind="reg",
                          xlim=(df.model.min(), df.model.max()),
                          ylim=(df.observed.min(), df.observed.max()),
                          color=color, size=7)
    plt.suptitle('Modeled vs. Observed {}: Linear Fit'.format(tidalStatClass.kind))

    if save:
        fig.savefig(out_f)
    else:
        fig.show()