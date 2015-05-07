#!/usr/bin/python2.7
# encoding: utf-8
import cPickle as pickle
import numpy as np
import pandas as pd
import sys

# ALTERNATE VERSION FOR ANDY

def valTable(struct, filename, vars, debug=False, debug_plot=False):
    '''
    Takes validation data from the struct and saves it into a .csv file .

    Takes a single argument, a dictionary
    '''
    # initialize  lists
    val_dict = {}
    type, name, ovORuv, RMSE, CF, SD, POF, NOF, MDPO, MDNO, skill, r2, phase = \
    [], [], [], [], [], [], [], [], [], [], [], [], []
    num_tg = 1

    # append to the lists the stats from each site for each variable
    for var in vars:
        (type, name, ovORuv, RMSE, CF, SD, POF, NOF, MDPO, MDNO, skill, r2, phase) \
        = siteStats(struct, var, type, name,ovORuv, RMSE, CF, SD, POF,
                    NOF, MDPO, MDNO, skill, r2, phase, debug=False, debug_plot=False)

    # put stats into dict and create dataframe
    val_dict = {'Type':type, 'ovORuv':ovORuv, 'RMSE':RMSE, 'CF':CF, 'SD':SD, 'POF':POF,
                'NOF':NOF, 'MDPO':MDPO, 'MDNO':MDNO,  'skill':skill, 'r2':r2, 'phase':phase}

    table = pd.DataFrame(data=val_dict, index=name, columns=val_dict.keys())

    # export as .csv file
    out_file = '{}_val.csv'.format(filename)
    table.to_csv(out_file)
    return table

def siteStats(site, variable, type, name, ovORuv, RMSE, CF, SD, POF, NOF, MDPO, MDNO,
              skill, r2, phase, debug=False, debug_plot=False):
    '''
    Takes in the run (an array of dictionaries) and the type of the run (a
    string). Also takes in the list representing each statistic.
    '''
    if debug: print "siteStats..."
    # check if it's a tidegauge site
    if ((site['type'] != 'TideGauge') and (variable != 'tg')):
        stats = site['{}_val'.format(variable)]
        type.append(variable)
        name.append(site['name'].split('/')[-1].split('.')[0])

    elif ((site['type'] == 'TideGauge') and (variable == 'tg')):
        stats = site['tg_val']
        type.append('elev')
        name.append(site['name'].split('/')[-1].split('.')[0])

    # do nothing if a tidegauge is encountered but variable isn't tg
    else:
        print "---The variable tg is missing---"
        sys.exit()
   
    # add the statistics to the list, round to 2 decimal places
    ovORuv.append(stats['ovORuv'])
    RMSE.append(round(stats['RMSE'], 2))
    CF.append(round(stats['CF'], 2))
    SD.append(round(stats['SD'], 2))
    POF.append(round(stats['POF'], 2))
    NOF.append(round(stats['NOF'], 2))
    MDPO.append(stats['MDPO'])
    MDNO.append(stats['MDNO'])
    skill.append(round(stats['skill'], 2))
    r2.append(round(stats['r_squared'], 2))
    phase.append(stats['phase'])

    if debug: print "...siteStats done."

    return (type, name, ovORuv, RMSE, CF, SD, POF, NOF, MDPO, MDNO, skill, r2, phase)
