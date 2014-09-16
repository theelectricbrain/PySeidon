from __future__ import division
import numpy as np
#from rawADCPclass import rawADCP
from datetime import datetime
from datetime import timedelta
import scipy.io as sio
import scipy.interpolate as sip
import matplotlib.pyplot as plt
import seaborn

def date2py(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + \
        timedelta(days=matlab_datenum%1) - timedelta(days = 366)

    return python_datetime


def py2date(dt):
   mdn = dt + timedelta(days = 366)
   frac_seconds = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
   return mdn.toordinal() + frac_seconds + frac_microseconds

def calc_ensemble(x, ens, ens_dim):

    #initialize input
    ens = int(ens)
    #x = x[:, None]

    if ens_dim == 1:
        ens_size = np.floor(x.shape[0]/60)
    else:
        pass

    #x_ens = np.empty((ens_size, 1, ens))
    x_ens = np.empty((ens_size, ens))
    x_ens[:] = np.nan

    for j in xrange(ens):
        if ens_dim == 1:
            ind_ens = np.arange(j, x.shape[0] - (ens - j), ens)
            #x_ens[..., j] = x[ind_ens]
            x_ens[..., j] = x[ind_ens]

        else:
            pass

    #x_ens = np.nanmean(x_ens, axis=2)
    x_ens = np.nanmean(x_ens, axis=1)
    return x_ens


def rotate_coords(x, y, theta):
    '''
    Similar to "rotate_to_channelcoords.m" code,
    theta is now the angle
    between the old axis and the new x-axis (CCw is positive)
    '''

    xnew = x * np.cos(theta) + y * np.sin(theta)
    ynew = -x * np.sin(theta) + y * np.cos(theta)

    return xnew, ynew

def rotate_to_true(X, Y, theta=-19):
    '''
    % X,Y are the X and Y coordinates (could be speeds) relative to magnetic
    % north -- inputs can be vectors
    % x,y are the coordinates relative to true north
    % This function assumes the measured location is Nova Scotia where the
    % declination angle is -19 degrees.
    %
    % Sept 29, 2012: Changed print statement
    %
    % Sept 20, 2012: Modified the function to allow for theta to be input.
    % Default will remain at -19 degrees, but this may not be accurate for all
    % places in Nova Scotia.
    '''

    print 'Rotating velocities to be relative to true north (declination = {0})'.format(theta)

    Theta = theta * np.pi / 180

    x = X * np.cos(Theta) + Y * np.sin(Theta)
    y = -X * np.sin(Theta) + Y * np.cos(Theta)

    return x, y


def get_DirFromN(u,v):
    '''
    #This function computes the direction from North with the output in degrees
    #and measured clockwise from north.
    #
    # Inputs:
    #   u: eastward component
    #   v: northward component
    '''

    theta = np.arctan2(u,v) * 180 / np.pi

    ind = np.where(theta<0)
    theta[ind] = theta[ind] + 360
    return theta

def sign_speed(u_all, v_all, s_all, dir_all, flood_heading):

    if type(flood_heading)==int:
        flood_heading += np.array([-90, 90])

    s_signed_all = np.empty(s_all.shape)
    s_signed_all.fill(np.nan)

    PA_all = np.zeros(s_all.shape[-1])
    for i in xrange(s_all.shape[-1]):
        u = u_all[:, i]
        v = v_all[:, i]
        dir = dir_all[:, i]
        s = s_all[:, i]

        #determine principal axes - potentially a problem if axes are very kinked
        #   since this would misclassify part of ebb and flood
        PA, _ = principal_axis(u, v)
        PA_all[i] = PA

        # sign speed - eliminating wrap-around
        dir_PA = dir - PA

        dir_PA[dir_PA < -90] += 360
        dir_PA[dir_PA > 270] -= 360

        #general direction of flood passed as input argument
        if flood_heading[0] <= PA <= flood_heading[1]:
            ind_fld = np.where((dir_PA >= -90) & (dir_PA<90))
            s_signed = -s
            s_signed[ind_fld] = s[ind_fld]
        else:
            ind_ebb = np.where((dir_PA >= -90) & (dir_PA<90))
            s_signed = s
            s_signed[ind_ebb] = -s[ind_ebb]

        s_signed_all[:, i] = s_signed

    return s_signed_all, PA_all

def principal_axis(u, v):

    #create velocity matrix
    U = np.vstack((u,v)).T
    #eliminate NaN values
    U = U[~np.isnan(U[:, 0]), :]
    #convert matrix to deviate form
    rep = np.tile(np.mean(U, axis=0), [len(U), 1])
    U -= rep
    #compute covariance matrix
    R = np.dot(U.T, U) / (len(U) - 1)

    #calculate eigenvalues and eigenvectors for covariance matrix
    lamb, V = np.linalg.eig(R)
    #sort eignvalues in descending order so that major axis is given by first eigenvector
    # sort in descending order with indices
    ilamb = sorted(range(len(lamb)), key=lambda k: lamb[k], reverse=True)
    lamb = sorted(lamb, reverse=True)
    # reconstruct the eigenvalue matrix
    lamb = np.diag(lamb)
    #reorder the eigenvectors
    V = V[:, ilamb]

    #rotation angle of major axis in radians relative to cartesian coordiantes
    ra = np.arctan2(V[0,1], V[1,1])
    #express principal axis in compass coordinates
    # WES_COMMENT: may need to change this, cause in original is -ra
    PA = ra * 180 / np.pi + 90
    #variance captured by principal
    varxp_PA = np.diag(lamb[0]) / np.trace(lamb)

    return PA, varxp_PA


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


def save_FlowFile_BPFormat(fileinfo, adcp, rbr, params, options, debug=False):

    comments = ['data is in Polagye Tools format',
                'data.east_vel and data.north_vel are relative to true north',
               'The parameters were set by ' + fileinfo['paramfile']]

    day1 = date2py(adcp['mtime'][0][0])
    print day1
    #date_time = [date2py(tval[0]) for tval in adcp.mtime[:]]
    datenum = datetime(day1.year,1,1) + timedelta(365)
    datenum = datenum.toordinal()

    yd = adcp['mtime'][:].flatten() - datenum
    tind = np.where((yd > params['tmin']) & (yd < params['tmax']))[0]

    pres = {}
    time = {}
    time['mtime'] = adcp['mtime'][:].flatten()[tind]
    dt = np.nanmean(np.diff(time['mtime']))

    if not rbr:
        print 'Depths measured by ADCP not yet coded.'
        comments.append('Depths as measured by ADCP')
    else:
        print 'Ensemble averaging rbr data'
        comments.append('Depths as measured by RBR sensor')

        nens = round(dt/(rbr.mtime[1] - rbr.mtime[0]))
        temp = np.arange(rbr.mtime[nens/2-1], rbr.mtime[-1-nens/2], dt)
        #temp2 = np.r_[rbr.mtime[nens/2-1]: rbr.mtime[-1-nens/2]: dt]

        mtimeens = np.arange(rbr.mtime[nens/2-1], rbr.mtime[-1-nens/2], dt)
        mtimeens = mtimeens + params['rbr_hr_offset'] / 24
        depthens = calc_ensemble(rbr.depth, nens, 1)

        temp = sip.interp1d(mtimeens, depthens, kind='linear')

        pres['surf']= temp(time['mtime']) + params['dabPS']

        if debug:
            # Load in matlab values for testing
            filename = './140703-EcoEII_database/scripts_examples/mtime.mat'
            mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
            matTimes = mat['mtimeens']
            filename = './140703-EcoEII_database/scripts_examples/dt.mat'
            mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
            matdt = mat['dt']


            filename = './140703-EcoEII_database/scripts_examples/depthens.mat'
            mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
            matdepthens = mat['depthens']

            filename = './140703-EcoEII_database/scripts_examples/time.mat'
            mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
            matmtime = mat['mtime']

            print matTimes.shape
            print temp - matTimes
            print temp2 - matTimes
            print dt - matdt
            print depthens - matdepthens
            print 'time'
            print time['mtime'] - matmtime

    ## zlevels
    data = {}
    z = adcp['config']['ranges'][:] + params['dabADCP']
    z = z.flatten()
    zind = np.where((z > params['zmin']) & (z < params['zmax']))[0]
    data['bins'] = z[zind]

    ## Currents
    data['vert_vel'] = adcp['vert_vel'][:][tind][:, zind]
    data['error_vel'] = adcp['error_vel'][:][tind][:, zind]

    # If compass wasn't calibrated
    if 'hdgmod' in params:
        adcp['east_vel'][:], adcp['north_vel'][:] = rotate_coords(adcp['east_vel'][:],
                                                              adcp['north_vel'][:],
                                                              params['hdgmod'])

        comments.append('East and north velocity rotated by params.hdgmod')

    # Rotate east_vel and north_vel to be relative to true north
    data['east_vel'], data['north_vel'] = \
        rotate_to_true(adcp['east_vel'][:][tind][:, zind],
                       adcp['north_vel'][:][tind][:, zind],
                       params['declination'])

    # Direction
    data['dir_vel'] = get_DirFromN(data['east_vel'],data['north_vel'])

    # Signed Speed
    spd_all = np.sqrt(data['east_vel']**2+data['north_vel']**2)

    # Determine flood and ebb based on principal direction (Polagye Routine)
    print 'Getting signed speed (Principal Direction Method) -- used all speeds'
    s_signed_all, PA_all = sign_speed(data['east_vel'], data['north_vel'],
                                      spd_all, data['dir_vel'], params['flooddir'])

    data['mag_signed_vel'] = s_signed_all

    if options['showRBRavg'] or debug:
        print 'Plotting RBR vs average'
        plt.plot(rbr.mtime + params['rbr_hr_offset'] / 24, rbr.depth+params['dabPS'],
                    label='RBR')
        plt.plot(time['mtime'], pres['surf'], 'r', label='AVG')
        plt.xlabel('Time')
        plt.ylabel('Elevation')
        plt.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=plt.gcf().transFigure)

        plt.show()

    if options['showPA'] or debug:
        print 'Plotting PA vs mean'
        plt.plot(PA_all, data['bins'], label='PA')
        plt.plot(np.array([PA_all[0], PA_all[-1]]),
                 np.array([np.mean(pres['surf']), np.mean(pres['surf'])]),
                 label='mean')

        plt.xlabel('Principal Axis Direction\n(clockwise from north)')
        plt.ylabel('z (m)')
        plt.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=plt.gcf().transFigure)
        plt.show()

    ## save
    lon = params['lon']
    lat = params['lat']

    outfile = fileinfo['outdir'] + fileinfo['flowfile']
    print 'Saving data to {0}'.format(outfile)

    saveDict = {'data':data, 'pres':pres, 'time':time, 'lon':lon, 'lat':lat,
                'params':params, 'comments':comments}
    #save(outfile,'data','pres','time','lon','lat','params','Comments')

    ## Save metadata
    #metadata.progname=[mfilename('fullpath')];
    #metadata.date = datestr(now);
    #metadata.paramfile = fileinfo.paramfile;
    #save(outfile,'metadata','-append')
    return saveDict



if __name__ == '__main__':
    filename = '140703-EcoEII_database/data/GP-120726-BPd_raw.mat'
    data = rawADCP(filename)
    rawdata = rawADCP(filename)
    #adcp = Struct(**data.adcp)
    #rawADCP = data.adcp
    adcp = data.adcp
    #params = Struct(**data.saveparams)
    params = data.saveparams
    rbr = Struct(**data.rbr)

#    save_FlowFile_BPFormat(data.fileinfo, data.adcp, data.rbr,
#                           data.saveparams, data.options)

    saveDict = \
    save_FlowFile_BPFormat(data.fileinfo, adcp, rbr,
                           params, data.options)
