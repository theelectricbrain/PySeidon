%Author: Brian Polagye
%Date: April 14, 2010

%Change Log:
% 5/2/10, BLP: integrated all subroutines into master function
% 5/3/10, BLP: fixed glitch in histogram plot routine
% 5/19/10, BLP: added histogram routine for vertical fit
% 6/9/10, BLP: separated plotting functions from data generating functions
% 6/14/10, BLP: added routine to generate signed speed if not included
% 7/21/10, BLP: added free surface to exported data in _stat
% 8/1/11, BLP: updated all routines, added new routines to calculate
%               probability distributions
% 8/23/11, BLP: moved all inputs to options structure input
% 11/6/11, BLP: added logarithmic profile function, modified power law
%               function
% 11/11/11, BLP: moved all functions to separate files so that they can be
%   called independently
% 11/11/11, BLP: streamlined output variables - particularly for profiles
% 05/31/12, JMM: added twoD option to take in 2D model data.
% 10/01/12, JMM: now calculating ensemble times within calc_ensemble routine
% 04/23/13, JMM: calculating ensemble direction differently (using
%                get_DirFromN function)
% 27/10/13, JMM: -Changed some of the variable names, used more structures, 
%                when computing ensembles, std dev is also output
%               - removed ensemble averaging part of script. (do this
%               before running this)
 


%Description: analysis toolpack for velocity information collected at
%   tidal energy sites by ADCP

%ADCP input must be in matlab format with the following structures -
%   assumed to be t rows (each corresponding to a time index) and z columns 
%   (each corresponding to specific depth bin)

%   time.mtime (t) - start time for ensemble

%   data.bins (z,1) - mid points of each measurement bin
%   data.dir_vel (t,z) - direction of velocity ensembles (north = 0 degrees, cw positive)
%   data.north_vel (t,z) - north-south component of ensemble velocity
%   data.east_vel (t,z) - east-west component of ensemble velocity
%   data.vert_vel (t,z) - vertical component of ensemble velocity
%   data.mag_signed_vel (t,z) - signed speed ensembles (ebb negative, flood
%       positive) (optional - if not present, will be automatically generated) 

%   pres.surf (t) - free surface elevation (m) (optional - if not present,
%   output value set to zero)

%Input
%   ADCP_file - name of file with ADCP data in above format
%   ADCP_path - full or relative path with ADCP data file
%   out_file - file to store results of from this routine
%   t_ens - ensemble durations for cycle statistics (s) - must be an integer multiple of the underlying series
%   z_target - height(s) relative to seabed (m) to perform calculations - an
%       empty value carries out the calculation for all vertical bins
%   max_bin - highest vertical bin to include - a value of zero uses NaN
%       fraction to determine this. AWAC with waves data will automatically
%       fail this test, requring that a maximum bin be specified.
%   options: structure containing configuration for characterization -
%       missing fields have defaults specified below

function exit_flag = characterize_resource_v2(ADCP_file, ADCP_path, out_file, t_ens, z_target, ...
    max_bin, options)

options = set_defaults(options);    %set default options
exit_flag = 0; 
%% Load data
disp(['Loading ' ADCP_path ADCP_file '...'])
load([ADCP_path ADCP_file])

%populate signed speed if not already present
if ~isfield(data, 'mag_signed_vel')
   flood_heading = input('Enter approximate direction of flood tide (0 = north, clockwise positive): ');
   flood_heading = flood_heading + [-90, +90];
   data.mag_signed_vel = sign_speed(data.east_vel, data.north_vel, data.mag_vel, data.dir_vel, flood_heading); 
end

% vertical bins
z = data.bins;

%convert target depths for analysis bin numbers
if ~isempty(z_target)
    for i = 1:length(z_target)
        bins(i) = find(abs(z-z_target(i))==min(abs(z-target(i))),1,'first');
    end
else
    bins = 1:length(z);
end

%Determine highest full bin and truncate bins above
if max_bin == 0
    nan_frac = sum(isnan(data.mag_signed_vel),1)/size(data.mag_signed_vel,1);
    z_max = find(nan_frac<=options.nan_max,1,'last');
else
    z_max = max_bin;
end
bins = bins(bins<=z_max);


% populate local variables -- ensembles
ens.t = time.mtime;             %time index
ens.d = data.dir_vel(:,1:z_max);           %current direction
ens.s = data.mag_signed_vel(:,1:z_max);   %signed speed
ens.u = data.east_vel(:,1:z_max);          %east velocity
ens.v = data.north_vel(:,1:z_max);         %north velocity
ens.w = data.vert_vel(:,1:z_max);          %vertical velocity
if exist('pres','var')
    ens.h = pres.surf;          %surface elevation
%elseif isfield(config,'depth')
%    h = config.depth;
else
    ens.h = 0;
end

% populate local variables -- standard deviations
stddev.d = data.dir_vel_std(:,1:z_max);           %current direction
stddev.s = data.mag_signed_vel_std(:,1:z_max);    %signed speed
stddev.u = data.east_vel_std(:,1:z_max);          %east vel
stddev.v = data.north_vel_std(:,1:z_max);         %north vel
stddev.h = pres.surf_std;                         %depth

%clear original variables
clear data time pres


% Check dt to ensure it is the same as t_ens
dt = mean(diff(ens.t))*3600*24;   %average measurement interval in seconds 
t_tol = 5; %tolerance of 5 s
if abs(dt-t_ens)>t_tol 
    error('This version does not compute ensembles. (Should update)')
    return
end

%% Directional characteristics
% - principal axes
% - standard deviation
% - directional asymmetry

%principal axis direction (ebb, flood, composite)
[fld.d_mean, ebb.d_mean, all.d_mean, ens.d_PA] = ...
    dir_PrincipalAxis(ens.d(:,bins), ens.s(:,bins), options.dir_min_s, ens.u(:,bins), ens.v(:,bins));  
%standard deviation from principal axis
[fld.d_sigma, ebb.d_sigma, all.d_sigma] = ...
    dir_std(ens.s(:,bins), options.dir_min_s, ens.d_PA, ens.u(:,bins), ens.v(:,bins), ...                     
    fld.d_mean, ebb.d_mean, options.show_scatter);
%directional asymmetry between ebb and flood
all.d_asym = abs(abs((ebb.d_mean - fld.d_mean))-180);   

%% Velocity characteristics
% - mean
% - max velocity
% - asymmetry 
% - vertical profiles - power law and logarithmic law
% - vertical shear
% - velocity probability distributions
% - residual currents

%mean velocity
[fld.s_mean, ebb.s_mean, all.s_mean] = s_mean(options.vel_min_s, ens.s(:,bins));
%max velocity
[fld.s_max, ebb.s_max, all.s_max] = s_max(ens.s(:,bins));           
%velocity asymmetry - ratio of mean ebb/mean flood, not residual current
all.s_asym = abs(ebb.s_mean./fld.s_mean);                           

% power law profile
if options.calc_profile_power
  [profile_power] = fitpower_profile(z, ens.s, options.profile_power_min_s);  
else
   profile_power = NaN;
end

% logarithmic profile - shear velocity and bottom roughness
if options.calc_profile_log
    [profile_log] = fitlog_profile(z, ens.s, options.rho, options.profile_log_min_s,...
        options.profile_log_z_start, options.profile_log_z_end);
else
   profile_log = NaN; 
end

%vertical shear
if options.twoD~=1
    [fld.ds_dz, ebb.ds_dz, all.ds_dz] = vertical_shear(z, options.vel_min_s, ens.s, bins);
end
%residual currents
if options.calc_residuals
    res_s = residual_currents(ens.s, ens.t, options.residual_half_amp_T);  
else
    res_s = NaN;
end
%probability density functions - speed, joint speed and direction
[pdf_v,jpdf_vd,v_bins,d_bins] = pdf_currents(ens.s, ens.d, options.hist_ds, options.hist_dd);

%Power characteristics (tabular investigation)
[fld.P_mean, ebb.P_mean, all.P_mean] = p_mean(options.power_min_s, options.rho, ens.s(:,bins));  
%ratio of mean ebb/mean flood power density
all.P_asym = abs(ebb.P_mean./fld.P_mean);                               

%Site characteristics (tabular)
if length(ens.h)>1
    all.h = nanmean(ens.h);
else
    all.h = h;
end

%% Misc. resource characteristics
% - current ellipses
% - cycle statistics (peak flood, peak ebb, cycle duration, cycle direction)

%isolate individual tidal cycles by identifying slack waters
[slack_ind, slack_t] = get_slacks(ens.s, ens.t, options.cycle_T_threshold);

[ellipse] = res_Ellipse(ens.u, ens.v, ens.w, slack_ind);        %calculate u,v,w ellipses over pairs of ebb/flood cycles
[cycle] = res_Cycle(slack_ind, slack_t, ens.s, ens.d_PA);       %calculate maximum amplitude, duration, and direction of tidal cycles

%% Generate summary plot of resource characteristics
if options.show_profiles
    plot_profiles(all, ebb, fld, z(bins))
    set(gcf,'name',ADCP_file)
end

%% Store results in characterization file

save([ADCP_path out_file],...
    'cycle', 'ellipse',...                                      %cycle statistics
    'z',...                                                     %vertical bins
    't_ens','bins','slack_t',...                                
    'profile*',...                                              %detailed statistics for vertical profiles
    'ens*','std*',...                                           %ensemble average properties, standard deviations
    'res_s','pdf_v','jpdf_vd','*bins',...                       %residual currents, probability distribution of currents
    'all','ebb','fld',...                                       %velocity, power, direction statitics    
    'lon','lat')                                                %latitude and longitude
end

