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
% 04/23/12, JMM: calculating ensemble direction differently (using
%                get_DirFromN function)

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

function exit_flag = characterize_resource(ADCP_file, ADCP_path, out_file, t_ens, z_target, ...
    max_bin, options)

options = set_defaults(options);    %set default options
  
%% Load data
disp(['Loading ' ADCP_path ADCP_file '...'])
load([ADCP_path ADCP_file])

%populate signed speed if not already present
if ~isfield(data, 'mag_signed_vel')
   flood_heading = input('Enter approximate direction of flood tide (0 = north, clockwise positive): ');
   flood_heading = flood_heading + [-90, +90];
   data.mag_signed_vel = sign_speed(data.east_vel, data.north_vel, data.mag_vel, data.dir_vel, flood_heading); 
end

%populate local variables
z = data.bins;              %vertical bins
t = time.mtime;             %time index
d = data.dir_vel;           %current direction
s = data.mag_signed_vel;    %signed speed
u = data.east_vel;          %east velocity
v = data.north_vel;         %north velocity
w = data.vert_vel;          %vertical velocity
if exist('pres','var')
    h = pres.surf;          %surface elevation
%elseif isfield(config,'depth')
%    h = config.depth;
else
    h = 0;
end


%clear original variables
clear data config QA anc time pres

%% Process raw data to working quantities for characterization

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
    nan_frac = sum(isnan(s),1)/size(s,1);
    z_max = find(nan_frac<=options.nan_max,1,'last');
else
    z_max = max_bin;
end

s = s(:,1:z_max);
d = d(:,1:z_max);
u = u(:,1:z_max);
v = v(:,1:z_max);
w = w(:,1:z_max);

bins = bins(bins<=z_max);

%Ensemble average measured quantities (generally used to minimize influence
%of turbulent length scales)

dt = (t(2)-t(1))*3600*24;   %measurement interval in seconds

[ens_t, ens_s, ens_d, ens_u, ens_v, ens_w, ens_h, exit_flag] = EnsembleData(t_ens, dt, t, s, d, u, v, w, h);

if exit_flag == 1, 
    error('Exiting, not completing resource characterization')
    return
end

%% Directional characteristics
% - principal axes
% - standard deviation
% - directional asymmetry

%principal axis direction (ebb, flood, composite)
[fld.d_mean, ebb.d_mean, all.d_mean, ens_d_PA] = ...
    dir_PrincipalAxis(ens_d(:,bins), ens_s(:,bins), options.dir_min_s, ens_u(:,bins), ens_v(:,bins));  
%standard deviation from principal axis
[fld.d_sigma, ebb.d_sigma, all.d_sigma] = ...
    dir_std(ens_s(:,bins), options.dir_min_s, ens_d_PA, ens_u(:,bins), ens_v(:,bins), ...                     
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
[fld.s_mean, ebb.s_mean, all.s_mean] = s_mean(options.vel_min_s, ens_s(:,bins));
%max velocity
[fld.s_max, ebb.s_max, all.s_max] = s_max(ens_s(:,bins));           
%velocity asymmetry - ratio of mean ebb/mean flood, not residual current
all.s_asym = abs(ebb.s_mean./fld.s_mean);                           

% power law profile
if options.calc_profile_power
  [profile_power] = fitpower_profile(z, ens_s, options.profile_power_min_s);  
else
   profile_power = NaN;
end

% logarithmic profile - shear velocity and bottom roughness
if options.calc_profile_log
    [profile_log] = fitlog_profile(z, ens_s, options.rho, options.profile_log_min_s,...
        options.profile_log_z_start, options.profile_log_z_end);
else
   profile_log = NaN; 
end

%vertical shear
if options.twoD~=1
    [fld.ds_dz, ebb.ds_dz, all.ds_dz] = vertical_shear(z, options.vel_min_s, ens_s, bins);
end
%residual currents
res_s = residual_currents(ens_s, ens_t, options.residual_half_amp_T);                        
%probability density functions - speed, joint speed and direction
[pdf_v,jpdf_vd,v_bins,d_bins] = pdf_currents(ens_s, ens_d, options.hist_ds, options.hist_dd);

%Power characteristics (tabular investigation)
[fld.P_mean, ebb.P_mean, all.P_mean] = p_mean(options.power_min_s, options.rho, ens_s(:,bins));  
%ratio of mean ebb/mean flood power density
all.P_asym = abs(ebb.P_mean./fld.P_mean);                               

%Site characteristics (tabular)
if length(ens_h)>1
    all.h = nanmean(ens_h);
else
    all.h = h;
end

%% Misc. resource characteristics
% - current ellipses
% - cycle statistics (peak flood, peak ebb, cycle duration, cycle direction)

%isolate individual tidal cycles by identifying slack waters
[slack_ind, slack_t] = get_slacks(ens_s, ens_t, options.cycle_T_threshold);

[ellipse] = res_Ellipse(ens_u, ens_v, ens_w, slack_ind);        %calculate u,v,w ellipses over pairs of ebb/flood cycles
[cycle] = res_Cycle(slack_ind, slack_t, ens_s, ens_d_PA);       %calculate maximum amplitude, duration, and direction of tidal cycles

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
    'ens_*',...                                                 %ensemble average properties
    'res_s','pdf_v','jpdf_vd','*bins',...                       %residual currents, probability distribution of currents
    'all','ebb','fld',...                                       %velocity, power, direction statitics    
    'lon','lat')                                                %latitude and longitude
end

%% EnsembleData
function [ens_t, ens_s, ens_d, ens_u, ens_v, ens_w, ens_h, exit_flag] = EnsembleData(t_ens, dt, t, s, d, u, v, w, h)

%Description: ensemble average (time and space) measurements to smooth effect of turbulence

%Inputs:
%   - t_ens: ensemble time in seconds
%   - dt: time interval in seconds
%   - t: time stamps (matlab format)
%   - s: horizontal velocity (signed ebb and flood)
%   - d: horozintal velocity direction
%   - u: east velocity
%   - v: north velocity
%   - w: vertical velocity
%   - h: free surface

%Outputs: ensemble average inputs, time stamps for start of ensemble

exit_flag = 0;  %default to successful execution
ens_tol = 0.1;  %maximum allowable mismatch between number of integer samples in ensemble

ens = t_ens/dt;   %number of points in ensemble average

if abs(round(ens)-ens) > ens_tol || ens+ens_tol < 1
    exit_flag = 1;
    disp('Error: mismatch between measurement time step and specified ensemble averaging time greater than tolerance')
    ens_t = NaN;
    ens_s = NaN;
    ens_d = NaN;
    ens_u = NaN;
    ens_v = NaN;
    ens_w = NaN;
    ens_h = NaN;
else
    ens = round(ens);
    disp(['Number of points in ensemble average: ',num2str(ens)])
    
    if ens == 1
        disp('Not computing ensembles')
        ens_t = t;
        ens_s = s;
        ens_d = d;
        ens_u = u;
        ens_v = v;
        ens_w = w;
        ens_h = h;
    else
        %ens_t = t(1):t_ens/(60*60*24):t(end);    %time stamp is start of averaging period
        t_ens_d = t_ens/(60*60*24);
        ens_t = t(1)+t_ens_d/2:t_ens_d:t(end);%-t_ens_d/2; %time stamps are the middle of the interval
        %temp = datevec(ens_t);
        %ens_t = datenum(temp(:,1),temp(:,2),temp(:,3),temp(:,4),round(temp(:,5)+temp(:,6)/60),0);   %clean up time stamps to integer minutes
        %ens_t = datenum(temp(:,1),temp(:,2),temp(:,3),temp(:,4),temp(:,5),round(temp(:,6)));   %clean up time stamps to integer seconds
     
        ens_s = calc_ensemble(s,ens,1);     %horizontal velocity
        ens_u = calc_ensemble(u,ens,1);     %east velocity
        ens_v = calc_ensemble(v,ens,1);     %west velocity
        %ens_d = calc_ensemble(d,ens,1);     %direction
        ens_d = get_DirFromN(ens_u,ens_v); %Justine Modification (gives weird angles otherwise)
        ens_w = calc_ensemble(w,ens,1);     %vertical velocity
        if length(h)>1
            ens_h = calc_ensemble(h,ens,2);     %free surface elevation
        else
            ens_h = h;
        end
        
        %prune time series for incomplete ensembles
        if length(ens_t) > size(ens_s,1)
            ens_t = ens_t(1:size(ens_s,1));
        end
    end
end

end
