%Brian Polagye
%May 2, 2010

%Change Log:
%   6/10/10 (BLP) Moved statistics table to separate routine
%   8/31/11 (BLP) Moved all inputs to batch file, rather than embedding in
%       the sub-function
%   11/11/11 (BLP) Revamped options specification so that all function
%       options are exposed in BatchResource.m

%Description: routine to characterize resource from ADCP data

clear
close all

%% Setup
batch_file = 'Data_Resource.xlsx';   %batch file with resource file locations
char_proc = 1;                       %indices to process from resource file

addpath('functions');

t_ens = 5*60;   %duration for statistics calculation in seconds - end of eddy scale, beginning of aharmonic scale - must be an integer multiple of the underlying series

bins = [];      %depth relative to seabed (m) to calculate statistics - [] processes all depths (recommended)

%options: see discussion of parameters in characterize_ADCP.m
options = [];
% options.rho = 1024;                 %default = 1024 kg/m3
% options.calc_profile_power = 1;     %default = 0 (do not fit profile)
% options.profile_power_min_s = 0;    %default = 0 m/s (include all velocities)
% options.calc_profile_log = 1;       %default = 0 (do not fit profile)
% options.profile_log_min_s = 0;      %default = 0 m/s (include all velocities)
% options.profile_log_z_start = 10;   %default = 10 (starting bin for profile calculation)
% options.profile_log_z_end = [];     %default = [] (calculate over all bins)
% options.nan_max = 2e-2;             %default = 2e-2 (max fraction of NaNs in top layer)
% options.dir_min_s = 0.5;            %default = 0.5 m/s (minimum speed for direction calculation)
% options.vel_min_s = 0;              %default = 0 m/s (minimum speed to include in velocity calculations)
% options.power_min_s = 0;            %default = 0 m/s (minimum speed to include in power density calculations)
% options.cycle_T_threshold = 1;      %default = 1 hr (minimum cycle duration)
% options.hist_ds = 0.1;              %default = 0.1 m/s (speed bins for histograms)
% options.hist_dd = 1;                %default = 1 degree (degree bins for histograms)
% options.residual_half_amp_T = 40;   %default = 40 hrs (half-amplitude time for residual current filter)
% options.show_scatter = 1;           %default = 1 (show scatter plot)
% options.show_profiles = 1;          %default = 1 (show resource profiles)

%% Load batch data

%load batch data
[~,~,batch_list] = xlsread(batch_file);

for i = 2:size(batch_list,1)
    ADCP_files{i-1,1} = batch_list{i,2};        %names of ADCP input files (in .mat format)
    stat_files{i-1,1} = batch_list{i,4};        %names of output stat files
    ADCP_paths{i-1,1} = batch_list{i,3};        %paths to ADCP files
    max_bins{i-1,1} = batch_list{i,5};          %maximum bin to process - bypasses max bin identification routine
end

%% Run characterization routine
for i = char_proc
    exit_flag = characterize_resource(ADCP_files{i},ADCP_paths{i},stat_files{i},t_ens,bins,...
        max_bins{i},options);
    if exit_flag == 1, break; end    
end