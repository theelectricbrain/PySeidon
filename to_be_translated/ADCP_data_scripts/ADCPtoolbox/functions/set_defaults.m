%Brian Polagye
%November 11, 2011
%
% Updates: 
%   - Oct 28, 13,JMM: Added option to skip the calculation of residual
%   currents

%Description: set defaults for all user-programmable options

function [options] = set_defaults(options)

%% Seawater properties

%nominal seawater density (kg/m3)
if ~isfield(options,'rho'), options.rho = 1024; end

%% 2D
if ~isfield(options,'twoD'), options.twoD = 0; end

%% Vertical profiles

%Power law profile calculation
%   1 = calculate best fit of power law to profile
%   0 = do not calculate 
% Depending on series legnth, routine may be time intensive
if ~isfield(options,'calc_profile_power'), options.calc_profile_power = 0; end

%Minimum speed to calculate power law fit - will speed up calculations
if ~isfield(options,'profile_power_min_s'), options.profile_power_min_s = 0; end

%Log law profile calculation
%   1 = calculate best fit of log layer to profile
%   0 = do not calculate  
% Routine is time intensive - for 90 days series, expect several hours
if ~isfield(options,'calc_profile_log'), options.calc_profile_log = 0; end

%Minimum speed to calculate log law fit - will speed up calculations
if ~isfield(options,'profile_log_min_s'), options.profile_log_min_s = 0; end

%Starting depth bin for log profile fit
if ~isfield(options,'profile_log_z_start'), options.profile_log_z_start = 10; end

%Ending depth bin for log profile fit
% defaults to highest complete bin in velocity profile
if ~isfield(options,'profile_log_z_end'), options.profile_log_z_end = []; end

%% Valid bins

%maximum allowable fraction of NaN values for bin to still be considered complete
if ~isfield(options,'nan_max'), options.nan_max = 2e-2; end

%% Minimum speeds

%minimum speed to be considered ebb or flood, rather than slack
%   - used to exclude directional aberations around slack water
if ~isfield(options,'dir_min_s'), options.dir_min_s = 0.5; end
            
%mininum speed to be included in mean speed analysis
%   - default to include all velocities
if ~isfield(options,'vel_min_s'), options.vel_min_s = 0; end

%minimum speed to be included in power density analysis
%   - default to include all velocities
if ~isfield(options,'power_min_s'), options.power_min_s = 0; end

%% Cycle analysis
if ~isfield(options,'slacks'), options.slacks = 1; end
if ~isfield(options,'ellipse'), options.ellipse = 1; end
if ~isfield(options,'cyclestats'), options.cyclestats = 1; end

%minimum cycle duration (in hours)
%   - used to screen for "cycles" that are merely oscillations around slack water
if ~isfield(options,'cycle_T_threshold'), options.cycle_T_threshold = 1; end

%% Histogram bins

%interval (m/s) for probability distribution of speed
if ~isfield(options,'hist_ds'), options.hist_ds = 0.1; end

%interval (degrees) for joint probability distributions of speed and direction
if ~isfield(options,'hist_dd'), options.hist_dd = 1; end

%% Residual currents

%PL64 filter half-amplitude period (hours) 
%   - default value = 33 - per Beardsley et al.
%   - a value of 40 gives better results for strong tidal currents
if ~isfield(options,'calc_residuals'), options.calc_residuals = 0; end
if ~isfield(options,'residual_half_amp_T'), options.residual_half_amp_T = 40; end

%% Visualization

%Show scatter plots with principal axis fits
%   1 = show scatter plots (debugging option)
%   0 = hide these
if ~isfield(options,'show_scatter'), options.show_scatter = 1; end

%Show profile plots for characteristics
%   1 = at end of routine, show plot profiles for characteristics
%   0 = do not generate plot
if ~isfield(options,'show_profiles'), options.show_profiles = 1; end

end

