%Brian Polagye
%June 11, 2010

%Description: routine to generate plots of resource characteristics

%Change Log:
%   -7/12/10 (BLP) Fixed plot_cycle_histogram.m to display depth in title instead of bin number
%   -8/15/11 (BLP) Updated routines with modifications to characterize_ADCP
%   -8/15/11 (BLP) Added residual current plotting
%   -8/15/11 (BLP) Added joint pdf plotting for speed and direction
%   -9/8/11 (BLP) Added pcolor plotting of magnitude and direction
%   -11/7/11 (BLP) Added logarithmic profile plotting

clear
close all

addpath('functions');

%% Setup

z_target = 15;     %depth for table values (<1 = normalized z/D, > 1 = absolute height (m))

batch_file = 'Data_Resource.xlsx';  %batch file with resource file locations
plot_case = 1;                      %index in batch file to plot

%Choose plots to generate
show.current = 1;           %plot a trace of currents for a particular bin
show.ellipse = 1;           %plot a fortnight of current ellipses
show.cycle = 0;             %plot histograms of cycle properties
show.residual = 0;          %plot residual current contour
show.pdf = 1;               %plot probability distribution functions
show.stats = 0;             %plot statistical properties for all depths
show.power_profile = 0;     %plot vertical power law statistics
show.log_profile = 1;       %plot logarithmic law statistics
show.rose = 0;              %plot 4-point metric rose
show.pcolor = 0;            %plot pcolor of magnitude and direction

%% Open data set

if length(plot_case) > 1, error('May only plot single cases'), end

%load batch data
[~,~,batch_list] = xlsread(batch_file);

for i = 2:size(batch_list,1)    
    stat_files{i-1,1} = batch_list{i,4};        %names of output stat files
    ADCP_paths{i-1,1} = batch_list{i,3};        %paths to ADCP files
end

disp(['Loading ' stat_files{plot_case} '...'])
load([ADCP_paths{plot_case} stat_files{plot_case}])

if z_target <= 1
    z_target = z_target * all.h;
end

bin = find(abs(z-z_target)==min(abs(z-z_target)),1);      %determine bin closest to target height

%% Generate plots

%velocity trace at hub height
if show.current, plot_current(ens_s(:,bin), ens_t, slack_t, ens_w(:,bin)); end

%cycle ellipses over a fortinight (with vertical velocity)
if show.ellipse, plot_ellipse(ellipse, 1, bin, 1); end

%cycle histograms of duration, intensity, and direction
if show.cycle, plot_cycle_histogram(cycle, bin, z, 0.2, 0.5, 5); end

%residual currents
if show.residual, plot_residual(res_s,z,bin,10,10); end

%probability distribution functions
if show.pdf, plot_pdf(pdf_v, jpdf_vd, v_bins, d_bins, bin, z(bin)); end

%primary characteristics rose
if show.rose, plot_rose(all,bin,z(bin)); end

%statistics over water column
if show.stats, plot_profiles(all, ebb, fld, z(bins)), end

%power law profile details
if show.power_profile && isstruct(profile_power), plot_powerlaw_stats(profile_power,[],[],[]), end

%logarithmic law profile details
if show.log_profile && isstruct(profile_log), plot_loglaw_stats(profile_log, ens_t, [], [], [], []), end

%pcolor of magnitude and direction
if show.pcolor, plot_pcolor(ens_s, ens_t, ens_h, ens_d, z, 15), end