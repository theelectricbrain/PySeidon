%Brian Polagye
%August 15, 2011

%Description: routine to plot operational metrics

clear
close all

%% Setup

batch_file_site = 'Data_Resource.xlsx';             %batch file with operational parameters
batch_file_operation = 'Data_Operation.xlsx';       %batch file with resource files
file_path = 'Results\';                             %directory containing results from operational estimates
plot_case = 1;                                      %operational case index in operation batch file

%plots to generate
show.operation = 1;         %plot trace of currents overlaying times operating and idle
show.operation_CDF = 1;     %plot CDF of operation, velocity, and power generation

%% Load data set

if length(plot_case) > 1, error('May only plot single cases'), end

[~,~,batch_list_site] = xlsread(batch_file_site);

for i = 2:size(batch_list_site,1)
    site_files{i-1,1} = batch_list_site{i,4};       %names of site files
    site_paths{i-1,1} = batch_list_site{i,3};       %names of site paths
end

[~,~,batch_list_operation] = xlsread(batch_file_operation);
for i = 2:size(batch_list_operation,1)
    site(i-1) = batch_list_operation{i,2};
end

file_name = [site_files{site(plot_case)}(1:end-4) '_turbine_' num2str(plot_case)];
load([file_path file_name])                                         %load operation data
disp(['Loading ' site_files{site(plot_case)} '...'])

%% Generate plots

%operational periods
if show.operation, plot_operation_period(operation, turbine.cutin_u, s, t); end

%operational cumulative distribution functions
if show.operation_CDF, plot_operation_CDF(generation, pdf_s, s_bins); end