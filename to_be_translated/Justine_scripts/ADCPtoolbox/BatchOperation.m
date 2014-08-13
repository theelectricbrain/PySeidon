%Brian Polagye
%August 1, 2011

%Change Log:
%  9/14/11 (BLP) Additional commenting in line with code

%Description: routine to estimate turbine performance based on resource
%   characteristics at a site and defined operational parameters

clear
close all

%% Setup

batch_file_operation = 'Data_Operation.xlsx';    %batch file with operational parameters
batch_file_site = 'Data_Resource.xlsx';          %batch file with resource files
operation_proc = [1:3];                          %cases to process from operation file
out_path = 'Results\';                           %directory to store results

addpath('functions');

prop.rho = 1024;    %nominal sea water density (kg/m3)

%load batch data
[~,~,batch_list_site] = xlsread(batch_file_site);

for i = 2:size(batch_list_site,1)
    stat_files{i-1,1} = batch_list_site{i,4};        %names of stat files
    stat_paths{i-1,1} = batch_list_site{i,3};        %paths to stat files
end

[~,~,batch_list_operation] = xlsread(batch_file_operation);
for i = 2:size(batch_list_operation,1)
    operation(i-1).ID = batch_list_operation{i,1};
    operation(i-1).site = batch_list_operation{i,2};
    operation(i-1).hub = batch_list_operation{i,3};
    operation(i-1).D = batch_list_operation{i,4};
    operation(i-1).eta = batch_list_operation{i,5};
    operation(i-1).cutin_u = batch_list_operation{i,6};
    operation(i-1).rated_u = batch_list_operation{i,7};
    operation(i-1).yaw = batch_list_operation{i,8};
end

%% Run characterization routine

for i = operation_proc
    exit_flag = characterize_operation(stat_files{operation(i).site},...
        stat_paths{operation(i).site},operation(i),prop,out_path);
end
