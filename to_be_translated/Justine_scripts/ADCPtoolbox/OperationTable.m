%Brian Polagye
%August 2, 2011

%Description: make operational data table

clear
close all

%% Setup

batch_file_site = 'Data_Resource.xlsx';
batch_file_operation = 'Data_Operation.xlsx';
operation_proc = [1:3];   %turbine selection relative to batch file indices

file_path = 'Results\';

%% Load batch data

[~,~,batch_list_site] = xlsread(batch_file_site);

for i = 2:size(batch_list_site,1)
    site_files{i-1,1} = batch_list_site{i,4};        %names of site files
end

[~,~,batch_list_operation] = xlsread(batch_file_operation);
for i = 2:size(batch_list_operation,1)
    site(i-1) = batch_list_operation{i,2};
end

%% Load data and extract for table

for i = operation_proc
    file_name = [site_files{site(i)}(1:end-4) '_turbine_' num2str(i)];
    load([file_path file_name])
    out_table(i,2) = generation.power_avg*1000;
    out_table(i,4) = generation.cf;
    out_table(i,3) = generation.power_max*1000;
    out_table(i,5) = operation.f_op;
    out_table(i,1) = turbine.angle;
end
