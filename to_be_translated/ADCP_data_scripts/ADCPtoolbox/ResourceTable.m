%Description: create resource characteristics table
%
% Modified from Brian Polagye's ResourceTable.m


clear

%% Setup                         

batch_file = 'Data_Resource.xlsx';      %batch file with resource file locations
char_proc = [1];                        %indices from resource file to include in table

z_target = 15;       %depth for table values (<1 = normalized z/D, > 1 = absolute height (m))

%% Load batch data

[~,~,batch_list] = xlsread(batch_file);

for i = 2:size(batch_list,1)
    stat_files{i-1,1} = batch_list{i,4};        %names of stat files
    stat_paths{i-1,1} = batch_list{i,3};        %paths to stat files
    lats(i-1,1) = batch_list{i,6};              %latitude
    lons(i-1,1) = batch_list{i,7};              %longitude
end

%% Extract information for table

cnt = 1;
for i = char_proc
   disp(['Loading ' stat_files{i} '...'])
   load([stat_paths{i} stat_files{i}])
   
   %set target depth if not already specified
   if z_target <= 1
       z_target = z_target * all.h;
   end
   
   %identify bin closest to target depth
   bin = find(abs(z-z_target)==min(abs(z-z_target)),1,'first');
   
   %set reference location for comparative values
   if cnt == 1
      [x1, y1, ~] = wgs2utm(lats(i),lons(i)); 
   end
   
   [x,y,~] = wgs2utm(lats(i),lons(i));
   
   table(cnt,1) = (( (x-x1)^2 + (y-y1)^2 )^0.5)/1000;   %distance from reference location
   table(cnt,2) = z_target;                             %target depth (m)
   table(cnt,3) = all.P_mean(bin);                      %mean power density (kW/m^2) of all tides
   table(cnt,4) = all.P_asym(bin);                      %mean ebb power density / mean flood power density
   table(cnt,5) = all.d_asym(bin);                      %direction of ebb - direction of flood - 180 degrees
   table(cnt,6) = all.d_sigma(bin);                     %standard deviation from principal axes on ebb and flood
   table(cnt,7) = all.s_max(bin);                       %maximum observed velocity (m/s)
   
   cnt = cnt + 1;
end