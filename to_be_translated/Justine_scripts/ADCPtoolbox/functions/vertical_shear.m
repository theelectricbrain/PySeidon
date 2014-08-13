%Brian Polagye
%November 11, 2011

function [ds_dz_fld, ds_dz_ebb, ds_dz_all, ds_dz] = vertical_shear(z, min_s, s, bins)

%Description: calculate mean vertical shear (du/dz) for ebb and flood

%Inputs:
%   - z: vertical bin centers
%   - min_s: minimum speed to include for shear calculations - can be used to exclude shear anomalies around slack water
%   - s: horizontal velocity
%   - bins: vertical bins to perform calculation over

%Outputs:
%   - ds_dz_fld: mean vertical shear on flood
%   - ds_dz_ebb: mean vertical shear on ebb
%   - ds_dz_all: mean vertical shear over all tides
%
% Modified by Justine McMillan (Dec 5th, 2013)
%   - now outputs ds_dz (unmeaned)

z_max = size(s,2);
dz = z(2) - z(1);

%initialize variables
ds_dz_all = zeros(z_max,1);
ds_dz_ebb = zeros(z_max,1);
ds_dz_fld = zeros(z_max,1);

ds_dz_all([1,end]) = NaN;
ds_dz_ebb([1,end]) = NaN;
ds_dz_fld([1,end]) = NaN;

ds_dz = zeros(size(s,1),1);
ds_dz(:,1) = NaN;
ds_dz(:,z_max) = NaN;

%calculate vertical shear using central difference
for i = 2:z_max-1;
    ds_dz(:,i) = abs(s(:,i+1)-s(:,i-1))/(2*dz);
    
    ind_fld = find(s(:,i)>=min_s & ~isnan(ds_dz(:,i)));
    ind_ebb = find(s(:,i)<=-min_s & ~isnan(ds_dz(:,i)));

    ds_dz_ebb(i) = mean(ds_dz(ind_ebb,i));      %mean vertical shear on ebb (m/s per m)
    ds_dz_fld(i) = mean(ds_dz(ind_fld,i));      %mean vertical shear on flood (m/s per m)
    ds_dz_all(i) = (ds_dz_ebb(i)*length(ind_ebb)+ ...
        ds_dz_fld(i)*length(ind_fld))/(length(ind_fld)+length(ind_ebb));   %mean shear over full cycle   
end

%truncate to bins of interest
ds_dz_ebb = ds_dz_ebb(bins);
ds_dz_fld = ds_dz_fld(bins);
ds_dz_all = ds_dz_all(bins);
ds_dz     = ds_dz(:,bins);

end