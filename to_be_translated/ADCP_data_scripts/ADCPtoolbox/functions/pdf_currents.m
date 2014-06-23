%Brian Polagye
%November 11, 2011

function [pdf_v, jpdf_vd, v_bins, d_bins] = pdf_currents(s, d, dv, dd)

%Description: calculate probability distribution of horizontal velocity and
%joint distribution of horizontal velocity and direction

%Inputs:
%   - s: horizontal velocity (m/s)
%   - d: direction (degrees - relative to principal axes)
%   - dv: horizontal velocity bin intervals (m/s)
%   - dd: direction bin intervals (degrees)

%Outputs:
%   - pdf_v: probability distribution of horizontal velocity (s x z)
%   - jpdf_vd: joint probability distribution of horizontal velocity and direction (s x d x z) 
%   - v_bins: horizontal velocity bins (m/s)
%   - d_bins: directional bins (deg)

%determine range of values for probability distributions
v_min = floor(min(min(s)));
v_max = ceil(max(max(s)));

v_bins = v_min:dv:v_max;

%set minimum and maximum values to span complete set of compass coordinates
d_min = 0;
d_max = 360;

d_bins = d_min:dd:d_max;

%initialize output variables
pdf_v = zeros(length(v_bins),size(s,2));
jpdf_vd = zeros(length(v_bins),length(d_bins),size(s,2));

%generate distributions for each bin
ctrs{1} = v_bins;
ctrs{2} = d_bins;
for i = 1:size(s,2)
    %calculate counts by bins
    jpdf_vd(:,:,i) = hist3([s(:,i) d(:,i)],ctrs);
    pdf_v(:,i) = sum(jpdf_vd(:,:,i),2);
    
    %convert to frequency distribution
    jpdf_vd(:,:,i) = jpdf_vd(:,:,i)/sum(sum(jpdf_vd(:,:,i)));
    pdf_v(:,i) = pdf_v(:,i)/sum(pdf_v(:,i));    
end

end