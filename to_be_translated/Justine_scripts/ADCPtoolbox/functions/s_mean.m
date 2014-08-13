% Brian Polagye
% November 11, 2011

function [s_mean_fld, s_mean_ebb, s_mean_all] = s_mean(min_speed, s)

%Description: calculate mean speed for ebb, flood, and all conditions
%   can, optionally, exclude "slack" conditions

%Inputs:
%   - min_speed: minimum speed to be included in mean value
%   - s: horizontal velocity

%Ouputs:
%   - s_mean_fld: mean flood speed
%   - s_mean_ebb: mean ebb speed
%   - s_mean_all: mean speed over both ebb and flood

s_mean_all = zeros(size(s,2),1);
s_mean_ebb = zeros(size(s,2),1);
s_mean_fld = zeros(size(s,2),1);

%loop through all bins
for i = 1:size(s,2)
    s_mean_all(i) = nanmean(abs(s(abs(s(:,i))>=min_speed,i)));    
    s_mean_fld(i) = nanmean(s(s(:,i)>=min_speed,i));
    s_mean_ebb(i) = nanmean(s(s(:,i)<=-min_speed,i));
end

end