function [field_out, bin_max] = nan_AboveSurf(field_in,z,depth,varargin)
% function [field_out, bin_max] = nan_AboveSurf(field_in,z,depth,varargin)
%
% Nans all the values of field_in that are above the threshold 
% (fraction of surface height) 
%
% Justine McMillan
% Dec 5th, 2013

if nargin == 3
    threshold = 0.95;
elseif nargin == 4;
    threshold = varargin{1};
end

field_out = field_in;

for ii = 1:length(depth)
    zmax = threshold*depth(ii);
    [~ ,bin_max(ii)] = min(abs(z-zmax));
    if z(bin_max(ii))>zmax
        bin_max(ii) = bin_max(ii) - 1;
    end
    field_out(ii,bin_max(ii)+1:length(z)) = NaN;
end