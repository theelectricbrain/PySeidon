%Brian Polagye
%November 11, 2011

function [res] = residual_currents(s, t, half_amp_T)

%Description: using pl66 filter, calculate residual currents

%Inputs:
%   - s: horizontal velocity (m/s)
%   - t: time
%   - half_amp_T: PL64 filter half-amplitude period (hours)

%Outputs:
%   - res: structure containing cropped time and velocity (filtered and unfiltered)

dt = (t(2)-t(1))*24;    %time interval of input series in hours

%calculate low pass filtered velocities
s_filter = NaN(size(s));
for i = 1:size(s_filter,2)
    s_filter(:,i) = pl66tn(s(:,i),dt,half_amp_T);
end

%eliminate end regions
crop_days = half_amp_T/24;    %number of days to crop at start and end of time series - junk data created by filter noise
indices = find(t-t(1) >= crop_days & t(end)-t >= crop_days);

%store cropped time and velocity (filtered and unfiltered) for future plotting
res.t = t(indices);
res.s = s(indices,:);
res.s_filter = s_filter(indices,:);

end