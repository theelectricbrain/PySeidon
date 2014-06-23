%Brian Polagye
%November 11, 2011

%Descrption: function to generate statistics about each tidal cycle within measurement period (duration and maximum amplitude)

function [cycle] = res_Cycle(slack_ind, slack_t, s, d)
%Inputs:
%   slack_ind: indices in time series corresponding to slack water
%   slack_t: time of slack water
%   s: time series of horizontal velocity
%   d: time series of horizontal velocity direction
%       - note: d_PA_ens removes risk of distortions around 0/360, but must
%           be corrected back to earth coordinates from principal axis
%           coordinates

%Outputs:
%   cycle: structure containing duration, intensity, and direction for each
%       cycle in input time series

%initialize storage
T_cycle = zeros(size(slack_ind,1)-1,size(s,2));
u_cycle = zeros(size(T_cycle));
davg_cycle = zeros(size(T_cycle));
dstd_cycle = zeros(size(T_cycle));

%populate cycle statistics
for j = 1:size(s,2) 
   
    for i = 1:size(slack_ind,1)-1
        
        T_cycle(i,j) = (slack_t(i+1,j)-slack_t(i,j))*24;     %cycle duration (hr)
        
        %peak velocity during ebb or flood tide
        if nanmean(s(slack_ind(i,j):slack_ind(i+1,j),j))>0
            u_cycle(i,j) = nanmax(s(slack_ind(i,j):slack_ind(i+1,j),j));
        else
            u_cycle(i,j) = nanmin(s(slack_ind(i,j):slack_ind(i+1,j),j));
        end
        
        %mean direction and direction deviation
        davg_cycle(i,j) = nanmean(d(slack_ind(i,j):slack_ind(i+1,j),j));
        dstd_cycle(i,j) = nanstd(d(slack_ind(i,j):slack_ind(i+1,j),j));
        
    end
end

%store underlying values
cycle.T_cycle = T_cycle;
cycle.u_cycle = u_cycle;
cycle.davg_cycle = davg_cycle;
cycle.dstd_cycle = dstd_cycle;

end

