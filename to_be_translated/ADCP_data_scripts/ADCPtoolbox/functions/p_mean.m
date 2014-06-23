%Brian Polagye
%November 11, 2011

%Description: calculate mean kinetic power density

function [P_mean_fld, P_mean_ebb, P_mean_all] = p_mean(min_s, rho, s)

%Inputs:
%   - min_s: minimum speed to include in calculations (m/s)
%   - rho: seawater density (kg/m^3)
%   - s: horizontal speed (m/s)

%Outputs:
%   - P_mean_fld: mean power density on flood (kW/m^2)
%   - P_mean_ebb: mean power density on ebb (kW/m^2)
%   - P_mean_all: mean power density over all tides (kW/m^2)

P = zeros(size(s));
P_mean_all = zeros(size(s,2),1);
P_mean_ebb = zeros(size(s,2),1);
P_mean_fld = zeros(size(s,2),1);

for i = 1:size(s,2)
   P(:,i) = abs(1/2*rho*(s(:,i)).^3/1000);      %kinetic power density kW/m2
   
   P_mean_all(i) = nanmean(P(:,i));             %mean power density over all tidal stages
   P_mean_ebb(i) = mean(P(s(:,i)<=-min_s,i));   %mean power density on ebb tide
   P_mean_fld(i) = mean(P(s(:,i)>=min_s,i));    %mean power density on flood tide
end

end