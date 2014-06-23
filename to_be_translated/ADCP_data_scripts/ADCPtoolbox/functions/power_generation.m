%% Power generation metrics
function [power_avg, cf, power_max, energy_vpdf, energy_spdf, pdf_s, s_bins] = ...
    power_generation(jpdf_vd, pdf_v, turbine, v_bins, d_bins, prop)

%Inputs:
%   - jpdf_vd: joint probability distribution of horizontal velocity with direction
%   - turbine: turbine parameters (structure)
%   - v_bins: horizontal velocity bins (m/s)
%   - d_bins: velocity direction bins (degrees)
%   - prop: seawater properties

%Outputs:
%   - power_avg: average power generated (MW)
%   - cf: capacity factor
%   - power_max: maximum power generated (MW)
%   - energy_vpdf: probability distribution of energy generation as function
%       of horizontal velocity
%   - energy_spdf: probability distribution of energy generation as function
%       of horizontal speed

%calculate off-axis angle for each bin center
if turbine.yaw ~=0
    d_theta = (d_bins-turbine.angle)*pi/180;
else
    d_theta = zeros(size(d_bins));
end

%initialize power
jpdf_power = zeros(size(jpdf_vd));

%For each bin in distribution, calclate power generated - accounting for
%off-axis flow angles. Off-axis flow is treated as a reduction in inflow
%velocity (per. R. Thresher - Aug 2011)
for i = 1:size(jpdf_vd,1)       %speed
    for j = 1:size(jpdf_vd,2)   %direction
        
        v_eff = abs(v_bins(i)*cos(d_theta(j)));
        
        %between cut-in and rated speed
        if v_eff>=turbine.cutin_u && v_eff<turbine.rated_u
            jpdf_power(i,j) = 1/2*prop.rho*v_eff.^3*turbine.A*turbine.eta/1e6;
        
        %greater than rated speed
        elseif v_eff>=turbine.rated_u
            jpdf_power(i,j) = turbine.rated_p;
        
        %below cut-in speed
        else
           jpdf_power(i,j) = 0; 
        end
    end
    
end

%calculate convolution of power at each stage with probability of occurence
jpdf_c = jpdf_vd;
jpdf_c(jpdf_c>0)=1;

power_avg = sum(sum(jpdf_vd.*jpdf_power));      %average power generation (MW)
cf = power_avg./max(max(jpdf_power.*jpdf_c));   %capacity factor - using power occuring at least once
power_max = power_avg./cf;                      %maximum power generation (MW)

%probability distribution of energy generated as a function of horizontal velocity
energy_vpdf = sum(jpdf_vd.*jpdf_power,2)/sum(sum(jpdf_vd.*jpdf_power));

%probability distribution of energy generated as a function of speed
s_bins = unique(abs(v_bins));
energy_spdf = zeros(size(s_bins));
pdf_s = zeros(size(s_bins));

for i = 1:length(s_bins)
   ind = find(abs(v_bins)==s_bins(i));
   energy_spdf(i) = sum(energy_vpdf(ind));
   pdf_s(i) = sum(pdf_v(ind));
end

end