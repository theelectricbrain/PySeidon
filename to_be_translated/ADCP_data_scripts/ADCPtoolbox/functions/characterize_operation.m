%Brian Polagye
%August 1, 2011

%Description: routine to characterize turbine operational parameters

%Change Log:
% July 23, 2012 -- JMM: Included u and v (at hub  height in the outputs)

%Inputs:
%   - stat_file: file containing processed statistics (output from
%       characterize_ADCP.m)
%   - stat_path: path for stat_file
%   - turbine: structure containing turbine properties (diameter,
%       efficiency, etc.)
%   - prop: seawater properties
%   - out_path: path to store operational characterization

%Outputs:
%   - generation: structure containing statistics about power generation
%   - operation: structure containing statistics about turbine operation
%   - s_bins: horizontal speed bins for probability distributions
%   - v_bins: horizontal velocity bins for probability distributions
%   - d_bins: horizontal velocity direction bins for probability
%       distributions
%   - z: vertical coordinate for analysis relative to seabed (m)
%   - bin: vertical bin for analysis
%   - turbine: structure containing turbine properties
%   - s: underlying horizontal velocity series
%   - u: underlying eastward velocity series
%   - v: underlying northward velocity series
%   - d: unerlying direction series
%   - P: time series of turbine power ouput
%   - t: underlying time series accompanying 's'
%   - pdf_v: probability distribution of horizontal velocity
%   - pdf_s: probability distribution of horizontal speed
%   - jpdf_vd: joint probability distribution of horizontal velocity and
%       direction

function exit_flag = characterize_operation(stat_file, stat_path,turbine,prop, out_path)

%% Load site data

disp(['Loading ' stat_path stat_file '...'])
load([stat_path stat_file])


% Create output file name
out_file = [stat_file(1:end-4) '_TurbineID' num2str(turbine.ID)]; 
%disp(out_file)

%% Extract hub height information

%find bin closest to specified hub height
if turbine.hub > 1
    bin = find(abs(z-turbine.hub)==min(abs(z-turbine.hub)),1,'first');
else
    bin = find(abs(z-turbine.hub*all.h)==min(abs(z-turbine.hub*all.h)),1,'first');
end

%extract parameters
z = z(bin);
s = ens_s(:,bin);
u = ens_u(:,bin);
v = ens_v(:,bin);
t = ens_t;
d = ens_d(:,bin);
jpdf_vd = jpdf_vd(:,:,bin);
pdf_v = pdf_v(:,bin);

%cleanup input variables not used in operational calculations
clear ens_* slack_t res_s cycle ellipse bins

%% Determine turbine parameters

%calculate rated power
turbine.A = pi*turbine.D^2/4;
turbine.rated_p = 1/2*prop.rho*turbine.rated_u.^3*turbine.A*turbine.eta/1e6;    %rated power (MW)

%fixed, optimized yaw angle
if turbine.yaw == -999
   turbine.angle = optimize_angle(jpdf_vd, v_bins, d_bins, turbine);
%fixed, arbitrary yaw angle
else
   turbine.angle = turbine.yaw;
end

%% Power generation

%general statistics based on joint probability distribution
[generation.power_avg,generation.cf,generation.power_max, ...
    generation.energy_vpdf, generation.energy_spdf, pdf_s, s_bins] = ...
    power_generation(jpdf_vd, pdf_v, turbine, v_bins, d_bins, prop);

%power generation time series
P = power_generation_series(s,d,t,prop,turbine);

%% Operation timing
[operation.windows, operation.f_op, operation.f_on, operation.f_off, operation.T_bins] ...
    = operation_time(s, t, turbine.cutin_u);

%% Save results
savefile = [out_path out_file '.mat'];
disp(['Output saved to: ', savefile])
save(savefile,'generation','operation','*_bins','z','bin','turbine',...
    's','u','v','t','d','P','pdf_*','jpdf_*')

%set exit flag
exit_flag = 1;

end

%% Yaw angle optimization
function [angle] = optimize_angle(jpdf_vd, v_bins, d_bins, turbine)
disp('Optimizing yaw angle')
%determine optimal deployment angle to maximize power generation

%for each bin in distribution, calclate power generated for a passive yaw system
jpdf_power = zeros(size(jpdf_vd));
for i = 1:size(jpdf_vd,1)       %speed
    for j = 1:size(jpdf_vd,2)   %direction
        
        %between cut-in and rated speed
        if abs(v_bins(i))>=turbine.cutin_u && abs(v_bins(i))<turbine.rated_u
            jpdf_power(i,j) = abs(v_bins(i)).^3;
        %greater than rated speed
        elseif abs(v_bins(i))>=turbine.rated_u
            jpdf_power(i,j) = (turbine.rated_u).^3;
        %below cut-in speed
        else
           jpdf_power(i,j) = 0; 
        end
    end
    
end

power_avg_max = sum(sum(jpdf_vd.*jpdf_power));  %average power generation

%set minimization options
options = optimset('MaxFunEvals',1e6,'MaxIter',1e6, 'TolFun', 1e-6, 'TolX', 1e-6);

%find optimal turbine direction - angle coming closest to a free yaw system
[out, fval, exitflag, output]  = ...
    fminsearch('power_angle',[0],options,jpdf_vd, v_bins, d_bins, turbine, power_avg_max);

angle = out(1);

end

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

%% Power generation time series
function [P] = power_generation_series(s,d,t,prop,turbine)

%Inputs:
%   - s: horizontal current velocity time series
%   - d: horizontal current direction time series
%   - t: time stamps
%   - prop: seawater properties
%   - turbine: turbine properties

%Outputs:
%   - P: power generation (MW)

%Initialize output
P = zeros(size(s));

%calculate effective velocity
if turbine.yaw ~=0
   d_theta = (d-turbine.angle)*pi/180;
else
   d_theta = zeros(size(d)); 
end
v_eff = abs(s.*cos(d_theta));

%power at greater than rated speed
pts = find(v_eff>turbine.rated_u);
P(pts) = turbine.rated_p;

%power between rated and cut-in speed
pts = find(v_eff>=turbine.cutin_u & v_eff<turbine.rated_u);
P(pts) = 1/2*prop.rho*v_eff(pts).^3*turbine.A*turbine.eta/1e6;

end

%% Operating time metrics
function [windows, f_op, f_on, f_off, T_bins] = operation_time(s, t, cutin)

%Inputs:
%   -s: horizontal current velocity time series
%   -t: time stamps accompanying s
%   -cutin: turbine cut-in speed

%Outputs:
%   -windows: duration of operating state (:,1), operating state (:,2) (1=operating, 0=idle)
%   -f_op: fraction of time operating
%   -f_on: probability distribution of operating durations
%   -f_off: probability distribution of idle durations
%   -T_bins: durations for probability distributions

%determine state of turbine at start of time series
if abs(s(1))<cutin
    ind_off = 1;
    current_state = 0;
else
    ind_on = 1;
    current_state = 1;
end

%loop through all operational and idle periods
i = 1;
while 1
    
    %if turbine is operating
    if current_state == 1
        %find index turbine cuts out
        ind_off = ind_on + find(abs(s(ind_on:end))<cutin,1,'first')-1;
        current_state = 0;
        if isempty(ind_off), break, end
        
        T(i,:) = [(t(ind_off)-t(ind_on))*24,1];     %operating period (h)
        i = i + 1;
        
    %if turbine is idle
    else
        %find index turbine cuts back in
        ind_on = ind_off+find(abs(s(ind_off:end))>=cutin,1,'first')-1;
        current_state = 1;
        if isempty(ind_on), break, end
        
        T(i,:) = [(t(ind_on)-t(ind_off))*24,0];     %idle period (h)
        i = i + 1;
    end
end

%generate histograms of operational and idle durations (h)
T_bins = floor(min(T(:,1))):1:ceil(max(T(:,1)));

[n_off] = hist(T(T(:,2)==0),T_bins);
[n_on] = hist(T(T(:,2)==1),T_bins);

windows = T;
f_on = n_on/sum(n_on);
f_off = n_off/sum(n_off);

%determine fraction of time operating
ind = find(abs(s)>=cutin);
f_op = length(ind)/length(s);

end
