%% Yaw angle optimization
function [angle] = optimize_angle(jpdf_vd, v_bins, d_bins, turbine)
%disp('Optimizing yaw angle')
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

