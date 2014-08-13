%Brian Polagye
%November 11, 2011

%% Speed: vertical profile power law
function [profile] = fitpower_profile(z, s, min_s)

%Description: for each profile, calculate a best fit power law of the form
%   u = u0*(z/z0)^(1/alpha)
%   Note: this function is time consuming to execute

%Inputs:
%   - z: depth bins
%   - mins_s: minimum velocity to calculate profile - around slack water
%       reversing tides will not be described by a power law
%   - s: velocity time series

%Outputs:
%   profile: depth_average speed (s), power law exponent (exp), and
%       goodness of fit (R2) for each profile

debug = 0;

disp('Fitting vertical profile to power law...')

t_timer = now();

%initialize output variables
exp = NaN(length(s),1);
R2 = NaN(length(s),1);
s_avg = NaN(length(s),1);

%append zero velocity at seabed
u_data = [zeros(size(s,1),1), s];
z_data = [0, z(1:size(s,2))];

%fit each vertical profile to a power law
options = optimset('MaxFunEvals',1e6,'MaxIter',1e3, 'TolFun', 1e-6, 'TolX', 1e-6, 'Display','none');

if debug, figure, end

for i = 1:size(s,1)   %process each profile blocks
    
    if abs(nanmean(s(i,:),2))>=min_s
        
        [out,~,exitflag,~]  = fminsearch('fitpowerlaw',...
            [u_data(i,fix(size(u_data,2)/2))/z_data(fix(length(z_data/2)))^(1/5); 5],...
            options, u_data(i,:), z_data);
        
        %compute R^2 statistic for power law fit
        a = out(1);
        b = out(2);
        R2(i) = gfit2(u_data(i,2:end),a*(z_data(2:end)).^b,'8');
        s_avg(i) = mean(s(i,:),2);
        
        %store power law exponent
        if exitflag == 1
            exp(i) = 1/b;
        else
            exp(i) = NaN;
        end
        
        if debug
            clf
            hold on
            plot(u_data(i,:),z_data,'ob')
            plot(a*(z_data).^b,z_data,'--r')
        end
        
    end
    
end
disp(['Elapsed time: ' num2str(round((now()-t_timer)*24*60*10)/10) ' minutes'])

%store individual profile information (speed, exponent, goodness of fit)
profile.s = s_avg;
profile.exp = exp;
profile.R2 = R2;

end
