%Brian Polagye
%November 11, 2011

function [profile] = fitlog_profile(z, s, rho, min_s, z_start, z_end)

%Description: for each profile, calculate a best fit for velocities of the form:
%   u(z) = ushear/K*ln(z/z0) 
%Fitting is carried out for all possible ranges of bins to identify the
%   height of the log layer starting from the first z_start bins 

%   Note: this function is time consuming to execute (hours for full
%   profile)

%Inputs:
%   - z: depth bins
%   - z_start: first bin to attempt log fit
%   - z_end: last bin to attempt log fit (if blank, defaults to full
%   profile)
%   - mins_s: minimum velocity to calculate profile - around slack water
%       reversing tides will not be described by a power law
%   - s: velocity time series

%Outputs:
%   profile: depth_average speed (s), bottom roughness (z0), shear 
%       velocity (ushear), shear stress (tau), log layer thickness (z_log)
%       and goodness of fit (R2) for each profile

debug = 0;  %0 to supress graphical depiction of each fit

K = 0.4;    %Von Karman's constant

disp('Fitting vertical profile to logarithmic law...')

t_timer = now();

u_data = s;
z_data = z(1:size(u_data,2));

if isempty(z_end), z_end = length(z_data); end

%initialize output variables
s_avg = NaN(length(s),1);       %depth averaged velocity

z0 = NaN(length(s),length(z_start:z_end)+1);   %bottom roughness
ushear = NaN(size(z0));         %shear velocity
R2 = NaN(size(z0));             %coefficient of determination for fit
tau = NaN(size(z0));            %shear stress
z_log = NaN(size(z0));          %log layer height
exitflag = zeros(size(z0));     %exit flag from fitting routine

%fit each vertical profile to a logarithmic profile
options = optimset('MaxFunEvals',1e6,'MaxIter',1e3, 'TolFun', 1e-6, 'TolX', 1e-6, 'Display','none');

if debug
    figure
    set(gcf,'position',[62 378 1267 420])
    hold on
end

for i = 1:size(s,1)   %process each velocity profile
    
    s_avg(i) = nanmean(s(i,:),2);  %store depth average velocity
    
    %exlcuding profiles around slack water will speed up the algorithm
    %since these generally have low-quality fits
    if abs(nanmean(s(i,:),2))>=min_s
        
        %calculate log law fit for all specified depth ranges
        for j = z_start:z_end
            index = j + 2 - z_start;    %storage index
            
            %calculate log law fit for first j bins
            [out,~, exitflag(i,index),~] = ...
                fminsearch('fitloglaw',[0.05; 7e-5],options, u_data(i,1:j), z_data(1:j));
            
            ushear(i,index) = out(1);     %shear velocity (m/s)
            z0(i,index) = out(2);         %bottom roughness (m)
            z_log(i,index) = z_data(j);   %assumed log law height (m)
            
            %compute R^2 statistic for log law fit
            if z0(i,index) <= 0
                R2(i,index) = 0;  %exclude non-physical results for bottom roughness 
            else
                R2(i,index) = gfit2(u_data(i,1:j),ushear(i,index)/K*log(z_data(1:j)/z0(i,index)),'8'); 
            end
        
            if debug && z0(i,index)>0, plot(ushear(i,index)/K*log(z_data(1:j)/z0(i,index)),z_data(1:j),'--r'), end          
        end
        
        %find index corresponding to maximum R^2 value on a valid exit
        max_R2 = find(R2(i,exitflag(i,2:end)>0)==max(R2(i,exitflag(i,2:end)>0)),1,'last');
            
        if ~isempty(max_R2)
            R2(i,1) = R2(i,max_R2);
            z0(i,1) = z0(i,max_R2);
            ushear(i,1) = ushear(i,max_R2);
            tau(i,:) = ushear(i,:).^2*rho;
            z_log(i,1) = z_data(max_R2);
            exitflag(i,1) = exitflag(i,max_R2);
            
            if debug
                plot(u_data(i,:),z_data,'ob')
                hold on
                
                if R2(i,1)>0.95
                    plot(ushear(i,1)/K*log(z_data(1:max_R2-1+z_start)/z0(i,1)),z_data(1:max_R2-1+z_start),'-g','linewidth',3)
                end
                set(gca,'XLim',[-3 3])
                text(-2.5, 25, ['R^2 = ' num2str(round(R2(i,1)*1e3)/1e3)]);
                text(-2.5, 15, ['\tau = ' num2str(round(tau(i,1)*100)/100) ...
                    ' (' num2str(round(nanmin(tau(i,2:end))*100)/100) '-' num2str(round(nanmax(tau(i,2:end))*100)/100) ')'])
                pause(0.2)
                
                clf
                hold on
            end
        end
    end
    
end
disp(['Elapsed time: ' num2str(round((now()-t_timer)*24*60*10)/10) ' minutes'])

%store individual profile information (speed, roughness, shear, goodness of fit)
profile.s = s_avg;
profile.z0 = z0;
profile.ushear = ushear;
profile.R2 = R2;
profile.tau = tau;
profile.z_log = z_log;
profile.exitflag = exitflag;

end
