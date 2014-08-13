function [ens_t, ens_s, ens_d, ens_u, ens_v, ens_w, ens_h, exit_flag] = ensemble_data(t_ens, dt, t, s, d, u, v, w, h)

%Description: ensemble average (time and space) measurements to smooth effect of turbulence

%Inputs:
%   - t_ens: ensemble time in seconds
%   - dt: time interval in seconds
%   - t: time stamps (matlab format)
%   - s: horizontal velocity (signed ebb and flood)
%   - d: horozintal velocity direction
%   - u: east velocity
%   - v: north velocity
%   - w: vertical velocity
%   - h: free surface

%Outputs: ensemble average inputs, time stamps for start of ensemble

exit_flag = 0;  %default to successful execution
ens_tol = 0.1;  %maximum allowable mismatch between number of integer samples in ensemble

ens = t_ens/dt;   %number of points in ensemble average

if abs(round(ens)-ens) > ens_tol || ens+ens_tol < 1
    exit_flag = 1;
    disp('Error: mismatch between measurement time step and specified ensemble averaging time greater than tolerance')
    ens_t = NaN;
    ens_s = NaN;
    ens_d = NaN;
    ens_u = NaN;
    ens_v = NaN;
    ens_w = NaN;
    ens_h = NaN;
else
    ens = round(ens);
    disp(['Number of points in ensemble average: ',num2str(ens)])
    
    if ens == 1
        disp('Not computing ensembles')
        ens_t = t;
        ens_s = s;
        ens_d = d;
        ens_u = u;
        ens_v = v;
        ens_w = w;
        ens_h = h;
    else
        ens_t = t(1):t_ens/(60*60*24):t(end);    %time stamp is start of averaging period
       
        %temp = datevec(ens_t);
        %ens_t = datenum(temp(:,1),temp(:,2),temp(:,3),temp(:,4),round(temp(:,5)+temp(:,6)/60),0);   %clean up time stamps to integer minutes
        %ens_t = datenum(temp(:,1),temp(:,2),temp(:,3),temp(:,4),temp(:,5),round(temp(:,6)));   %clean up time stamps to integer seconds
     
        ens_s = calc_ensemble(s,ens,1);     %horizontal velocity
        ens_d = calc_ensemble(d,ens,1);     %direction
        ens_u = calc_ensemble(u,ens,1);     %east velocity
        ens_v = calc_ensemble(v,ens,1);     %west velocity
        ens_w = calc_ensemble(w,ens,1);     %vertical velocity
        if length(h)>1
            ens_h = calc_ensemble(h,ens,2);     %free surface elevation
        else
            ens_h = h;
        end
        
        %prune time series for incomplete ensembles
        if length(ens_t) > size(ens_s,1)
            ens_t = ens_t(1:size(ens_s,1));
        end
    end
end

end
