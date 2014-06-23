function [time_out,data_out,pres_out] = EnsembleData_FlowFile(t_ens, time, data,pres,options)

%Description: ensemble average (time and space) measurements that are in
% BrianPolagye's format
% Modification of Brian Polagye's ensemble_data.m code
% Here the time vector is first set up and then looped through (this allows
% for nonuniform delta t values which tends to happen when using high
% sampling rates). The standard deviation is also computed and outputted
%
%Inputs:
%   - t_ens: ensemble time in seconds
%   - time, data, pres: fields of Flow Files
%   - options: structure with some options (i.e. which method to use -- TimeSearch or NumPoints)
%
% Outputs:
%    - time_out, data_out, pres_out: same format as input files (with std
%    dev included)
%
% Justine McMillan
% Oct 28, 2013
%
% Changes:
% - Dec 14, 2013: Now computes Ualong and Ucross for each ensemble where
%                 direction is determined based on the mean direction of 
%                 the flow during that ensemble
% - Apr 21, 2014: Modified a few parts so that it would compute averages of
%                  bursts correctly (See code used for DG140115BPb ADCP

%%
if ~isfield(options,'method'), options.method = 'TimeSearch'; end
if ~isfield(options,'showoutput'), options.showoutput = 1; end

if strcmp(options.method,'TimeSearch') == 1
    if options.showoutput
    disp('[COMPUTING ENSEMBLES] - Method: Time Search')
    end
    ens_tol = 0.1;  % allow for 10% of the points to be missing, and still compute the ensemble
    
    if round(t_ens/60)>60
        warning(['Check to make sure ensemble averaging is done correctly.',...
            'Originally created for much shorter ensembles'])
    end
    
    %% Set up time vector -- want timestamps to be on even minutes
    SEC2DAY = 1/(24*3600);
    dt = mean(diff(time.mtime))/SEC2DAY; % delta t in seconds
    
    dt_one = (time.mtime(2)-time.mtime(1))*24*3600;
    
    [yr,mm,dd,hr,mins,sec] = datevec(time.mtime(1));
    
    time_out.mtime = [datenum(yr,mm,dd,hr,0,0):t_ens/3600/24:time.mtime(end)]; %start on the hour -- leading values will be removed at the end
    
    
    if options.showoutput
        if find(diff(time.mtime))/SEC2DAY>dt_one*2
            disp(['[COMPUTING ENSEMBLES] - Burst mode?'])
            nensapprox = (t_ens/dt_one);
        else
            
            nensapprox = (t_ens/dt);
        end
   % disp(['[COMPUTING ENSEMBLES] - first delta t: ',num2str(dt_one), ' seconds']);
    disp(['[COMPUTING ENSEMBLES] - Ensemble interval: ',num2str(t_ens/60),' mins',...
        ' (Approx ',num2str(round(nensapprox)),' points)'])
    end
    %% Initialize
    pres_out.surf = NaN*ones(size(time_out.mtime));
    pres_out.surf_std = NaN*ones(size(time_out.mtime));
    fields = {'mag_signed_vel','east_vel','north_vel','vert_vel','error_vel','dir_vel','Ualong','Ucross'};
    nbins = size(data.east_vel,2);
    for ii = 1:length(fields)
        fieldvar =  fields{ii};
        fieldstd = [fields{ii},'_std'];
        data_out.(fieldvar) = NaN*ones(length(time_out.mtime),nbins);
        data_out.(fieldstd) = NaN*ones(length(time_out.mtime),nbins);
    end
    
    
    %% Compute means and standard deviations
    for ii = 1:length(time_out.mtime)
        tstart = time_out.mtime(ii) - 1/2*(t_ens*SEC2DAY);
        tend   = time_out.mtime(ii) + 1/2*(t_ens*SEC2DAY);
        
        
        ind_t = find(time.mtime>=tstart & time.mtime<tend);
%         if ii<20
%             datestr(tstart)
%             datestr(tend)
%             disp(['numpoints = ',num2str(length(ind_t))])
%             pause
%         end
        
        if length(ind_t)<nensapprox*(1-ens_tol)
            time_out.mtime(ii) = NaN;
            pres_out.surf(ii)  = NaN;
            pres_out.surf_std(ii) = NaN;
        else
            pres_out.surf(ii) = nanmean(pres.surf(ind_t));
            pres_out.surf_std(ii) = nanstd(pres.surf(ind_t));
        end
        % loop through data fields
        for jj = 1:length(fields)
            fieldvar =  fields{jj};
            fieldstd = [fields{jj},'_std'];
            if length(ind_t)<nensapprox*(1-ens_tol)
                data_out.(fieldvar)(ii,:) = NaN;
                data_out.(fieldstd)(ii,:) = NaN;
            else
                if strcmp(fields{jj},'Ucross') == 0 && strcmp(fields{jj},'Ualong') == 0
                    eval(['field = data.',fields{jj},';']);
                    %Compute mean
                    if strcmp(fields{jj},'dir_vel') %Treat direction separately
                        data_out.dir_vel(ii,:) = get_DirFromN(data_out.east_vel(ii,:),data_out.north_vel(ii,:));
                    else
                        data_out.(fieldvar)(ii,:)    = nanmean(field(ind_t,:));
                    end
                    %Compute std dev
                    data_out.(fieldstd)(ii,:) =  nanstd(field(ind_t,:));
                end
            end %if
        end %jj
        
        % Ualong and Ucross - use average angle for the ensemble
        if data_out.mag_signed_vel(ii,:)>0 %Flood
            theta = -data_out.dir_vel(ii,:);
        else %ebb
            theta = 180-data_out.dir_vel(ii,:);
        end
        thetamat = ones(length(ind_t),1)*theta;
        [Ucross, Ualong] = rotate_to_channelcoords(data.east_vel(ind_t,:),...
                data.north_vel(ind_t,:),thetamat);
        
        data_out.Ualong(ii,:) = nanmean(Ualong);
        data_out.Ualong_std(ii,:) = nanstd(Ualong);
        data_out.Ucross(ii,:) = nanmean(Ucross);
        data_out.Ucross_std(ii,:) = nanstd(Ucross);
    
        clear ind_t Ucross Ualong
    end
    
    %% Display some info to the screen
    
    ind_t2 = find(~isnan(time_out.mtime));

    if options.showoutput
        t_drop1 = (time_out.mtime(ind_t2(1)) - time.mtime(1) - 1/2*(t_ens*SEC2DAY))*24*60;
        t_drop2 = (time.mtime(end) - time_out.mtime(ind_t2(end)) - 1/2*(t_ens*SEC2DAY))*24*60;
    
        disp(['[COMPUTING ENSEMBLES] - Dropped ',num2str(t_drop1),' mins of data at the start of file'])
        disp(['[COMPUTING ENSEMBLES] - Dropped ',num2str(t_drop2),' mins of data at the end of file'])
    end
    %% Trim Nans from the front and end
    time_out.mtime = time_out.mtime(ind_t2(1):ind_t2(end));
    pres_out.surf = pres_out.surf(ind_t2(1):ind_t2(end));
    pres_out.surf_std = pres_out.surf_std(ind_t2(1):ind_t2(end));
    for ii = 1:length(fields)
        fieldvar =  fields{ii};
        fieldstd = [fields{ii},'_std'];
        data_out.(fieldvar) = data_out.(fieldvar)(ind_t2(1):ind_t2(end),:);
        data_out.(fieldstd) = data_out.(fieldstd)(ind_t2(1):ind_t2(end),:);
    end
    
elseif strcmp(options.method,'NumPoints')
    error('Method not set up yet -- just copy frome ensemble_data.m?')
else
    error(['Invalid Method'])
end

