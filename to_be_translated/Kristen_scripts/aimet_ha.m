function ha = aimet_ha(ts)
% Perform harmonic analysis using t_tide on a time series, ts.
% This function assumes that the first dimension is the time series, with
% additional dimenions representing space.
% Input:
% If ts is a vector, it is a time series from one point in space
% If ts has two dimensions, [txi], there are time series for i points.
% If ts has three dimensions, [txjxi], there are time series for j*i points
%   in space (j in y and i in x).
% If ts has four dimensions, [txkxjxi], there are time series for k*j*i
%   points in space (k in z, j in y, and i in x).
% Note that if the time series numbers can be real or complex, and complex
% should be used for velocity.
%
% Output:
%  ha  Structure in t_tide format containing fields: 'names', 'freq'
%      and 'tidecon'

%% File series structure identifying input files available
% disp('temp change to dt')
% % SD = roms_createSeriesDef('~/Desktop/ai','ocean_his_');
% % SD = roms_createSeriesDef('~/roms/projects/ai100','ocean_his_');
SD = roms_createSeriesDef('~/roms/projects/ai65/OUT','ocean_his_');
int = SD.nctime(2)-SD.nctime(1); %in days
int = int*24; %in hours
% int =  0.25;%0.2667;
% disp('int has been set directly')

%% Harmonic Analysis  
N = ndims(ts);

if N == 1
    if ~isnan(nansum(ts)) % make sure not on land
        [ha,~]=t_tide(ts,'interval',int,'latitude',48,...
            'start time',SD.nctime(1),'output','none');
    end
elseif N == 2
    for i = 1:length(ts(1,:))
        if ~isnan(nansum(ts(:,i))) % make sure not on land
            [ha(i),~]=t_tide(ts(:,i),'interval',int,'latitude',48,...
                'start time',SD.nctime(1),'output','none');
        else
            ha(i).name=[]; ha(i).freq=[]; ha(i).tidecon=[];
        end
    end
elseif N == 3
    for i = 1:length(ts(1,1,:))
        for j = 1:length(ts(1,:,1))
            if sum(~isnan(ts(:,j,i))) > 72 % make sure not on land
                [ha(j,i),~]=t_tide(ts(:,j,i),'interval',int,'latitude',48,...
                    'start time',SD.nctime(1),'output','none');
            else
                ha(j,i).name=[]; ha(j,i).freq=[]; ha(j,i).tidecon=[];
            end
        end
    end
elseif N == 4
    for i = 1:length(ts(1,1,1,:))
        for j = 1:length(ts(1,1,:,1))
            for k= 1:length(ts(1,:,1,1))
                if ~isnan(nansum(ts(:,k,j,i))) % make sure not on land
                    [ha(k,j,i),~]=t_tide(ts(:,k,j,i),'interval',int,'latitude',48,...
                        'start time',SD.nctime(1),'output','none');
                else
                    ha(k,j,i).name=[]; ha(k,j,i).freq=[]; ha(k,j,i).tidecon=[];
                end
            end
        end
    end
end