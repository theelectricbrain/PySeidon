function save_FlowFile_BPFormat(fileinfo,adcp,rbr,params,options)
% This function takes a processed ADCP file and saves the flow (velocities,
% pressures, bins, times, etc) as a .mat file. 
% The format of the output is such that the file can be used as input for
% the Polagye toolbox. 
% Inputs:
%   - adcp : A structure containing the raw adcpdata
%   - rbr: A structure containing the rbr data (use [] if no rbr data)
%   - params: A structure containing the parameters (time values, depths,
%   etc)
%   - options: A sructure containing optional input
%
% Justine McMillan
% Sept 29, 2012
%
% Change log
% Jan 3, 2013  - Average rbr data before interpolating
% Aug 7, 2013  - Cut out rbr part of the code
% Dec 14, 2013 - Now also computes along and cross channel flows
% Jan 13, 2014 - Added part to do extra rotation when compass wasn't
%                   calibrated
% Apr 15, 2014 - Set default "depth" measurement as adcp.pressure instead
%                of adcp.depth because adcp.depth is quantized

%% Comments
Comments{1} ='data is in Polagye Tools format';
Comments{2} = 'data.east_vel and data.north_vel are relative to true north';
Comments{3} = ['The parameters were set by ',fileinfo.paramfile];
%% Time
day1 = datevec(adcp.mtime(1));
yd = adcp.mtime - datenum(day1(1),0,0);
tind = find(yd>params.tmin & yd<params.tmax);
time.mtime = adcp.mtime(tind);
dt = mean(diff(time.mtime));


%% Depth
if isempty(rbr)
    disp('Depth as measured by ADCP')
    Comments{end+1} = 'Depth as measured by ADCP';
    if isfield(adcp,'pressure')
        pres.surf = adcp.pressure(tind)+params.dabPS;
        disp('Using pressure from ADCP, should check to see if depth is better')
        Comments{end+1} = 'pres.surf from adcp.pressure';
    else
        pres.surf = adcp.depth(tind)+params.dabPS;
        Comments{end+1} = 'pres.surf from adcp.depth (may be quantized)';
    end
        if sum(pres.surf)<1
        pres.surf = params.approxdepth+pres.surf;
    end
else
    Comments{end+1} = 'Depth as measured by RBR sensor';
    error('NEED TO CONFORM WHETHER RBR PART OF THE CODE WORKS')
    disp('Ensemble averaging rbr data')
    % Calculate ensemble average
    nens = round(dt/(rbr.mtime(2) - rbr.mtime(1)));
    mtimeens = rbr.mtime(nens/2):dt:rbr.mtime(end-nens/2);
    mtimeens = mtimeens+params.rbr_hr_offset/24;
    depthens = calc_ensemble(rbr.depth,nens,1);

    disp('Interpolating the ensembled data')
    pres.surf = interp1(mtimeens,depthens,time.mtime,'linear')+params.dabPS;
    
    if isfield(options,'showRBRavg')
        if options.showRBRavg == 1
            figure
            hold all
            plot(rbr.mtime,rbr.depth+params.dabPS)
            plot(time.mtime,pres.surf,'r','linewidth',2)
            legend('Raw data', 'Averaged data')
        end
    end
end

%% zlevels
z = adcp.config.ranges + params.dabADCP;
zind = find(z>params.zmin & z<params.zmax);
data.bins = z(zind);

%% Currents
data.vert_vel = adcp.vert_vel(zind,tind)';
data.error_vel = adcp.error_vel(zind,tind)';

% If compass wasn't calibrated
if isfield(params,'hdgmod')
    [adcp.east_vel,adcp.north_vel]=rotate_coords(adcp.east_vel,adcp.north_vel,params.hdgmod);
    Comments{end+1} = 'East and north velocity rotated by params.hdgmod';
end

% Rotate east_vel and north_vel to be relative to true north
[data.east_vel,data.north_vel] = rotate_to_true(adcp.east_vel(zind,tind)',...
    adcp.north_vel(zind,tind)',params.declination);

% Direction
data.dir_vel = get_DirFromN(data.east_vel,data.north_vel);

% Signed Speed
spd_all = sqrt(data.east_vel.^2+data.north_vel.^2);
% Determine flood and ebb based on principal direction (Polagye Routine)
disp('Getting signed speed (Principal Direction Method) -- used all speeds')
[s_signed_all, PA_all] = sign_speed(data.east_vel, data.north_vel, spd_all, data.dir_vel, params.flooddir);
data.mag_signed_vel = s_signed_all;

if isfield(options,'showPA')
    if options.showPA == 1
        figure,clf
        plot(PA_all, data.bins)
        hold on
        plot([PA_all(1) PA_all(end)],[mean(pres.surf) mean(pres.surf)],'k')
        xlabel('principal axis direction (clockwise from north?)')
        ylabel('z (m)')
        legend('PA','mean depth')
    end
end
if isfield(options,'checkSS')
    if options.checkSS == 1
        disp('Checking the signed speed at middepth')
        data
        [~,midind] = min(abs(data.bins-params.approxdepth/2));
        dt=mean(diff(time.mtime));
        nhour = floor(1/(24*(dt)))+1;
        
        figure,clf
        hold on
        x = 0;
        y = 0;
        max_east = max(data.east_vel(:,midind));
        max_north = max(data.north_vel(:,midind));
        max_spd = sqrt(max_east.^2+max_north.^2);
        xlim(max_east*[-1 1])
        ylim(max_north*[-1 1])
        xlabel('east vel (m/s)')
        ylabel('north vel (m/s)')
        title('First tidal cycle')
        plot(sin(params.flooddir*pi/180)*max_spd*[1 0],...
            cos(params.flooddir*pi/180)*max_spd*[1 0],'k') 
        legend('flood dir')
        
        axis('equal')
         for ii=1:nhour:nhour*15 %For 15hours
            vvec = [data.east_vel(ii,midind),data.north_vel(ii,midind)];
            
            if data.mag_signed_vel(ii,midind)<0
                Col='b';
            else
                Col = 'r';
            end
            warning off
            arrow([x, y], [x, y]+vvec, 'FaceColor',Col,'EdgeColor',Col);
            warning on
            title(['time = ',num2str((time.mtime(ii) - time.mtime(1))*24),' hours'])
        pause(0.5)
         end
 
    end
end

%% save 
lon = params.lon;
lat = params.lat;
outfile = [fileinfo.outdir fileinfo.flowfile];
disp(['Saving data to ',outfile])
save(outfile,'data','pres','time','lon','lat','params','Comments')

%% Save metadata
metadata.progname=[mfilename('fullpath')];
metadata.date = datestr(now);
metadata.paramfile = fileinfo.paramfile;
save(outfile,'metadata','-append')







