function aimet_plot(x,y,z,color,varargin)%(ijz,x,y,z,color)

grid = roms_get_grid('OUT/ocean_his_0002.nc','OUT/ocean_his_0002.nc');
% grid = roms_get_grid('results/ocean_his_0002.nc','results/ocean_his_0002.nc');

% File series structure identifying input files available
% SD = roms_createSeriesDef('~/Desktop/ai','ocean_his_');
% SD = roms_createSeriesDef('~/roms/projects/ai65/OUT','ocean_his_');
% SD = roms_createSeriesDef('/drive2/troc/Headland_farm/results/','ocean_his_');
SD = roms_createSeriesDef('./OUT','ocean_his_');
font = 20;

fname = 'figure';

if ~isempty(varargin)
    iarrows = strcmp(varargin,{'arrows'});
    if sum(iarrows) 
        u = varargin{find(iarrows)+1}; 
        v = varargin{find(iarrows)+2}; 
        xu = varargin{find(iarrows)+3};
        yu = varargin{find(iarrows)+4};
    end
    ifname = strcmp(varargin,{'fname'});
    if sum(ifname)
        fname = varargin{find(ifname)+1};
    end
    imovie = strcmp(varargin,{'movie'});
    iclims = strcmp(varargin,{'clims'});
    if sum(iclims)
        clims = varargin{find(iclims)+1};
    end
    ixlims = strcmp(varargin,{'xlims'});
    if sum(ixlims)
        xlims = varargin{find(ixlims)+1};
    end
    iylims = strcmp(varargin,{'ylims'});
    if sum(iylims)
        ylims = varargin{find(iylims)+1};
    end
    nofig = strcmp(varargin,{'nofig'});
    ishowsignal = strcmp(varargin,{'showsignal'}); % for movie snapshots when have 1 timestep at a time
    if sum(ishowsignal)
        tind = varargin{find(ishowsignal)+1};
%         load matfiles/signal.mat %roc
        load('matfiles/metsh10byt.mat','zsave'); signal=zsave; %ai65
    end
    % Axis labels
    ixlabel = strcmp(varargin,{'xlabel'});
    if sum(ixlabel)
        xname = varargin(find(ixlabel)+1);
    end
    iylabel = strcmp(varargin,{'ylabel'});
    if sum(iylabel)
        yname = varargin(find(iylabel)+1);
    end
    izlabel = strcmp(varargin,{'zlabel'});
    if sum(izlabel)
        zname = varargin(find(izlabel)+1);
    end
    %%%
    % Contour options
    icontour = strcmp(varargin,{'contour'});
    if sum(icontour)
        conoption = varargin(find(icontour)+1);
    end
else
    iarrows = 0;
    imovie = 0;
    iclims = 0;
    nofig = 0;
    ishowsignal = 0;
    ixlabel = 0;
    iylabel = 0;
    izlabel = 0;
    icontour = 0;
    ixlims = 0;
    iylims = 0;
%     ifname = 0;
end

figure
% % FOR SLICES
% set(gcf,'position',[0         275        1228         409]) % slices
% % subaxis(1,1,1,'ml',.085,'mb',.21,'mr',.075,'mt',.1) % slices
% FOR x-y
set(gcf,'position',[148    31   731   634]) 
% set(gcf,'position',[148    31   731   634]) 
% subaxis(1,1,1,'ml',.09,'mb',.06,'mr',.105,'mt',.06)
subaxis(1,1,1,'ml',.105,'mb',.06,'mr',.105,'mt',.06)

if ~sum(imovie) % Not a movie
%     contourf(x,y,z)
    pcolorjw(x,y,z)
    colormap(color), lcon = size(color,1); d = floor(lcon/2);
    colorbar
    if sum(iclims)
        caxis(clims);
        cmax = clims(2); cmin = clims(1);
    else
        cmin = min(min(z)); cmax = max(max(z)); cmax = max(cmax,abs(cmin));
%         caxis([-cmax cmax])
%         caxis([cmin cmax])
    end
    if sum(ixlims)
        xlim(xlims)
    end
    if sum(iylims)
        ylim(ylims)
    end
    %     cmin = 0;
%     cmax = 50;
    hold on
%     if cmin < 0
%         caxis([-cmax cmax])
    if sum(icontour)
%         set(gca,'clim',[cmin cmax])
        switch conoption{:}
            case 'specialdiv' % when I want contour at specific values for divergent info
                aa=2*cmax/lcon;
                contour(x,y,z,[-(aa/2+aa*(d:-1:0)) (aa/2+aa*(0:1:d))],'k')
            case 'specialseq' % better contours for sequential data
                aa=(cmax-cmin)/lcon;
                contour(x,y,z,[cmin+aa*(0:lcon)],'k')                
%                 contour(x,y,z,[cmin+aa*(0:lcon)/2],'k')                
            case 'normal'
               [cs,h] = contour(x,y,z,'k'); % normal contour plot choice
%                [cs,h] = contour(x,y,z,cmin:(cmax-cmin)/10:cmax,'k'); % normal contour plot choice
        end
    end
%     [cs,h] = contour(x,y,z,[10 20 30 35 40 45],'k');
%     [cs,h] = contour(x,y,z,'k');
%     clabel(cs,h,'manual','color','m','fontsize',13,'fontweight','bold')

%     %%% Pilot Site
%     plot(-122.6877,48.1528,'r.','markersize',20)
%     %%% Contour for MKPD
%     [cs,h] = contour(x,y,z,[14.9 14.9],'k','linewidth',2);
%     clabel(cs,h,'manual','fontsize',15,'fontweight','bold')    
%     %%%
%     %%% Contour for mI
%     [cs,h] = contour(x,y,z,[10 10],'k','linewidth',2);
%     clabel(cs,h,'manual','fontsize',15,'fontweight','bold')    
%     %%%
%     %%% Contour for add
%     [cs,h] = contour(x,y,z,[15 40],'k');
%     clabel(cs,h,'manual','fontsize',15,'fontweight','bold')
%     %%%
%     else
%         caxis([cmin cmax])
%         aa=2*cmax/lcon;
%         contour(x,y,z,[(aa/2+aa*(0:8:d))],'k')
%         [cs,h] = contour(x,y,z,[10 10],'k');
% % %        [cs,h] = contour(x,y,z,'k'); % normal contour plot choice
%         clabel(cs,h,'fontsize',15,'color','k','rotation',0)
%     end
% [cs,h]=contour(x,y,z,[15 40],'k');
% clabel(cs,h,'fontsize',13,'fontweight','bold','color','k','labelspacing',150)
    freezeColors
    contour(grid.lon_rho,grid.lat_rho,grid.mask_rho,[0 0],'k')
    if sum(iarrows) % do arrows
        doarrows(xu,yu,u,v)
%         % dx=600 thalweg
%         % dy=10 thalweg
%         % slices1-3: dx=150, dy=4, slice4: dx=100, dy=4;
%         % slice5,7:dx=125,dy=6; slice6:dx=75,dy=5
%         dx = 600; dy = 10;
%         [xarrow,yarrow,u] = griddata(xu,yu,u,xu(1,1):dx:xu(1,end),...
%             (min(min(yu)):dy:max(max(yu)))');
%         [v] = griddata(xu,yu,v,xu(1,1):dx:xu(1,end),...
%             (min(min(yu)):dy:max(max(yu)))');
%         doarrows(xarrow,yarrow,u,v)
    end
    set(gca,'fontsize',font,'fontweight','bold')
    if sum(ixlabel)
        xlabel(xname)
    end
    if sum (iylabel)
        ylabel(yname)
    end
    if sum(izlabel)
        ylabel(colorbar,zname,'fontsize',font,'fontweight','bold') % colorbar label
    end
    set(gca,'fontsize',font,'fontweight','bold')
%     %%% For turbulence plots
%     plot(-122.6855,48.1515,'b.','markersize',20)
%     %%%
%     %%% Just for vort analysis plots
%     lon1 = -122.7013; lon2 = -122.6676; lat1 = 48.1492; lat2 = 48.1666; %try2
%     plot([lon1 lon1 lon2 lon2 lon1],[lat1 lat2 lat2 lat1 lat1],'k','linewidth',3)
%     lon1 = -122.7013; lon2 = -122.681; lat1 = 48.1634; lat2 = 48.1666; %box1
%     plot([lon1 lon1 lon2 lon2 lon1],[lat1 lat2 lat2 lat1 lat1],'k','linewidth',2)
%     lon1 = -122.7013; lon2 = -122.68; lat1 = 48.1602; lat2 = 48.1634; % box2
%     plot([lon1 lon1 lon2 lon2 lon1],[lat1 lat2 lat2 lat1 lat1],'k','linewidth',2)
%     lon1 = -122.7013; lon2 = -122.678; lat1 = 48.1570; lat2 = 48.1602; % box3
%     plot([lon1 lon1 lon2 lon2 lon1],[lat1 lat2 lat2 lat1 lat1],'k','linewidth',2)
%     lon1 = -122.7013; lon2 = -122.678; lat1 = 48.1538; lat2 = 48.1570; % box4
%     plot([lon1 lon1 lon2 lon2 lon1],[lat1 lat2 lat2 lat1 lat1],'k','linewidth',2)
%     lon1 = -122.7013; lon2 = -122.678; lat1 = 48.1492; lat2 = 48.1538; % box5
%     plot([lon1 lon1 lon2 lon2 lon1],[lat1 lat2 lat2 lat1 lat1],'k','linewidth',2)
%     lon1 = -122.678; lon2 = -122.6676; lat1 = 48.1492; lat2 = 48.1602; % box6
%     plot([lon1 lon1 lon2 lon2 lon1],[lat1 lat2 lat2 lat1 lat1],'k','linewidth',2)
%     %%%
%     %%% Just for velocity anomaly plots
%     plot([-122.6881 -122.688033333333],[48.1525 48.153],'k.','markersize',16)
%     %%%
    if sum(ishowsignal)
%         h = axes('Position', [.7 .21 .175 .19], 'Layer','top'); %roc
%         h = axes('Position', [.09 .21 .175 .16], 'Layer','top'); %sill1
%         h = axes('Position', [.055 .21 .175 .15], 'Layer','top'); %slices
%         h = axes('Position', [.55 .24 .33 .15], 'Layer','top'); %thalweg
        h = axes('Position', [.42 .77 .38 .15], 'Layer','top'); %x-y
%         h = axes('Position', [.52 .77 .28 .15], 'Layer','top'); %x-y zoom
%         plot(SD.nctime(1101:2880),signal(1101:2880),'k','linewidth',1), hold on, 
        plot(SD.nctime(1:1100),signal(1:1100),'k','linewidth',1), hold on, 
        plot(SD.nctime(tind),signal(tind),'r.','markersize',20)
        axis(h, 'off', 'tight')
%         set(gcf,'color','none')
    end
    set(gcf,'color','none')
%     axis equal 
%     axis tight
    fname=['figures/' fname];
    %export_fig(fname,'-pdf','-painters') % too high res and lines in adobe
    if ~sum(nofig)
        saveas(gcf,fname,'fig')
    end
    savefig(fname,'png','-r100')
    savefig(fname,'pdf','-r100')
else % a movie, assume 1st dimension is time
    % Load in surface signal
%     load 'signal.mat' % for ai100, Free surface at point off Admiralty Head
%     load('matfiles/metsh10byt.mat','zsave'); signal=zsave; %ai65
    mkdir(['figures/' fname]);
    if sum(iclims)
        caxis(clims);
        cmax = clims(2); cmin = clims(1);
    else
%         cmin = min(min(z)); cmax = max(max(z)); cmax = max(cmax,abs(cmin));
        cmin = min(min(min(min(z))));
        cmax = max(max(max(max(z))));
%         caxis([cmin cmax])
        if cmin < 0
    %         cmax = max(cmax,abs(cmin));
            clims = [-cmax cmax];
%             clim = [cmin cmax];
        else
    %         clim = [1020 cmax];
            clims = [cmin cmax];
        end
    end
%     clims = [cmin cmax]; %TEMP
%     clim = [0 .5];
    for i=1:size(x,1)
        xt = squeeze(x(i,:,:,:));
        yt = squeeze(y(i,:,:,:));
        zt = squeeze(z(i,:,:,:));
        if sum(iarrows)
            xtu = squeeze(xu(i,:,:,:));
            ytu = squeeze(yu(i,:,:,:));
        end
% subaxis(1,1,1,'ml',.085,'mb',.21,'mr',.075,'mt',.1) % slices
% subaxis(1,1,1,'ml',.05,'mb',.21,'mr',.045,'mt',.1) % slices
% subaxis(1,1,1,'ml',.05,'mb',.15,'mr',.075,'mt',.1) % slices
subaxis(1,1,1,'ml',.085,'mb',.21,'mr',.075,'mt',.1) % slices
        pcolorjw(xt,yt,zt)
        colormap(color), lcon = size(color,1); d = floor(lcon/2);
        colorbar
        caxis(clims)
        if sum(ixlims)
            xlim(xlims)
        end
        if sum(iylims)
            ylim(ylims)
        end
        hold on
    if sum(icontour)
%         set(gca,'clim',[cmin cmax])
        switch conoption{:}
            case 'specialdiv' % when I want contour at specific values for divergent info
                aa=2*cmax/lcon;
                contour(xt,yt,zt,[-(aa/2+aa*(d:-1:0)) (aa/2+aa*(0:1:d))],'k')
            case 'specialseq' % better contours for sequential data
                aa=(cmax-cmin)/lcon;
                contour(xt,yt,zt,[cmin+aa*(0:lcon)],'k')                
%                 contour(x,y,z,[cmin+aa*(0:lcon)/2],'k')                
            case 'normal'
%                [cs,h] = contour(x,y,z,'k'); % normal contour plot choice
               [cs,h] = contour(xt,yt,zt,cmin:(cmax-cmin)/10:cmax,'k'); % normal contour plot choice
        end
    end
        %freezeColors
        contour(grid.lon_rho,grid.lat_rho,grid.mask_rho,[0 0],'k')
        if sum(iarrows) % do arrows
            ut = squeeze(u(i,:,:,:));
            vt = squeeze(v(i,:,:,:));
            doarrows(xtu,ytu,ut,vt)
%             % ai100:
%             % dx=600 thalweg
%             % dy=10 thalweg
%             % slices1-3: dx=150, dy=4, slice4: dx=100, dy=4;
%             % slice5,7:dx=125,dy=6; slice6:dx=75,dy=5
%             % ai65: sill1: dx=200, dy=6; ahead: dx=1; dy=1;
%             dx = 150; dy = 5;
%             [xarrow,yarrow,ut] = griddata(xtu,ytu,ut,xtu(1,1):dx:xtu(1,end),...
%                 (min(min(ytu)):dy:max(max(ytu)))');
%             [vt] = griddata(xtu,ytu,vt,xtu(1,1):dx:xtu(1,end),...
%                 (min(min(ytu)):dy:max(max(ytu)))');
%             doarrows(xarrow,yarrow,ut,vt)
        end
        set(gca,'fontsize',20,'fontweight','bold')
    if sum(ixlabel)
        xlabel(xname)
    end
    if sum (iylabel)
        ylabel(yname)
    end
    if sum(izlabel)
        ylabel(colorbar,zname,'fontsize',font,'fontweight','bold') % colorbar label
    end
        set(gca,'fontsize',20,'fontweight','bold')
    if sum(ishowsignal)
            h = axes('Position', [.755 .23 .1 .19], 'Layer','top'); %Roc
%         h = axes('Position', [.095 .22 .175 .16], 'Layer','top'); %sill1
%         h = axes('Position', [.055 .21 .175 .15], 'Layer','top'); %slices
%         h = axes('Position', [.55 .24 .33 .15], 'Layer','top'); %thalweg
%         h = axes('Position', [.42 .77 .38 .15], 'Layer','top'); %x-y
        plot(SD.nctime(tind),signal,'k','linewidth',2), hold on, 
        plot(SD.nctime(tind(i)),signal(i),'r.','markersize',20)
%         plot(SD.nctime(i),signal(i),'r.','markersize',20)
        axis(h, 'off', 'tight')
        set(gcf,'color','none')
    end
%         % Show tidal signal
%             h = axes('Position', [.7 .21 .175 .19], 'Layer','top'); %Roc
% %         h = axes('Position', [.055 .21 .175 .15], 'Layer','top'); %slices
% %         h = axes('Position', [.55 .24 .33 .15], 'Layer','top'); %thalweg
%         h = axes('Position', [.45 .77 .39 .15], 'Layer','top'); %x-y
%         plot(SD.nctime,signal,'k','linewidth',1), hold on, 
%         plot(SD.nctime(i),signal(i),'r.','markersize',20)
%         axis(h, 'off', 'tight')
%         set(gcf,'color','none')
        %%%%
        fname2=sprintf('figures/%s/%i',fname,i);
        savefig(fname2,'png','-r100')
        savefig(fname2,'pdf','-r100')
        clf reset
    end
%     op_image2mp4(['figures/' fname '/'])
end

end

function doarrows(x,y,u,v)
    %eliminate quiver dots outside of domain
%     skipx = 1; skipy = 1; % these are backwards
%     skipx = 8; skipy = 8; % ai100, horizontal plot
%     skipx = 6; skipy = 6; % ai65
%     skipx = 1; skipy = 36; % Roc jslice
%     skipx = 24; skipy = 14; % Roc zoom in in x zoom1
%     skipx = 24; skipy = 18; % Roc zoom out
%     skipx = 6; skipy = 4; % Roc zoom in 2
    skipx = 2; skipy = 2; % ai65 zoomed usually 2 2 zoom1
%     skipx = 3; skipy = 3; %ai65 less zoomed
    sx = 1:skipx:size(x,1);
    sy = 1:skipy:size(x,2);
    x = x(sx,sy);
    y = y(sx,sy);
    u = u(sx,sy);
    v = v(sx,sy);
    keep = find(isnan(u)~=1 & u~=0 & isnan(v)~=1 & v~=0 ...
    & isnan(x)~=1 & isnan(y)~=1);
    x = x(keep);
    y = y(keep);
    u = u(keep);
    v = v(keep);
    quiver(x,y,u,v,.75,'k')
% % h = quiver2(x,y,1*u,.5*v,.3,'n=',[],... % normal one
% h = quiver2(x,y,1*u,1.5*v,.3,'n=',[],...
% 'c=',[0 0 0],'l@',2,'w@',.25,'b@',45,'t@',60,'t=','quiver','filled');%'a@','fancy','filled'); %black
end