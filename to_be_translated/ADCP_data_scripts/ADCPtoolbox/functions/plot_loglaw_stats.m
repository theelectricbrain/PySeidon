%Brian Polagye
%November 7, 2011

%Description: plot logarithmic profile statistics

function plot_loglaw_stats(profile, t, ds, R2_min, min_profiles, min_percent_fit)

%Inputs:
%   - profile: structure containing depth average speed (s), bottom roughness (z0), 
%       shear velocity (ushear), and R2 statistic (R2) for each profile
%   - ds: horizontal velocity increment for binning profiles (0.2 m/s default)
%   - R2_min: minimum R2 value for an acceptable fit (0.9 default)
%   - min_profiles: minimum number of profiles per bin (30 default)
%   - min_percent_fit: minimum percentage of profiles fit by a log law in a
%   particular bin to consider the log law intepretation to have some
%   validity to describe condition (0.3 default)

%% Calculate statistics for each stage of the tide, as defined by depth-average velocity

%set default values
if isempty(R2_min), R2_min = 0.95; end
if isempty(ds), ds = 0.2; end
if isempty(min_profiles), min_profiles = 30; end
if isempty(min_percent_fit), min_percent_fit = 0.5; end

profile.s_bins = floor(nanmin(profile.s)):ds:ceil(nanmax(profile.s));

profile.s_mean = NaN(length(profile.s_bins)-1,1);
profile.z0_mean = NaN(size(profile.s_mean));
profile.z0_std = NaN(size(profile.s_mean));
profile.ushear_mean = NaN(size(profile.s_mean));
profile.ushear_std = NaN(size(profile.s_mean));
profile.tau_mean = NaN(size(profile.s_mean));
profile.tau_std = NaN(size(profile.s_mean));
profile.tau_min = NaN(size(profile.s_mean));
profile.tau_max = NaN(size(profile.s_mean));
profile.z_log_mean = NaN(size(profile.s_mean));
profile.z_log_std = NaN(size(profile.s_mean));
profile.fit = zeros(size(profile.s_mean));

for i = 1:length(profile.s_bins)-1
    ind_bin = find(profile.s>=profile.s_bins(i) & profile.s<profile.s_bins(i+1));
    ind_count = find(profile.s>=profile.s_bins(i) & profile.s<profile.s_bins(i+1) ...
        & profile.R2(:,1)>R2_min & ~isnan(profile.z0(:,1)));
    if length(ind_bin) > min_profiles && ~isempty(ind_count)
        profile.s_mean(i) = (profile.s_bins(i)+profile.s_bins(i+1))/2;
        
        profile.z0_mean(i) = nanmean(profile.z0(ind_count,1));
        profile.z0_std(i) = nanstd(profile.z0(ind_count,1));
        
        profile.ushear_mean(i) = nanmean(profile.ushear(ind_count,1));
        profile.ushear_std(i) = nanstd(profile.ushear(ind_count,1));
        
        profile.tau_mean(i) = nanmean(profile.tau(ind_count,1));
        profile.tau_std(i) = nanstd(profile.tau(ind_count,1));
        profile.tau_min(i) = nanmin(profile.tau(ind_count,1));
        profile.tau_max(i) = nanmax(profile.tau(ind_count,1));
        
        profile.z_log_mean(i) = nanmean(profile.z_log(ind_count,1));
        profile.z_log_std(i) = nanstd(profile.z_log(ind_count,1));
        
        profile.fit(i) = length(ind_count)/length(ind_bin);
    end
end

%% Plot time series

plot_init_figure([195 91 1101 686], []);

%R2 statistic
h(1) = subplot(5,1,1);
plot(t(profile.R2(:,1)>R2_min)-t(1),profile.R2(profile.R2(:,1)>R2_min,1),'ob')
hold on
plot(t(profile.R2(:,1)<=R2_min)-t(1),profile.R2(profile.R2(:,1)<=R2_min,1),'o','markeredgecolor',[192/255,192/255,192/255])
set(gca,'YLim',[0 1])
ylabel('R^2','fontweight','b')

%shear velocity
h(2) = subplot(5,1,2);
plot(t(profile.R2(:,1)>R2_min)-t(1),abs(profile.ushear(profile.R2(:,1)>R2_min,1)),'ob')
ylabel('u* (m s^{-1})','fontweight','b')

%bottom roughness
h(3) = subplot(5,1,3);
semilogy(t(profile.R2(:,1)>R2_min & profile.z0(:,1)>=0)-t(1),profile.z0(profile.z0(:,1)>=0 & profile.R2(:,1)>R2_min,1),'ob')
xlabel('Deployment time (days)','fontweight','b')
ylabel('z_0 (m)','fontweight','b')

%shear stress
h(4) = subplot(5,1,4);
plot(t(profile.R2(:,1)>R2_min)-t(1),abs(profile.tau(profile.R2(:,1)>R2_min,1)),'ob')
ylabel('\tau_b (N m^{-2})','fontweight','b')
linkaxes(h,'x')

%shear stress
h(5) = subplot(5,1,5);
plot(t(profile.R2(:,1)>R2_min)-t(1),abs(profile.z_log(profile.R2(:,1)>R2_min,1)),'ob')
ylabel('z_{log} (m)','fontweight','b')
linkaxes(h,'x')

%% Plot shear velocity, roughness, and shear stress as function of tidal current velocity

plot_init_figure([ 270 71 1101 735], []);
h = [];

%percentage fit
h(1) = subplot(5,1,1);
plot(profile.s_mean, profile.fit,'ok','linestyle','none')
hold on
plot(profile.s_mean,repmat(min_percent_fit,length(profile.s_mean),1),'--k')
plot(profile.s_mean(profile.fit<min_percent_fit),profile.fit(profile.fit<min_percent_fit),'xr','linestyle','none')
ylabel(['R^2 > ' num2str(R2_min)],'fontweight','b')

%shear velocity
h(2) = subplot(5,1,2);
errorbar(profile.s_mean,abs(profile.ushear_mean),(profile.ushear_std),'ok','linestyle','none');
hold on
plot(profile.s_mean(profile.fit<min_percent_fit),abs(profile.ushear_mean(profile.fit<min_percent_fit))...
    ,'xr','linestyle','none')
ylabel('u* (m s^{-1})','fontweight','b')
set(gca,'YLim',[0 max(get(gca,'YLim'))])

%bottom roughness
h(3) = subplot(5,1,3);
errorbar(profile.s_mean,(profile.z0_mean),(profile.z0_std),'ok','linestyle','none');
hold on
plot(profile.s_mean(profile.fit<min_percent_fit),(profile.z0_mean(profile.fit<min_percent_fit)),...
    'xr','linestyle','none')
ylabel('z_0 (m)','fontweight','b')
set(gca,'YLim',[0 max(get(gca,'YLim'))])

%shear stress
h(4) = subplot(5,1,4);
%errorbar(profile.s_mean,(profile.tau_mean),(profile.tau_std),'ok','linestyle','none');
errorbar(profile.s_mean,profile.tau_mean,profile.tau_std,'ok','linestyle','none')
hold on
errorbar(profile.s_mean,profile.tau_mean,profile.tau_mean-profile.tau_min,profile.tau_max-profile.tau_mean,'ok','linestyle','none')
plot(profile.s_mean(profile.fit<min_percent_fit),profile.tau_mean(profile.fit<min_percent_fit),'xr','linestyle','none')
ylabel('\tau_b (N m^{-2})','fontweight','b')
%xlabel('Current velocity (m s^{-1})','fontweight','b')
set(gca,'YLim',[0 max(get(gca,'YLim'))])

%log layer height
h(5) = subplot(5,1,5);
errorbar(profile.s_mean,(profile.z_log_mean),(profile.z_log_std),'ok','linestyle','none');
hold on
plot(profile.s_mean(profile.fit<min_percent_fit),profile.z_log_mean(profile.fit<min_percent_fit),'xr','linestyle','none')
ylabel('z_{log} (m)','fontweight','b')
xlabel('Current velocity (m s^{-1})','fontweight','b')
set(gca,'YLim',[0 max(get(gca,'YLim'))])

linkaxes(h,'x')

%% Debugging - test for normal distribution
% 
% plot_bin = 10;
% 
% ind_count = find(profile.s>=profile.s_bins(plot_bin) & profile.s<profile.s_bins(plot_bin+1) ...
%         & profile.R2>R2_min & ~isnan(profile.z0));
% 
% figure(3)
% clf
% normplot(profile.ushear(ind_count))
% title(['Velocity = ' num2str(profile.s_bins(plot_bin)) ' m/s (N = ' num2str(length(ind_count)) ')'])

%% Shear stress versus depth-averaged velocity

figure(3)
clf
ind_ebb = find(profile.s<=0 & profile.R2(:,1) >=R2_min);
ind_fld = find(profile.s>=0 & profile.R2(:,1) >=R2_min);

semilogy(abs(profile.s(ind_ebb)),profile.tau(ind_ebb,1),'ok','markerface','k','markersize',4)
hold on
semilogy(profile.s(ind_fld),profile.tau(ind_fld,1),'or','markerface','r','markersize',4)
legend('ebb','flood','location','southeast')
xlabel('Depth-average Velocity (m s^{-1})','fontweight','b')
ylabel('Shear stress (N m^{-2})','fontweight','b')

end

function add_colors(val_fit, s, val,add_colorbar)

%add colors denoting % fit
map=colormap('Copper');
map = map(end:-1:1,:);
colormap(map);
miv=0;
mav=1;
clrstep = (mav-miv)/size(map,1);

for nc=1:size(map,1)
    iv = find(val_fit>miv+(nc-1)*clrstep & val_fit<=miv+nc*clrstep) ;
    plot(s(iv),val(iv),...
        'o','color','k','markerfacecolor',map(nc,:),'markersize',8);
end

if add_colorbar
    %add colorbar
    h=colorbar('location','SouthOutside');
    
    ticks = 5;
    pstr = '%-4.2f';
    
    %tick marks
    yal=linspace(1,length(map),ticks);
    set(h,'xtick',yal);
    
    %tick labels
    ytl=linspace(miv,mav,ticks);
    s=char(ticks,4);
    for j=1:ticks
        B=sprintf(pstr,ytl(j));
        s(j,1:length(B))=B;
    end
    set(h,'xticklabel',s);
    
    %add colorbar label
    axes(h)
    xlabel('Fraction Good Fits','fontweight','b','fontname','Times','fontsize',8)
    
end

end