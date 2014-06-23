%Brian L. Polagye
%June 11, 2010

%Description: routine to plot statistics for accuracy of power law profile
%   fit to time series profiles

function plot_powerlaw_stats(profile, ds, R2_min, min_profiles)

%Inputs:
%   - profile: structure containing depth average speed (s), power law
%       exponent (exp), and R2 statistic (R2) for each profile
%   - ds: horizontal velocity increment for binning profiles (0.2 m/s default)
%   - R2_min: minimum R2 value for an acceptable fit (0.5 default)
%   - min_profiles: minimum number of profiles per bin (30 default)

%% Calculate statistics for each stage of the tide, as defined by depth-average velocity

if isempty(R2_min), R2_min = 0.5; end
if isempty(ds), ds = 0.2; end
if isempty(min_profiles), min_profiles = 30; end

profile.s_bins = floor(nanmin(profile.s)):ds:ceil(nanmax(profile.s));

profile.s_mean = zeros(length(profile.s_bins)-1,1);
profile.exp_mean = zeros(size(profile.s_mean));
profile.exp_std = zeros(size(profile.s_mean));
profile.exp_fit = zeros(size(profile.s_mean));

for i = 1:length(profile.s_bins)-1
    ind_bin = find(profile.s>=profile.s_bins(i) & profile.s<profile.s_bins(i+1));
    ind_count = find(profile.s>=profile.s_bins(i) & profile.s<profile.s_bins(i+1) ...
        & profile.R2>R2_min & ~isnan(profile.exp));
    if length(ind_bin) > min_profiles
        profile.s_mean(i) = (profile.s_bins(i)+profile.s_bins(i+1))/2;
        profile.exp_mean(i) = nanmean(profile.exp(ind_count));
        profile.exp_std(i) = nanstd(profile.exp(ind_count));
        profile.exp_fit(i) = length(ind_count)/length(ind_bin);
    end
end

%% Plot profile statistics

plot_init_figure([195 268 1101 415], []);

errorbar(profile.s_mean,(profile.exp_mean),(profile.exp_std),'k','linestyle','none');
hold on

%add colors denoting % fit
map=colormap('Copper');
map = map(end:-1:1,:);
colormap(map);
miv=0;
mav=1;
clrstep = (mav-miv)/size(map,1);

for nc=1:size(map,1)
    iv = find(profile.exp_fit>miv+(nc-1)*clrstep & profile.exp_fit<=miv+nc*clrstep) ;
    plot(profile.s_mean(iv),profile.exp_mean(iv),...
        'o','color','k','markerfacecolor',map(nc,:),'markersize',8);
end

ylabel('1/Power Law Exponent','fontweight','b')
xlabel('Depth-averaged Velocity (m/s)','fontweight','b')
set(gca,'XLim',[min(profile.s_bins) max(profile.s_bins)])
YLim = get(gca,'YLim');
set(gca,'YLim',[0 YLim(2)])

%add colorbar
h=colorbar;

ticks = 5;
pstr = '%-4.2f';

%tick marks
yal=linspace(1,length(map),ticks);
set(h,'ytick',yal);

%tick labels
ytl=linspace(miv,mav,ticks);
s=char(ticks,4);
for j=1:ticks
    B=sprintf(pstr,ytl(j));
    s(j,1:length(B))=B;
end
set(h,'yticklabel',s);

%add colorbar label
axes(h)
ylabel('Fraction Good Fits','fontweight','b','fontname','Times','fontsize',8)

end