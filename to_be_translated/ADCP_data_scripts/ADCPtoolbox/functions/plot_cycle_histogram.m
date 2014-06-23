%Brian L. Polagye
%June 11, 2010

%Description: plot histograms of cycle duration, ebb amplitude, flood
%   amplitude, and 

function plot_cycle_histogram(cycle, bin, z, ds, dT, dd)

%Inputs:
%   cycle: statistics about each tidal cycle
%   bin: vertical bin for plotting
%   z: coordinates relative to seabed of vertical bins
%   ds: bin increment for horizontal velocity (m/s)
%   dT: bin increment for cycle duration (h)
%   dd: bin increment for angle (degrees)

plot_init_figure([3 108 1060 840],[]);

%% Generate histogram data for duration and intensity

%cycle duration
T_bin_min = floor(nanmin(reshape(cycle.T_cycle(:,bin),numel(cycle.T_cycle(:,bin)),1)));
T_bin_max = ceil(nanmax(reshape(cycle.T_cycle(:,bin),numel(cycle.T_cycle(:,bin)),1)));
T_bins = T_bin_min:dT:T_bin_max;

%velocity
s_bin_min = floor(nanmin(reshape(cycle.u_cycle(:,bin),numel(cycle.u_cycle(:,bin)),1)));
s_bin_max = ceil(nanmax(reshape(cycle.u_cycle(:,bin),numel(cycle.u_cycle(:,bin)),1)));
ebb_bins = s_bin_min:ds:0;
flood_bins = 0:ds:s_bin_max;
 
%create histogram counts - probability distribution
T_count = hist(cycle.T_cycle(:,bin),T_bins);                                %cycle duration
ebb_count = hist(cycle.u_cycle(cycle.u_cycle(:,bin)<0,bin),ebb_bins);       %ebb amplitude
flood_count = hist(cycle.u_cycle(cycle.u_cycle(:,bin)>=0,bin),flood_bins);  %flood amplitude

%% Plot histograms - duration and intensity

subplot(3,1,1)
bar(T_bins,T_count/sum(T_count))
xlabel('Cycle Period (h)','fontweight','b')
ylabel('Probability','fontweight','b')
title([num2str(round(z(bin))) 'm Above Seabed'],'fontweight','b')
XLim = get(gca,'XLim');
set(gca,'XLim',[0 min(24,XLim(2))])

subplot(3,1,2)
bar(ebb_bins,ebb_count/sum(ebb_count))
xlim([min(s_bin_min), 0])
xlabel('Peak Ebb Velocity (m/s)','fontweight','b')
ylabel('Probability','fontweight','b')

subplot(3,1,3)
bar(flood_bins,flood_count/sum(flood_count))
xlim([0, max(s_bin_max)])
xlabel('Peak Flood Velocity (m/s)','fontweight','b')
ylabel('Probability','fontweight','b')

%% Joint probability distribution for cycle direction and speed

%determine limits of distribution for speed
smax = ceil(max(abs(cycle.u_cycle(:,bin))));

s_bins_fld = [0:ds:smax]';
s_bins_ebb = [-smax:ds:0]';

ind_fld = find(cycle.u_cycle(:,bin)>0);
ind_ebb = find(cycle.u_cycle(:,bin)<0);

%determine limits of distribution for direction
DLim_fld = [round(min(cycle.davg_cycle(ind_fld,bin))/10)*10, ...
    round(max(cycle.davg_cycle(ind_fld,bin))/10)*10];
DLim_ebb = [round(min(cycle.davg_cycle(ind_ebb,bin))/10)*10, ...
    round(max(cycle.davg_cycle(ind_ebb,bin))/10)*10];

d_bins_fld = [DLim_fld(1):dd:DLim_fld(2)]';
d_bins_ebb = [DLim_ebb(1):dd:DLim_ebb(2)]';

N_fld = hist3([cycle.u_cycle(ind_fld,bin), cycle.davg_cycle(ind_fld,bin)],{s_bins_fld,d_bins_fld});
N_ebb = hist3([cycle.u_cycle(ind_ebb,bin), cycle.davg_cycle(ind_ebb,bin)],{s_bins_ebb,d_bins_ebb});

%calculate % of time particular direction occurs for each cycle amplitude
%c_fld = N_fld./repmat(sum(N_fld,2),1,size(N_fld,2));
%c_ebb = N_ebb./repmat(sum(N_ebb,2),1,size(N_ebb,2));

%calculate % of time particular direction deviation occurs
c_fld = N_fld/sum(sum(N_fld));
c_ebb = N_ebb/sum(sum(N_ebb));

%blank out cells with zero count
c_fld(c_fld==0)=NaN;
c_ebb(c_ebb==0)=NaN;

%plot distribution
plot_init_figure([180          49        1176         768],[])

subplot(2,2,1)
pcolor(s_bins_fld,d_bins_fld',c_fld')
shading flat
xlabel('Cycle amplitude (m/s)','fontweight','b')
ylabel(['Average direction in cycle' char(10) '(deg rel. to principal axis)'],'fontweight','b')
title('Flood','fontweight','b')
h = colorbar;
caxis([0 ceil(max(max((c_fld)))*10)/10])
axes(h)
ylabel('Fraction occurrence','fontweight','b')

subplot(2,2,2)
pcolor(abs(s_bins_ebb),d_bins_ebb',c_ebb')
shading flat
h = colorbar;
XTickLabel = get(gca,'XTick')*-1;
set(gca,'XTickLabel',XTickLabel)
xlabel('Cycle amplitude (m/s)','fontweight','b')
ylabel(['Average direction in cycle' char(10) '(deg rel. to principal axis)'],'fontweight','b')
title('Ebb','fontweight','b')
caxis([0 ceil(max(max((c_fld)))*10)/10])
axes(h)
ylabel('Fraction occurrence','fontweight','b')

%% Joint probability distribution for cycle direction deviation and speed

%determine limits of distribution
DLim = round(max(cycle.dstd_cycle(:,bin))/10)*10;
d_bins_fld = [0:dd:DLim]';
d_bins_ebb = [0:dd:DLim]';

N_fld = hist3([cycle.u_cycle(ind_fld,bin), cycle.dstd_cycle(ind_fld,bin)],{s_bins_fld,d_bins_fld});
N_ebb = hist3([cycle.u_cycle(ind_ebb,bin), cycle.dstd_cycle(ind_ebb,bin)],{s_bins_ebb,d_bins_ebb});

%calculate % of time particular direction deviation occurs for each cycle amplitude
%c_fld = N_fld./repmat(sum(N_fld,2),1,size(N_fld,2));
%c_ebb = N_ebb./repmat(sum(N_ebb,2),1,size(N_ebb,2));

%calculate % of time particular direction deviation occurs
c_fld = N_fld/sum(sum(N_fld));
c_ebb = N_ebb/sum(sum(N_ebb));

c_fld(c_fld==0)=NaN;
c_ebb(c_ebb==0)=NaN;

%add to existing plot

subplot(2,2,3)
pcolor(s_bins_fld,d_bins_fld',c_fld')
shading flat
xlabel('Cycle amplitude (m/s)','fontweight','b')
ylabel('Standard deviation of direction in cycle (deg)','fontweight','b')
title('Flood','fontweight','b')
h = colorbar;
caxis([0 ceil(max(max((c_fld)))*10)/10])
axes(h)
ylabel('Fraction occurrence','fontweight','b')

subplot(2,2,4)
pcolor(abs(s_bins_ebb),d_bins_ebb',c_ebb')
shading flat
h = colorbar;
XTickLabel = get(gca,'XTick')*-1;
set(gca,'XTickLabel',XTickLabel)
xlabel('Cycle amplitude (m/s)','fontweight','b')
ylabel('Standard deviation of direction in cycle (deg)','fontweight','b')
title('Ebb','fontweight','b')
caxis([0 ceil(max(max((c_fld)))*10)/10])
axes(h)
ylabel('Fraction occurrence','fontweight','b')

end