%Brian L. Polagye
%June 11, 2010

%Description: routine to plot operational times and statistic histograms

function plot_operation_period(operation,cut_in, s, t)

%Inputs:
%   - operation: structure containing operating windows and fraction of
%       time idle/operating for a given duration
%   - cut_in: turbine cut-in speed
%   - s: horizontal current velocity time series
%   - t: accompanying time stamps

plot_init_figure([35 283 1330 422], []);

hold on
plot(t-t(1),s,'-b','linewidth',2)                 %plot current trace
plot([0,t(end)-t(1)],[cut_in, cut_in],'--k')      %superimpose cut-in speed (flood)
plot([0,t(end)-t(1)],[-cut_in, -cut_in],'--k')    %superimpose cut-in speed (ebb)
xlabel('Time (days)','fontweight','b')
ylabel('Velocity (m/s)','fontweight','b')

YLim = get(gca,'YLim');    %store y-axis limits

T_cum = cumsum(operation.windows(:,1))/24; %cumulative time in days
T_cum = vertcat(0, T_cum);
for i = 2:size(operation.windows,1)
    if operation.windows(i-1,2) == 0
        fill([T_cum(i-1) T_cum(i-1) T_cum(i) T_cum(i)],...
            [YLim,fliplr(YLim)],[205/255,0,0],'edgecolor','none','facealpha',0.25)
    else
        fill([T_cum(i-1) T_cum(i-1) T_cum(i) T_cum(i)],...
            [YLim,fliplr(YLim)],[0,205/255,0],'edgecolor','none','facealpha',0.25)       
    end
   
end

plot(t-t(1),s,'-b','linewidth',2)    %plot current trace over operational periods

%plot histograms of operational and idle durations
plot_init_figure([220 345 1060 377], []);
    
%XLimit definition
XLim = min(24, max(operation.T_bins));

subplot(1,2,1)
bar(operation.T_bins,operation.f_off)
title(['Turbine Idle (cut-in = ' num2str(cut_in) 'm/s)'],'fontweight','b')
xlabel('Idle Duration (h)','fontweight','b')
ylabel('Frequency','fontweight','b')
set(gca,'XLim',[-1 XLim])

subplot(1,2,2)
bar(operation.T_bins,operation.f_on)
title(['Turbine Operating (cut-in = ' num2str(cut_in) 'm/s)'],'fontweight','b')
xlabel('Operating Duration (h)','fontweight','b')
ylabel('Frequency','fontweight','b')
set(gca,'XLim',[-1 XLim])

end