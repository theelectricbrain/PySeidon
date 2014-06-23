function plot_residual(res, z, bin, n_contour, plot_interval)

dt = (res.t(2)-res.t(1))*24*60; %calculate time step
plot_skip = round(plot_interval/(dt));

plot_init_figure([281 192 1232 589],[])

%plot contours of residual currents
h(1) = subplot('Position',[0.07 0.55 0.9 0.33]);

contourf(res.t(1:plot_skip:end)'-res.t(1),z(1:size(res.s,2)),...
    res.s_filter(1:plot_skip:end,:)',n_contour);

set(gca,'XLim',[0, res.t(end)-res.t(1)])
ylabel('Distance from seabed (m)','fontweight','b')
title('Residual current (m/s)','fontweight','b')

h(2) = colorbar;
v = caxis;
v = max(abs(v));
caxis([-v v])

pcolor_pos = get(gca,'position');

%plot trace of currents at defined bin
%h(2) = subplot('Position',[0.07 0.13 0.815 0.33]);
h(2) = subplot('Position',[0.07 0.13 pcolor_pos(3) pcolor_pos(4)]);
hold on
plot(res.t(1:plot_skip:end)-res.t(1),res.s(1:plot_skip:end,bin),'-b','linewidth',2);
plot([0, res.t(end)-res.t(1)],[0, 0],'-k')
set(gca,'XLim',[0, res.t(end)-res.t(1)])
ylabel('Velocity (m/s)','fontweight','b')
xlabel('Deployment time (days)','fontweight','b')

linkaxes(h,'x')

%alternative plot
plot_init_figure([100   542   768   256],[])

colormap('jet')
%colormap('gray')

contourf(res.t(1:plot_skip:end)',z(1:size(res.s,2)),...
    res.s_filter(1:plot_skip:end,:)',n_contour,'edgecolor','none');

set(gca,'XLim',[0, res.t(end)-res.t(1)])
ylabel('Distance from seabed (m)','fontweight','b')
xlabel('Date','fontweight','b')
datetick('x')

h = colorbar;
v = caxis;
v = max(abs(v));
caxis([-v v])

axes(h)
ylabel('Residual current (m/s)','fontweight','b')

end