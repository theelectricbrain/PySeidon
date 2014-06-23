%% plot_pcolor

%Description: function to plot current magnitude and direction over measurement period

%Inputs:
%   - s: horizontal velocity
%   - t: time (days)
%   - h: free surface elevation
%   - d: current direction
%   - duration: length of record to plot (days)

function plot_pcolor(s, t, h, d, z, duration)

plot_init_figure([41 231 1349 533], []);

%colormap('gray')

if duration > t(end)-t(1)
    disp('Warning: specified duration for plot_pcolor exceeds record length')
    return
else
    ind_end = find(t-t(1)<=duration,1,'last');
    t = t(1:ind_end);
    s = s(1:ind_end,:);
    h = h(1:ind_end);
    d = d(1:ind_end,:);
end

%magnitude
subplot(2,1,1)

plot(t-t(1),h,'-k')
hold on
pcolor(t-t(1), z(1:size(s,2)), abs(s)');
shading flat
plot(t-t(1),h,'-k')

g = colorbar;
ylabel('Distance from seabed (m)','fontweight','b')
box off

axes(g)
ylabel('Velocity magnitude (m/s)','fontweight','b','fontname','times','fontsize',10)

%direction
subplot(2,1,2)
hold on

pcolor(t-t(1), z(1:size(s,2)), d');
shading flat
plot(t-t(1),h,'-k')

g = colorbar;
ylabel('Distance from seabed (m)','fontweight','b')
xlabel('Observation time (days)','fontweight','b')

axes(g)
ylabel('Velocity direction (^o)','fontweight','b','fontname','times','fontsize',10)

end