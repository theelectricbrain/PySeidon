%Brian Polagye
%May 2, 2010

%Description: routine to plot depth profiles of velocity characteristics,
%   as measured by ADCP data.

function plot_profiles(all, ebb, fld, z)

plot_init_figure([107 132 1145 766],[])

all_color = 'b';
ebb_color = 'r';
fld_color = 'g';

z_max = size(all.d_mean,1);
z_data = z(1:z_max);

%plots:
%   1. mean speed
%   2. speed asymmetry
%   3. mean power density
%   4. mean direction
%   5. direction standard deviation
%   6. direction asymmetry

subplot(2,3,1)  %mean speed
hold on
plot(all.s_mean, z_data,'-','color',all_color,'linewidth',2)
plot(abs(ebb.s_mean), z_data,'-','color',ebb_color)
plot(abs(fld.s_mean), z_data,'-','color',fld_color)
xlabel('Mean speed (m/s)','fontweight','b')
ylabel('Distance from seabed (m)','fontweight','b')

subplot(2,3,2)  %speed asymmetry
hold on
plot(all.s_asym, z_data,'-','color',all_color,'linewidth',2)
xlabel('Ebb/flood speed asymmetry','fontweight','b')

subplot(2,3,3)  %mean power density
hold on
plot(all.P_mean, z_data,'-','color',all_color,'linewidth',2)
plot(ebb.P_mean, z_data,'-','color',ebb_color)
plot(fld.P_mean, z_data,'-','color',fld_color)
xlabel('Mean power density','fontweight','b')

subplot(2,3,4)  %mean direction
hold on
if mean(all.d_mean)>=0 
    plot(all.d_mean, z_data,'-','color',all_color,'linewidth',2)
else
    plot(all.d_mean+180, z_data,'-','color',all_color,'linewidth',2)
end
if mean(ebb.d_mean)>=0
    plot(ebb.d_mean, z_data,'-','color',ebb_color)
else
   plot(ebb.d_mean+180,z_data,'-','color',ebb_color) 
end
if mean(fld.d_mean)>=0
    plot(fld.d_mean, z_data,'-','color',fld_color)
else
    plot(fld.d_mean+180, z_data,'-','color',fld_color)
end
xlabel('Principal direction','fontweight','b')
ylabel('Distance from seabed (m)','fontweight','b')

subplot(2,3,6)  %direction standard deviation
hold on
plot(all.d_sigma, z_data,'-','color',all_color,'linewidth',2)
plot(ebb.d_sigma, z_data,'-','color',ebb_color)
plot(fld.d_sigma, z_data,'-','color',fld_color)
xlabel('Direction deviation','fontweight','b')
legend('All','Ebb','Flood','location','east')

subplot(2,3,5)  %direction asymmetry
plot(all.d_asym, z_data,'-','color',all_color,'linewidth',2)
xlabel('Ebb/flood direction asymmetry','fontweight','b')

end