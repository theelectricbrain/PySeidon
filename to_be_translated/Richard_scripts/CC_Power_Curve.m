%u=0:0.01:7;
u=speed;
data=find(u>=1.0 & u<=4.5);

nem = 0.0016*u.^4-0.0324*u.^3+ 0.1369*u.^2-0.1534*u + 0.8396;

rho=1025;

A = 27.49;

TSR=4.3;

Cp = -0.0242*TSR^2 + 0.1963*TSR - 0.0049;

P = Cp*nem.*(0.5*rho*A*u.^3);

MeanP=mean(P);
MaxP=max(P);
AEP=MeanP*(365*24);

figure
subplot(1,2,1)
patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',MeanP'/1000)
shading flat;
axis(plot_range)
colorbar
caxis([0 40])
title('Mean Power (kW)')
drawnow

subplot(1,2,2)
patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',MaxP'/1000)
shading flat;
axis(plot_range)
colorbar
caxis([0 200])
title('Max Power (kW)')
drawnow

figure
subplot(1,2,1)
patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',AEP'/1e6)
shading flat;
axis(plot_range)
colorbar
title('AEP (MWh)')
caxis([0 350])

subplot(1,2,2)
patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',h)
shading flat;
axis(plot_range)
colorbar
title('Mean Water Depth')
caxis([0 40])

drawnow
% figure
% plot(u,P/1000,'linewidth',3)
% hold on
% plot(u(data),P(data)/1000,'r','linewidth',3)
% hold off
% xlabel('water speed (m/s)')
% ylabel('Power (kW)')
% axis([0 7 0 400])
% grid
% 
% figure
% plot(u,nem,'linewidth',3)
% hold on
% plot(u(data),nem(data),'r','linewidth',3)
% hold off
% xlabel('water speed (m/s)')
% ylabel('nem')
% axis([0 7 0 1])
% grid
