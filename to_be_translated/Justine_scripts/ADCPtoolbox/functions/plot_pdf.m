%Brian Polagye
%August 15, 2011

%Description: plots of probability distribution functions

function plot_pdf(pdf_s, jpdf_sd, s_bins, d_bins, bin, z)

%plot PDF and CPDF of horizontal speed
plot_init_figure([211   378   869   420],[])

subplot(1,2,1)
bar(s_bins,pdf_s(:,bin),'edgecolor','none','facecolor','b','barwidth',0.5)
set(gca,'XLim',[s_bins(1) s_bins(end)])
xlabel('Horizontal Velocity (m/s)','fontweight','b')
ylabel('Fraction Occurrence','fontweight','b')

subplot(1,2,2)
plot(s_bins,cumsum(pdf_s(:,bin)),'-b','linewidth',2)
set(gca,'XLim',[s_bins(1) s_bins(end)])
xlabel('Horizontal Velocity (m/s)','fontweight','b')
ylabel('Cumulative Fraction Occurrence','fontweight','b')

%plot PDF and CPDF of horizontal velocity
s_abs_bins = unique(abs(s_bins));
for i = 1:length(s_abs_bins)
   ind = find(abs(s_bins)==s_abs_bins(i));
   pdf_s_abs(i) =  sum(pdf_s(ind,bin));
end

plot_init_figure([211   378   869   420],[])

subplot(1,2,1)
bar(s_abs_bins,pdf_s_abs)
set(gca,'XLim',[s_abs_bins(1) s_abs_bins(end)])
xlabel('Horizontal Speed (m/s)','fontweight','b')
ylabel('Fraction Occurrence','fontweight','b')

subplot(1,2,2)
plot(s_abs_bins,cumsum(pdf_s_abs),'-b','linewidth',2)
set(gca,'XLim',[s_abs_bins(1) s_abs_bins(end)])
xlabel('Horizontal Speed (m/s)','fontweight','b')
ylabel('Cumulative Fraction Occurrence','fontweight','b')

%plot JPDF of direction with horizontal velocity

plot_init_figure([520   291   663   507],[])

x = jpdf_sd(:,:,bin);
x(x==0)=NaN;

pcolor(d_bins,s_bins',x)
shading flat
xlabel('Direction (degrees true)','fontweight','b')
ylabel('Horizontal Velocity (m/s)','fontweight','b')
title([num2str(round(z)) 'm Above Seabed'],'fontweight','b')

hold on
plot([0 360],[0 0],'--k')

h = colorbar;
axes(h)
ylabel('Fraction Occurrence','fontweight','b')

%plot cumulative distributions of velocity and power density

p_bins = (s_abs_bins.^3);
cpdf_p = cumsum(pdf_s_abs.*p_bins/sum(pdf_s_abs.*p_bins));
cpdf_s = cumsum(pdf_s_abs);

plot_init_figure([211   378   869   420],[])

subplot(1,2,1);
plot(s_abs_bins,cpdf_s,'-k','linewidth',2)
xlabel('Horizontal velocity (m/s)','fontweight','b')
ylabel('Velocity CPDF','fontweight','b')
set(gca,'Ylim',[0 1])
grid on

subplot(1,2,2)
plot(s_abs_bins,cpdf_p,'-k','linewidth',2)
xlabel('Horizontal velocity (m/s)','fontweight','b')
ylabel('Power density CPDF','fontweight','b')
set(gca,'Ylim',[0 1])
grid on

end