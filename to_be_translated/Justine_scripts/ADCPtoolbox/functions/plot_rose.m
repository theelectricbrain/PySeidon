%Brian Polagye
%August 15, 2011

%Description: generate a metric rose at specified depth

function plot_rose(all,bin,z)

plot_init_figure([520   275   643   523],[])

%map ordinates for rose
ord(1) = all.P_mean(bin);   %mean power density - larger is better
ord(2) = all.d_asym(bin);   %directional asymmetry - smaller is better
ord(3) = all.P_asym(bin);   %mean power asymetry - 1 is ideal
ord(4) = all.d_sigma(bin);  %directional variability - smaller is better

ord_title{1} = ['Mean power density' char(10) '(kW/m^2)'];
ord_title{2} = ['Directional' char(10) 'asymmetry (deg)'];
ord_title{3} = 'Power density asymmetry';
ord_title{4} = ['Directional' char(10) 'variability (deg)'];

fill([0 ord(2) 0 -ord(4)],[ord(1) 0 -ord(3) 0],'b')
XLim = get(gca,'XLim');
YLim = get(gca,'YLim');
hold on
plot(XLim, [0 0],'--k')
plot([0 0],YLim,'--k')

xtick = get(gca,'XTick');
set(gca,'XTickLabel',abs(xtick))

ytick = get(gca,'YTick');
set(gca,'YTickLabel',abs(ytick))

%generate axes labels
dx = abs(xtick(2)-xtick(1));
dy = abs(ytick(2)-ytick(1));

y_upper = max(ytick)/2;
y_lower = min(ytick)/2;
x_right = max(xtick)/2;
x_left = min(xtick)/2;

text(min(xtick)-dx/1.5,y_upper,ord_title{1},'rotation',90,'horizontalalignment','center')
text(min(xtick)-dx/1.5,y_lower,ord_title{3},'rotation',90,'horizontalalignment','center')
text(x_left,min(ytick)-dy/1.5,ord_title{2},'horizontalalignment','center')
text(x_right,min(ytick)-dy/1.5,y_lower,ord_title{4},'horizontalalignment','center')

title([num2str(round(z)) 'm Above Seabed'],'fontweight','b')

end