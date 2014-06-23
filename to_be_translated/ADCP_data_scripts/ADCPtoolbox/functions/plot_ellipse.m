%Brian Polagye
%6/9/2010

%Description: routine to generate plot of cycle ellipses at a specified
%   depth

%Inputs:
%   - ellipse: structure containing u,v,w for all bins, separated by cycle
%   - plot_w: flag to plot vertical velocity as color on ellipse
%   - bin: bin to plot
%   - offset: number of cycles offset from start of record

%color coding modified from plotclr.m - available for download from
%   Mathworks file exchange

function plot_ellipse(ellipse, plot_w, bin, offset)

num_cycles = 15;

plot_init_figure([52   116   801   832], 10);

%assign color codes for vertical velocities
if plot_w ~= 0
    
    min_w = 0;
    max_w = 0;
    
    for i = offset:offset+num_cycles
        min_w = min(min_w, min(ellipse(i).w(:,bin)));
        max_w = max(max_w, max(ellipse(i).w(:,bin)));
    end
    
    max_w = ceil(max_w*10)/10;
    min_w = floor(min_w*10)/10;
    
    map=colormap;
    clrstep = (max_w-min_w)/size(map,1);
end

%plot num_cycles ellipses
h = zeros(num_cycles,1);
ULim = [0 0];
VLim = [0 0];

for i = 1:num_cycles
    
    h(i)=subplot(5,3,i);    %update if changing number of cycles to plot
    
    if plot_w ~=0
        plot(ellipse(i+offset).u(:,bin),ellipse(i+offset).v(:,bin),'-k')
        hold on
        for nc=1:size(map,1)
            iv = find(ellipse(i+offset).w(:,bin)>min_w+(nc-1)*clrstep & ...
                ellipse(i+offset).w(:,bin)<=min_w+nc*clrstep) ;
            plot(ellipse(i+offset).u(iv,bin),ellipse(i+offset).v(iv,bin),...
                'o','color',map(nc,:),'markerfacecolor',map(nc,:),'markersize',2);
        end
    else
        plot(ellipse(i+offset).u(:,bin),ellipse(i+offset).v(:,bin),'o-b')
    end
    
    ulim = get(gca,'XLim');
    vlim = get(gca,'YLim');
    
    ULim = [min(ulim(1), ULim(1)), max(ulim(2), ULim(2))];
    VLim = [min(vlim(1), VLim(1)), max(vlim(2), VLim(2))];
    
    title(['Cycle ' num2str(i)],'fontweight','b')
    grid on
end

%link axes on all plots
linkaxes(h)
set(h(1),'XLim',ULim)
set(h(1),'YLim',VLim)

%add colorbar
if plot_w ~= 0
    set(gcf,'currentaxes',h(1));
    h=colorbar;
    set(h,'ylim',[1 length(map)]);
    
    ticks = 5;
    pstr = '%-4.2f';
    
    %tick marks
    yal=linspace(1,length(map),ticks);
    set(h,'ytick',yal);
    
    %tick labels
    ytl=linspace(min_w,max_w,ticks);
    s=char(ticks,4);
    for j=1:ticks
        B=sprintf(pstr,ytl(j));
        s(j,1:length(B))=B;
    end
    set(h,'yticklabel',s);
    set(get(h,'YLabel'),'String','Vertical velocity (m/s)','fontweight','b')
    
end

end