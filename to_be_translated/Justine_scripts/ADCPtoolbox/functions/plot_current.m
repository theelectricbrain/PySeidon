%% plot_current
%
%Description: function to plot currents over measurement period, ensembled
%   around a specified vertical bin.
%
%Inputs:
%   - s: horizontal velocity
%   - t: time (days)
%   - slack_t: times of slack water (days, same reference as t)
%   - w: vertical velocity (optional)

function plot_current(s, t, slack_t, w)

if size(s,2) > 1 || size(w,2) > 1
   disp('plot_current: horizontal and vertical velocity must be for a single bin')
   return
end

plot_init_figure([24 499 1349 437], []);

%plot signed speed as a function of time with slacks marked
hold on
plot(t-t(1), s,'-b')
plot(slack_t-t(1),zeros(size(slack_t)),'og')

%plot vertical velocity over tidal cycle
if ~isempty(w)
    map=colormap;
    miv=min(w);
    mav=max(w);
    clrstep = (mav-miv)/size(map,1);
    
    for nc=1:size(map,1)
        iv = find(w>miv+(nc-1)*clrstep & w<=miv+nc*clrstep) ;
        plot(t(iv)-t(1),s(iv),...
            'o','color',map(nc,:),'markerfacecolor',map(nc,:),'markersize',4);
    end
    
    xlabel('Time (days)','fontweight','b')
    ylabel('Speed (m/s)','fontweight','b')
    
    %add colorbar
    h=colorbar;
    set(h,'ylim',[1 length(map)]);
    
    ticks = 10;
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
    set(get(h,'YLabel'),'String','Vertical velocity (m/s)','fontweight','b')
    
end

end