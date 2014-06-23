%Brian Polagye
%April 23, 2010

%Description: initialize default properties for new figures

function plot_init_figure(figure_size, font_size)

if isempty(font_size);
    font_size = 12;
end

figure
clf
set(gcf,'color','w')
set(gcf,'defaulttextfontsize',font_size,'defaulttextfontname','Arial')
set(gcf,'defaultaxesfontsize',font_size,'defaultaxesfontname','Arial')

if isempty(figure_size)
    set(gcf,'position',[3 108 1060 840])
else
    set(gcf,'position',figure_size)
end