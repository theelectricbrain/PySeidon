%Brian L. Polagye
%June 11, 2010

%Description: Generate velocity and power cumulative distribution functions

function plot_operation_CDF(generation, pdf_s, s_bins)

%Inputs:
%   - generation: structure containing probability distribution functions
%       of energy generation as a function of horizontal velocity and speed
%   - pdf_s: probability distribution function of horizontal speed
%   - s_bins: horizontal speed bins for probability distribution function

plot_init_figure([164 328 1180 446], []);

subplot(1,2,1)
hold on
plot(s_bins,cumsum(pdf_s),'-b','linewidth',2)
xlabel('Horizontal Speed (m/s)','fontweight','b')
ylabel('Speed CPDF','fontweight','b')
grid on
set(gca,'YLim',[0 1])

subplot(1,2,2)
plot(s_bins,cumsum(generation.energy_spdf),'-b','linewidth',2)
xlabel('Horizontal Speed (m/s)','fontweight','b')
ylabel('Generation CPDF','fontweight','b')
grid on
set(gca,'YLim',[0 1])

end

