% Brian Polagye
% November 11, 2011

function [std_fld, std_ebb, std_all] = dir_std(s, min_speed, d_PA, u, v, d_fld, d_ebb, show)

%Description: compute the standard deviation of the direction of the ebb
%   and flood tides relative to the principal axes
%   The composite tide uses a weighted average of ebb and flood values.
%   Calculation is only for speeds in excess of cut-in floor to avoid
%   weighting statistics towards meaningless variations around slack water.

%Input
% - z: vertical bins to process
% - d_PA: current direction aligned to principal axes (all bins)
% - s: signed current speed (all bins)
% - min_speed: minimum speed to include for analysis
% - u: x-component of horizontal velocity
% - v: y-component of horizontal velocity
% - d_fld, d_ebb: principal axis direction for ebb and flood

std_all = zeros(size(s,2),1);
std_fld = zeros(size(s,2),1);
std_ebb = zeros(size(s,2),1);

%plot scatter plots to verify axes alignment visually
if show, plot_init_figure([30 127 950 821], []); end

for i = 1:size(s,2) %loop through all depth bins

    %determine standard deviation for ebb and flood
    ind_fld = find(s(:,i)>=min_speed);
    std_fld(i) = std(d_PA(ind_fld,i));

    ind_ebb = find(s(:,i)<=-min_speed);
    std_ebb(i) = std(d_PA(ind_ebb,i));

    %compute weighted average
    std_all(i) = (std_fld(i)*length(ind_fld) + std_ebb(i)*length(ind_ebb))...
        /(length(ind_ebb)+length(ind_fld));

    %generate plot to visually confirm axis assignment
    if show
        clf
        hold on
        plot(u(ind_ebb,i),v(ind_ebb,i),'ob','markersize',2)
        plot(u(ind_fld,i),v(ind_fld,i),'or','markersize',2)
        r_ebb = -min(s(:,i));
        r_fld = max(s(:,i));
        %principal axis
        plot([0 r_ebb*cos(-d_ebb(i)*pi/180+pi/2)],[0 r_ebb*sin(-d_ebb(i)*pi/180+pi/2)],'-k','linewidth',2)
        plot([0 r_fld*cos(-d_fld(i)*pi/180+pi/2)],[0 r_fld*sin(-d_fld(i)*pi/180+pi/2)],'-k','linewidth',2)
        %ebb standard deviation
        plot([0 r_ebb*cos(-d_ebb(i)*pi/180+pi/2+std_ebb(i)*pi/180)],...
            [0 r_ebb*sin(-d_ebb(i)*pi/180+pi/2+std_ebb(i)*pi/180)],'--g','linewidth',2)
        plot([0 r_ebb*cos(-d_ebb(i)*pi/180+pi/2-std_ebb(i)*pi/180)],...
            [0 r_ebb*sin(-d_ebb(i)*pi/180+pi/2-std_ebb(i)*pi/180)],'--g','linewidth',2)
        %flood standard deviation
        plot([0 r_fld*cos(-d_fld(i)*pi/180+pi/2+std_fld(i)*pi/180)],...
            [0 r_fld*sin(-d_fld(i)*pi/180+pi/2+std_fld(i)*pi/180)],'--g','linewidth',2)
        plot([0 r_fld*cos(-d_fld(i)*pi/180+pi/2-std_fld(i)*pi/180)],...
            [0 r_fld*sin(-d_fld(i)*pi/180+pi/2-std_fld(i)*pi/180)],'--g','linewidth',2)
        
        xlabel('u-velocity (m/s)','fontweight','b')
        ylabel('v-velocity (m/s)','fontweight','b')
        title(['bin number ', num2str(i)])
        axis equal
        
        pause(0.25)     %brief pause to view results 
       %pause
       
    end

end

end