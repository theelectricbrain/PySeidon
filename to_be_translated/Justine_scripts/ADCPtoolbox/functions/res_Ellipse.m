%Brian Polagye
%November 11, 2011

%Justine McMillan
%May 31,2012 -- modified i index -- was getting errors about sizes (need tp
%confirm that this does not introduce errors)

%Description: routine to generate data for progressive vector diagrams over multiple tidal cycles

function [ellipse] = res_Ellipse(u, v, w, slack_ind)

%Inputs:
%   - u: east velocity
%   - v: north velocity
%   - w: vertical velocity
%   - slack_ind: indices denoting slack water

%Outputs: ellipse structure containing u,v,w components over pairs of ebb/flood cycles

for i = 1:2:length(slack_ind)-2
    
 
    
    for j = 1:size(u,2)
        ellipse((i+1)/2).u(:,j) = u(slack_ind(i):slack_ind(i+2),j);
        ellipse((i+1)/2).v(:,j) = v(slack_ind(i):slack_ind(i+2),j);
        ellipse((i+1)/2).w(:,j) = w(slack_ind(i):slack_ind(i+2),j);
    end
end

if isempty(ellipse(end).u)
    ellipse = ellipse(1:end-1);
end

end