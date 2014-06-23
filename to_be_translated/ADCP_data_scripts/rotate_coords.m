function [xnew,ynew]=rotate_coords(x,y,Theta)
% Similar to "rotate_to_channelcoords.m" code, theta is now the angle
% between the old axis and the new x-axis (CCw is positive)
% 
% Jan 13, 2014

xnew = x.*cos(Theta)+y.*sin(Theta);
ynew = -x.*sin(Theta)+y.*cos(Theta);