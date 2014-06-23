function theta=get_DirFromN(u,v)
%% function theta=get_DirFromN(u,v)
%This function computes the direction from North with the output in degrees
%and measured clockwise from north. 
% 
% Inputs:
%   u: eastward component
%   v: northward component
%
% Justine McMillan
% April 2012

theta=atan2(u,v)*180/pi;

ind = find(theta<0);
theta(ind) = theta(ind)+360;