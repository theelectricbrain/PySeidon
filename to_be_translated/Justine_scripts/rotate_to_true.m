function [x,y]=rotate_to_true(X,Y,varargin)
% X,Y are the X and Y coordinates (could be speeds) relative to magnetic
% north -- inputs can be vectors
% x,y are the coordinates relative to true north
% This function assumes the measured location is Nova Scotia where the
% declination angle is -19 degrees.
%
% Sept 29, 2012: Changed print statement
%
% Sept 20, 2012: Modified the function to allow for theta to be input.
% Default will remain at -19 degrees, but this may not be accurate for all
% places in Nova Scotia. 
%
%
% May 1, 2012


if nargin == 3
    theta = varargin{1};
elseif nargin == 2
    theta = -19;
else
    error('Wrong number of inputs')
end
s=sprintf('Rotating velocities to be relative to true north (declination = %4.2f)',theta);
disp(s)

Theta = theta*pi/180;

x = X*cos(Theta)+Y*sin(Theta);
y = -X*sin(Theta)+Y*cos(Theta);