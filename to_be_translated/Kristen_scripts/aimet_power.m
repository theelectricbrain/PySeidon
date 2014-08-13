function [power,coords] = aimet_power(SD,tind,type,varargin)
% Calculate electrical power after calculating speed in aimet_s.
% Assume constant value for rho in this case. rho=1024 kg/m^3.
% sc, sr, and eta default to .7, 4, and .5, respectively, unless input with
% varargin, in which case they should be labeled then listed, e.g.
%  [power,coords] = aimet_power(SD,tind,'profile',y,x,'sc',1,'eta',.3)

%% Default turbine properties
sc = 0.7; sr = 4; eta = 0.5;

%% Check for input turbine properties
if ~isempty(varargin)
% for i=1:length(varargin)
    ind = strcmp(varargin,{'sc'});
    if sum(ind) sc = varargin{find(ind)+1}; end % cut-in speed
    ind = strcmp(varargin,{'sr'});
    if sum(ind) sr = varargin{find(ind)+1}; end % rated speed
    ind = strcmp(varargin,{'eta'});
    if sum(ind) eta = varargin{find(ind)+1}; end % turbine efficiency
% end
end

%% Constant values
rho = 1024; % kg/m^3

%% Calculate speed
[s,coords] = aimet_s(SD,tind,type,varargin{:});
% load 's.mat'
% s = data; clear data

%% Indices for speeds
inds = (s >= sc).*(s <= sr);
indr = s > sr;

%% Calculate electrical power
power = eta*.5*rho*((s.*inds).^3+(sr.*indr).^3)/1000; % In kW/m^4