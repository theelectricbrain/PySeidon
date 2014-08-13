function [preal,coords] = aimet_preal(SD,tind,type,varargin)
% Calculate actual electrical power based on turbine efficiency and cut-in 
% and rated speeds.
% Assume constant value for rho in this case. rho=1024 kg/m^3.

%% Power based on input turbine properties
[num,~] = aimet_power(SD,tind,type,varargin{:});
% load 'power.mat'
% num = data; clear data coords

% %% Power based on sc=0, sr=infinity
% % Check for input turbine properties
% ind = strcmp(varargin,{'sc'});
% if sum(ind)  % cut-in speed
%     varargin{find(ind)+1}=0; 
% else
%     varargin = cat(2,varargin,{'sc' 0});
% end
% ind = strcmp(varargin,{'sr'});
% if sum(ind)  % rated speed
%     varargin{find(ind)+1}=100; 
% else
%     varargin = cat(2,varargin,{'sr' 100});
% end
% ind = strcmp(varargin,{'eta'});
% if sum(ind) % turbine efficiency
%     varargin{find(ind)+1}=1; 
% else
%     varargin = cat(2,varargin,{'eta' .5});
% end
eta = .5;
% 
[denom,coords] = aimet_kpd(SD,tind,type,varargin{:});
% load 'kpd.mat' %kpd is 'data'
% coords = op_elimtdim(coords); % Fix coords to match
% denom = data/1000; clear data

%% preal
num = op_taverage(num,SD,tind); % Calc average
denom = op_taverage(denom,SD,tind); % Calc average
preal = num./(eta*denom);