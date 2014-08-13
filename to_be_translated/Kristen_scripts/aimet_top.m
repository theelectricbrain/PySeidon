function [top,coords] = aimet_top(SD,tind,type,varargin)
% Calculate turbine operation time.
% sc, sr, and eta default to .7, 4, and .5, respectively, unless input with
% varargin, in which case they should be labeled then listed, e.g.
%  [top,coords] = aimet_top(SD,tind,'profile',y,x,'sc',1,'eta',.3)

%% Default turbine properties
sc = 0.7;

%% Check for input turbine properties
% for i=1:length(varargin)
    ind = strcmp(varargin,{'sc'});
    if sum(ind)  % cut-in speed
        sc = varargin{find(ind)+1}; 
    end
% end

%% Calculate speed
[s,coords] = aimet_s(SD,tind,type,varargin{:});
% load 's.mat'
coords = op_elimtdim(coords); % Fix coords to match

%% Indices for speeds
ind = s(tind,:,:,:) >= sc; 
% ind = data >= sc;

%% Calculate turbine operation time
num = squeeze(sum(ind,1)); % number of indices speed > cut-in speed
denom = size(ind,1); % total number of time indices looked at
top = num/denom;
% top(coords.zm<=0) = nan;
