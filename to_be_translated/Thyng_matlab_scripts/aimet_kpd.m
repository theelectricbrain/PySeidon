function [kpd,datacoords] = aimet_kpd(SD,tind,type,varargin)
% Calculate kinetic power density after calculating speed in aimet_s.
% Assume constant value for rho in this case. rho=1024 kg/m^3.

%% For ai100
%% Constant density value
rho = 1024; % kg/m^3

%% Calculate speed
load 'matfiles/sz50.mat'
% [s,coords] = aimet_s(SD,tind,type,varargin{:});

%% Calculate kinetic power density
kpd = .5*rho*data.^3/1000;

% %% For ai65 and for mkpd
% kpd = 0; ms = 0; mrho = 0;
% ksave=zeros(1,length(tind)); mssave=zeros(1,length(tind)); mrhosave=zeros(1,length(tind));
% for i=1:length(tind)
%     %% Constant density value
%     [rho,~]=roms_extract(SD,'rho',SD.nctime(tind(i)),type,varargin{:});
%     rho = op_resize(op_resize(rho,1),2);
%     
%     %% Calculate speed
%     [s,coords] = aimet_s(SD,tind(i),type,varargin{:});
% 
%     %% Calculate kinetic power density
%     kpd = kpd + .5*(rho+1000).*s.^3/1000;
%     ksave(i) = kpd(297,121)/i; % just divide by i since dt is the same every step
%     
%     %% Save mean speed and mean density while I am at it
%     ms = ms + s;
%     mssave(i) = s(297,121)/i;
%     mrho = mrho + rho;
%     mrhosave(i) = rho(297,121)/i;
% end
% kpd = kpd/length(tind);
% save 'mkpdh10.mat' coords kpd ksave ms mssave mrho mrhosave