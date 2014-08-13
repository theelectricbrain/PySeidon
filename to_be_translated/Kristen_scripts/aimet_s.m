function [s,coords] = aimet_s(SD,tind,type,varargin)
% Calculate speed after reading in u and v

% load 'uvthalweg.mat'
% load 'uvh10.mat'
% load 'uvslice1.mat'
load 'matfiles/uz50.mat'; u=data; clear data
load 'matfiles/vz50.mat'; v=data; clear data

% Read in u and v
% [u,~]=roms_extract(SD,'u',SD.nctime(tind),type,varargin{:});
% [v,coords]=roms_extract(SD,'v',SD.nctime(tind),type,varargin{:});
% save matfiles/uvz50_13to57.mat u v coords
% grid=roms_get_grid('OUT/ocean_his_0002.nc','OUT/ocean_his_0002.nc');
% u=zeros(length(tind),size(grid.lat_u,1),size(grid.lat_u,2));
% v=zeros(length(tind),size(grid.lat_v,1),size(grid.lat_v,2));
% for i=1:length(tind)
%     [u(i,:,:),~]=roms_extract(SD,'u',SD.nctime(tind(i)),type,varargin{:});
%     [v(i,:,:),coords]=roms_extract(SD,'v',SD.nctime(tind(i)),type,varargin{:});
% end
% % save 'uh10.mat' u coords
% % save 'vh10.mat' v coords
% 
% load matfiles/uvbar.mat
% u=ubar; v=vbar;
% load matfiles/uvz50_13to57resized.mat %troc

%% Resize u and v to be on same grid for calculating s
% I think we can ignore the non-zslice cases since they are already
% interpolated
if ~sum(strcmp(type,{'profile';'point'}))
    [u,v,coords] = op_uv2psigrid(u,v,datacoords);
%     [u,v,coords] = op_uv2psigrid(u,v,coords);
end
save 'matfiles/uvz50resized.mat' u v coords
% load matfiles/uvsurfaceresized.mat
% load 'uvsiteprofiles.mat'
% save 'uvslice7.mat' u v coords
% save 'uh10.mat' u coords
% save 'vh10.mat' v coords
s = sqrt(u.^2+v.^2);