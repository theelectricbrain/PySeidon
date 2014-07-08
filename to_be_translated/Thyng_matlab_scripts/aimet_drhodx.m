function [dudx,coords] = aimet_drhodx(SD,tind,type,varargin)
% Calculate drho/dx after reading in u and v.
% Trying without central differencing since u and v are on different grids
% anyway.
% Since horizontal distance is in degrees, assume distance of 65 meters

dz = 3;
ds = 65;
grid = roms_get_grid('OUT/ocean_his_0002.nc','OUT/ocean_his_0002.nc');
dxd = (grid.lon_rho(1,2)-grid.lon_rho(1,1)); % delta x in degrees
dyd = (grid.lat_rho(2,1)-grid.lat_rho(1,1)); % delta y in degrees

switch type
%     case {'zslice';'depthslice';'surface'}
%         [u,coords] = roms_extract(SD,'rho',SD.nctime(tind),type,varargin{:});
%         N = ndims(u);
%         if isscalar(tind) %one time step
%             dudx = (u(:,1:end-1)-u(:,2:end))./ds;
%         else % multiple time steps
%             dudx = (u(:,:,1:end-1)-u(:,:,2:end))./ds;
%         end
%         coords = op_resize(coords,N);
    case 'profile' 
        % Size is kx[slice length] or txkx[slice length]
        % Read in neighboring points in necessary directions for
        % derivatives
        y = varargin{1}; x = varargin{2};
        [up,cp] = roms_extract(SD,'rho',SD.nctime(tind),type,y,x+dxd);
        [um,cm] = roms_extract(SD,'rho',SD.nctime(tind),type,y,x-dxd);
        dudx = (up-um)./(cp.xm*65/dxd-cm.xm*65/dxd);
        coords = cp;
        coords.ym = .5*(cp.ym+cm.ym);
        coords.zm = .5*(cp.zm+cm.zm);
%     case 'point'
%         z = varargin{1}; y = varargin{2}; x = varargin{3};
%         [up,coordsup] = roms_extract(SD,'rho',SD.nctime(tind),type,z,y,x+dxd);
%         [um,coordsum] = roms_extract(SD,'rho',SD.nctime(tind),type,z,y,x-dxd);
%         dudx = (um-up)/ds;
%         coords = coordsup;
%         coords.xm = .5*(coordsup.xm+coordsum.xm);
%     case 'full'
%         [u,coords] = roms_extract(SD,'rho',SD.nctime(tind),'full');
%         N = ndims(u);
%         if isscalar(tind) %one time step
%             dudx = (u(:,:,1:end-1)-u(:,:,2:end))./ds;
%         else % multiple time steps
%             dudx = (u(:,:,:,1:end-1)-u(:,:,:,2:end))./ds;
%         end
%         coords = op_resize(coords,N);
end