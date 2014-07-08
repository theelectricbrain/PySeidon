function [vort,coords] = aimet_vort(SD,tind,type,varargin)
% Calculate vertical vorticity after reading in u and v.
% Trying without central differencing since u and v are on different grids
% anyway.
% vort = dvdx - dudy
% Since horizontal distance is in degrees, assume distance of 100 meters
% (current grid size horizont

%% Constants
dz = 3;
ds = 65;
% dx = 30; dy = 10; % troc
dx = 65; dy = 65; % ai65
grid = roms_get_grid('OUT/ocean_his_0002.nc','OUT/ocean_his_0002.nc');
dxd = (grid.lon_rho(1,2)-grid.lon_rho(1,1)); %  delta x in degrees
dyd = (grid.lat_rho(2,1)-grid.lat_rho(1,1)); %  delta y in degrees
% dxd = .5*(grid.lon_rho(1,2)-grid.lon_rho(1,1)); % .5 delta x in degrees
% dyd = .5*(grid.lat_rho(2,1)-grid.lat_rho(1,1)); % .5 delta y in degrees

%% Vorticity
switch type
    case {'zslice';'depthslice';'surface'}
        [u,~] = roms_extract(SD,'u',SD.nctime(tind),type,varargin{:});
        [v,coords] = roms_extract(SD,'v',SD.nctime(tind),type,varargin{:});
        u = op_resize(u,1);
        v = op_resize(v,2);
        coords = op_resize(coords,2);
%         load matfiles/uvz50resized.mat
        [vx,~] = gradient(v,dx);
        [~,uy] = gradient(u,dy);
        vort = vx-uy; % vertical vorticity: vx-uy    
%         N = ndims(u);
%         if isscalar(tind) %one time step
%             dvdx = (v(:,1:end-1)-v(:,2:end))./dx;
%             dudy = (u(1:end-1,:)-u(2:end,:))./dy;
%         else % multiple time steps
%             dvdx = (v(:,:,1:end-1)-v(:,:,2:end))./dx;
%             dudy = (u(:,1:end-1,:)-u(:,2:end,:))./dy;
%         end
%         clear u v
%         vort = dvdx - dudy;
    case 'profile' 
        % Size is kx[slice length] or txkx[slice length]
        % Read in neighboring points in necessary directions for
        % derivatives
        y = varargin{1}; x = varargin{2};
        [up,coordsup] = roms_extract(SD,'u',SD.nctime(tind),type,y+dyd,x);
        [um,coordsum] = roms_extract(SD,'u',SD.nctime(tind),type,y-dyd,x);
        [vp,coordsvp] = roms_extract(SD,'v',SD.nctime(tind),type,y,x+dxd);
        [vm,coordsvm] = roms_extract(SD,'v',SD.nctime(tind),type,y,x-dxd);
        dvdx = (vp-vm)./(coordsvp.xm*65/dxd-coordsvm.xm*65/dxd);%(2*ds);
        dudy = (up-um)./(coordsup.ym*65/dyd-coordsum.ym*65/dyd);%(2*ds);
        vort = dvdx - dudy;
        coords = coordsup;
        coords.ym = .5*(coordsup.ym+coordsum.ym);
        coords.zm = .5*(coordsup.zm+coordsum.zm);
    case 'point'
        z = varargin{1}; y = varargin{2}; x = varargin{3};
        [up,coordsup] = roms_extract(SD,'u',SD.nctime(tind),type,z,y+dyd,x);
        [um,coordsum] = roms_extract(SD,'u',SD.nctime(tind),type,z,y-dyd,x);
        [vp,coordsvp] = roms_extract(SD,'v',SD.nctime(tind),type,z,y,x+dxd);
        [vm,coordsvm] = roms_extract(SD,'v',SD.nctime(tind),type,z,y,x-dxd);
        dvdx = (vm-vp)./(coordsvp.xm*65/dxd-coordsvm.xm*65/dxd);%ds;
        dudy = (um-up)./(coordsup.ym*65/dyd-coordsum.ym*65/dyd);
        vort = dvdx - dudy;
        coords = coordsup;
        coords.ym = .5*(coordsup.ym+coordsum.ym);
    case 'full'
        [u,~] = roms_extract(SD,'u',SD.nctime(tind),'full');
        [v,coords] = roms_extract(SD,'v',SD.nctime(tind),'full');
        N = ndims(u);
        if isscalar(tind) %one time step
            dvdx = (v(:,:,1:end-1)-v(:,:,2:end))./ds;
            dudy = (u(:,1:end-1,:)-u(:,2:end,:))./ds;
        else % multiple time steps
            dvdx = (v(:,:,:,1:end-1)-v(:,:,:,2:end))./ds;
            dudy = (u(:,:,1:end-1,:)-u(:,:,2:end,:))./ds;
        end
        vort = dvdx - dudy; clear dvdx dudy
        coords = op_resize(coords,N);
end