function [shear,coords] = aimet_shear(SD,tind,type,varargin)
% Calculate shear after reading in u and v.
% Shearu = (u_{up}-u_{down})./(z_{up}-z_{down})
% Shearv = (v_{up}-v_{down})./(z_{up}-z_{down})
% Shear = sqrt(u.^2+v.^2)

%% Constants
% dz = 0.5;
dz=2;
%% Shear
switch type
    case {'zslice';'depthslice';'hslice'}
        % Read in z-level of u and v above and below to calculate z derivative
        % uu,ud,vu,vd are either jxi or txjxi
        [uu,coordsuu] = roms_extract(SD,'u',SD.nctime(tind),type,varargin{1}+dz); %shallowest
%         save 'matfiles/uuhslice10.mat' uu coordsuu
%         clear uu coordsuu
        [ud,coordsud] = roms_extract(SD,'u',SD.nctime(tind),type,varargin{1}-dz); %deepest
%         save 'matfiles/udhslice10.mat' ud coordsud
%         clear ud coordsud
        [vu,coordsvu] = roms_extract(SD,'v',SD.nctime(tind),type,varargin{1}+dz); %shallowest
%         save 'matfiles/vuhslice10.mat' vu coordsvu
%         clear vu coordsvu
        [vd,coordsvd] = roms_extract(SD,'v',SD.nctime(tind),type,varargin{1}-dz); %deepest
%         save 'matfiles/vdhslice10.mat' vd coordsvd
%         clear vd coordsvd
%         load 'matfiles/uuhslice10.mat'; load 'matfiles/udhslice10.mat'; 
        shearu = (uu-ud)./(coordsuu.zm-coordsud.zm);
        clear uu coordsuu ud coordsud
%         load 'matfiles/vuhslice10.mat'; load 'matfiles/vdhslice10.mat';
        shearv = (vu-vd)./(coordsvu.zm-coordsvd.zm);
        clear vu vd 
        % Resize
        N = ndims(coordsvu.xm);
        shearu = op_resize(shearu,N-1);
        shearv = op_resize(shearv,N);
        coordsvu.zm = .5*(coordsvu.zm+coordsvd.zm);
        coords = op_resize(coordsvu,N);
        shear = sqrt(shearu.^2+shearv.^2);
    case 'profile' 
        % Size is kx[slice length] or txkx[slice length]
%         load matfiles/uvjslice.mat
%         load 'matfiles/usave.mat'; u=data;
%         load 'matfiles/vsave.mat'; v=data; coords=datacoords;
        [u,coords] = roms_extract(SD,'u',SD.nctime(tind),type,varargin{:});
        [v,~] = roms_extract(SD,'v',SD.nctime(tind),type,varargin{:});
        if isscalar(tind) % one time only
            shearu = (u(1:end-2,:)-u(3:end,:))./(coords.zm(1:end-2,:)-coords.zm(3:end,:));
            shearv = (v(1:end-2,:)-v(3:end,:))./(coords.zm(1:end-2,:)-coords.zm(3:end,:));
            coords.xm = coords.xm(2:end-1,:);
            coords.ym = coords.ym(2:end-1,:);
            coords.zm = coords.zm(2:end-1,:);
        else % multiple times
            shearu = (u(:,1:end-2,:)-u(:,3:end,:))./(coords.zm(:,1:end-2,:)-coords.zm(:,3:end,:));
            shearv = (v(:,1:end-2,:)-v(:,3:end,:))./(coords.zm(:,1:end-2,:)-coords.zm(:,3:end,:));
            coords.xm = coords.xm(:,2:end-1,:);
            coords.ym = coords.ym(:,2:end-1,:);
            coords.zm = coords.zm(:,2:end-1,:);
        end
        shear = sqrt(shearu.^2+shearv.^2);
%         save 'matfiles/sill1/shear.mat' shear coords
    case 'point'
        z = varargin{1}; y = varargin{2}; x = varargin{3};
        % Read in z-level of u and v above and below to calculate z derivative
        % uu,ud,vu,vd are either jxi or txjxi
        [uu,coordsuu] = roms_extract(SD,'u',SD.nctime(tind),type,z+dz,y,x); %shallowest
        [ud,coordsud] = roms_extract(SD,'u',SD.nctime(tind),type,z-dz,y,x); %deepest
        shearu = (uu-ud)./(coordsuu.zm-coordsud.zm);
        [vu,coordsvu] = roms_extract(SD,'v',SD.nctime(tind),type,z+dz,y,x); %shallowest
        [vd,coordsvd] = roms_extract(SD,'v',SD.nctime(tind),type,z-dz,y,x); %deepest
        shearv = (vu-vd)./(coordsvu.zm-coordsvd.zm);
        coords = coordsvd;
        coords.zm = .5*(coordsvu.zm+coordsvd.zm);
        shear = sqrt(shearu.^2+shearv.^2);
end