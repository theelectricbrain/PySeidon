function [theta,coords] = aimet_theta(SD,tind,type,varargin)
% Calculate angle between vertical velocity and horizontal speed.

%% Constants
dz = 0.5;

%% Angle
switch type
    case {'zslice';'depthslice';'hslice'}
        [s,coords] = aimet_s(SD,tind,type,varargin{:});
        [w,coords] = roms_extract(SD,'w',SD.nctime(tind),type,varargin{:});
        theta = atand(w./s);
        % Resize
        N = ndims(coordsvu.xm);
    case 'profile' 
        % Size is kx[slice length] or txkx[slice length]
        load 'matfiles/sill1/u.mat'; u=data;
        load 'matfiles/sill1/v.mat'; v=data;
        load 'matfiles/sill1/w.mat'; w=data;
        s = sqrt(u.^2+v.^2);
        s = op_resize(s,2);
%         [u,coords] = roms_extract(SD,'u',SD.nctime(tind),type,varargin{:});
%         [v,~] = roms_extract(SD,'v',SD.nctime(tind),type,varargin{:});
        theta = atand(w./s);
%     case 'point'
%         z = varargin{1}; y = varargin{2}; x = varargin{3};
%         % Read in z-level of u and v above and below to calculate z derivative
%         % uu,ud,vu,vd are either jxi or txjxi
%         [uu,coordsuu] = roms_extract(SD,'u',SD.nctime(tind),type,z+dz,y,x); %shallowest
%         [ud,coordsud] = roms_extract(SD,'u',SD.nctime(tind),type,z-dz,y,x); %deepest
%         shearu = (uu-ud)./(coordsuu.zm-coordsud.zm);
%         [vu,coordsvu] = roms_extract(SD,'v',SD.nctime(tind),type,z+dz,y,x); %shallowest
%         [vd,coordsvd] = roms_extract(SD,'v',SD.nctime(tind),type,z-dz,y,x); %deepest
%         shearv = (vu-vd)./(coordsvu.zm-coordsvd.zm);
%         coords = coordsvd;
%         coords.zm = .5*(coordsvu.zm+coordsvd.zm);
%         shear = sqrt(shearu.^2+shearv.^2);
end