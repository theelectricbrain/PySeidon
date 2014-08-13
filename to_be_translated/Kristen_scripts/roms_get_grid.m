function grd = roms_get_grid(grd_file,scoord,tindex,calc_zuv)
% $Id: roms_get_grid.m 391 2010-03-11 15:04:48Z wilkin $
% grd = roms_get_grid(grd_file)
% grd = roms_get_grid(grd_file,scoord);
% grd = roms_get_grid(grd_file,outfile);
% grd = roms_get_grid(grd_file,outfile,zeta_input,calc_zuv);
%
% Gets the lon,lat,mask,depth [and z coordinates] from netcdf roms netcdf
% grd_file or output file
%
% Input:
%     grd_file: The roms netcdf grid file name
%           or, an existing grd structure to which the vertical coordinates
%               are to be added or updated
%
% Optional inputs:
%     scoord:   ROMS his/rst/avg file from which the s-coord params can be
%               determined
%            or 4-element vector [theta_s theta_b Tcline N] where
%               Vtransform=1 and Vstretching=1 is assumed
%            or 6-element vector [theta_s theta_b Tcline N ...
%                                 Vtransform Vstretching]
%
%     zeta_in:  How to obtain zeta information to use
%               when including free surface height in calculating z-coordinates
%            0 implies assume zeta=0
%            integer implies use this time index into his/rst/avg
%               file to read zeta
%            2-d array of zeta values
%
%     calc_zuv: If present, this argument (any value) activates computing
%               the depths z_u and z_v on the u and v points of the
%               ROMS C-grid
%
% Output is a structure containing all the grid information
%
% NOTE:     All arrays are stored in standard Matlab row-major order
% ******    (C-language) as they appear in netcdf ncdump and
%           NOT in column-major order (Fortran) as assumed by Arango's
%           ROMS Matlab utilities
%           For example,
%           grd.h  (lat,lon)          (j,i)   order
%           grd_z_r(N,lat,lon)        (k,j,i) order
%
% John Wilkin
% Updated (Sept 2002) to correct scoordinate formulation and optionally
% include zeta in the calculation
%
% October 2009:
% Updated (thanks Hernan) to implement new ROMS options for vertical
% coordinate stretching - this change requires the functions
% stretching.m and set_depth.m from Hernan Arango's ROMS Matlab utilities
% available using svn from https://www.myroms.org/svn/src/matlab/utility
% Made do calc_zuv the default

if isstruct(grd_file)
  % if the first input is already a grd structure
  % the intention is to add vertical coordinates below
  grd = grd_file;
  
else
  % get the horizontal grid information from a ROMS grid file (or
  % history, average etc file)
  grd.grd_file = grd_file;
  
  varlist = ...
    { 'mask_rho','mask_psi','mask_u','mask_v','h','pm','pn','f','angle'};
  for v = varlist
    vname = char(v);
    try
      tmp = nc_varget(grd_file,vname);
      grd.(vname) = tmp; % grd = setfield(grd,vname,tmp);
    catch
      %warning('RomsGetGrid:NoMask',['Variable not found: ' vname])
      if strcmp(vname,'angle')
        grd.(vname) = zeros(size(grd.('h')));
      end
    end
  end
  
  varlist = {'x_rho','y_rho','x_u','y_u','x_v','y_v','x_psi','y_psi'};
  if nc_isvar(grd_file,'x_rho')
    for v = varlist
      vname = char(v);
      try
        tmp = nc_varget(grd_file,vname);
        grd.(vname) = tmp; %replaces grd = setfield(grd,vname,tmp);
      catch
        warning('RomsGetGrid:NoXY',['Variable not found: ' vname])
      end
    end
  end
  
  varlist = ...
    { 'lon_rho','lat_rho','lon_psi','lat_psi',...
    'lon_v','lat_v','lon_u','lat_u'};
  for v = varlist
    vname = char(v);
    try
      tmp = nc_varget(grd_file,vname);
      grd.(vname) = tmp; %replaces grd = setfield(grd,vname,tmp);     
    catch
      %warning('RomsGetGrid:NoLonLat',[vname ' not found. Substituting x/y coords instead'])
      if strcmp(vname(1:3),'lon')
        usevname = strrep(vname,'lon','x');
      else
        usevname = strrep(vname,'lat','y');
      end
      grd.(vname) = grd.(usevname);
      %replaces grd = setfield(grd,vname,getfield(grd,usevname));
    end
  end
  
  if isfield(grd,'mask_rho')
    grd.mask_rho_nan = grd.mask_rho;
    land = find(grd.mask_rho_nan==0);
    grd.mask_rho_nan(land) = NaN;
  else
    % there is no mask information in the file so create unit masks in case
    % code tries to use them
    grd.mask_rho = ones(size(grd.h));
    grd.mask_rho_nan = grd.mask_rho;
    grd.mask_u = ones(size(grd.h(:,2:end)));
    grd.mask_v = ones(size(grd.h(2:end,:)));
    grd.mask_psi = ones(size(grd.h(2:end,2:end)));
  end
  
  % If the grid file includes coastline data, such as a file being used
  % with the Rutgers version of editmask.m, load this too
  try
    grd.lon_coast = nc_varget(grd_file,'lon_coast');
  catch
  end
  try
    grd.lat_coast = nc_varget(grd_file,'lat_coast');
  catch
  end
  
end

if nargin > 1
  
  h = grd.h;
  % hmin=min(min(h));
  % hmax=max(max(h));
  [Mp Lp]=size(h);
  L=Lp-1;
  M=Mp-1;
  
  % get z_r and z_w for the given s-coordinate parameters
  
  if ~ischar(scoord)
    
    theta_s = scoord(1);
    theta_b = scoord(2);
    Tcline  = scoord(3);
    N       = scoord(4);
    if (length(scoord) < 5),
      Vtransform = 1;
      Vstretching = 1;
    else
      Vtransform = scoord(5);
      Vstretching = scoord(6);
    end
    hc = Tcline;
    
  else
    
    % input 'scoord' is a his/avg/rst file name or opendap url
    % attempt to get s-coord params from this file/url
    
    theta_s = nc_varget(scoord,'theta_s');
    theta_b = nc_varget(scoord,'theta_b');
    Tcline  = nc_varget(scoord,'Tcline');
    N       = length(nc_varget(scoord,'Cs_r'));
    hc      = nc_varget(scoord,'hc');
    if nc_isvar(scoord,'Vtransform')
      Vtransform = nc_varget(scoord,'Vtransform');
    else
      Vtransform = 1;
    end,
    if nc_isvar(scoord,'Vstretching')
      Vstretching = nc_varget(scoord,'Vstretching');
    else
      Vstretching = 1;
    end
    
  end
  
  if (Vtransform == 1)
    hc = min(Tcline,min(h(:)));
  end
  
  [s_rho,Cs_r]=stretching(Vstretching, theta_s, theta_b, hc, N, 0, 0);
  [s_w  ,Cs_w]=stretching(Vstretching, theta_s, theta_b, hc, N, 1, 0);
  
  % Np = N+1;
  sc_r = s_rho;
  sc_w = s_w;
  
  % zeta
  zeta = zeros(size(grd.h)); % default
  if nargin > 2 % option to include zeta in z calculation
    
    if tindex == 0 % if tindex==0 zeta defaults to zero
      % do nothing
      
    else % if tindex==0 zeta defaults to zero
      if length(tindex)==1
        % tindex is a single index to zeta in a roms output file
        if ~ischar(scoord)
          error([ 'Cannot process zeta from file in the case that ' ...
            ' scoord parameters are input as a vector'])
        end
        zeta = nc_varget(scoord,'zeta',[tindex-1 0 0],[1 -1 -1]);
        if isempty(zeta)
          warning([ 'zeta not found in ' scoord '. Assuming zeta=0.'])
          zeta = zeros(size(grd.h));
        end
      else
        % tindex should be a 2-d field of zeta values
        if any(size(tindex)-size(grd.h))
          % sizes of zeta and h don't match
          error('input tindex as zeta does not match grid dimensions')
        else
          zeta = tindex;
        end
      end
    end
    
  end
  grd.zeta = zeta;
  
  % rho-points
  
  h = grd.h;
  rgrid = 1;
  z_r = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
    rgrid, h', zeta', 0);
  grd.z_r = permute(z_r,[3 2 1]);
  clear z_r
  
  % w-points
  
  wgrid = 5;
  z_w = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
    wgrid, h', zeta', 0);
  grd.z_w = permute(z_w,[3 2 1]);
  
  if nargin > 3
    % if 1
    
    % compute the z depths on the velocity points as well
    % this used to be (before 2009/10/08) optional but there is little
    % reason not to do this and it's annoying when you've forgotten to
    
    % u-points (cell centres in vertical)
    ugrid = 3;
    z_u = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
      ugrid, h', zeta', 0);
    grd.z_u = permute(z_u,[3 2 1]);
    
    % u-points (cell edges in vertical)
    z_uw = 0.5.*(z_w(1:L,1:Mp,:)+z_w(2:Lp,1:Mp,:));
    grd.z_uw = permute(z_uw,[3 2 1]);
    
    % v-points (cell centres in vertical)
    vgrid = 4;
    z_v = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
      vgrid, h', zeta', 0);
    grd.z_v = permute(z_v,[3 2 1]);
    
    % v-points (cell edges in vertical)
    z_vw= 0.5.*(z_w(1:Lp,1:M,:)+z_w(1:Lp,2:Mp,:));
    grd.z_vw = permute(z_vw,[3 2 1]);
    
    clear z_u z_uw z_v z_vw
    
  end
  
  grd.Vtransform = Vtransform;
  grd.Vstretching = Vstretching;
  grd.theta_s = theta_s;
  grd.theta_b = theta_b;
  grd.Tcline = Tcline;
  grd.N = N;
  grd.hc = hc;
  grd.sc_w = sc_w;
  grd.Cs_w = Cs_w;
  grd.sc_r = sc_r;
  grd.Cs_r = Cs_r;
  grd.s_w = s_w;
  grd.s_rho = s_rho;
  
end

%--------------------------------------------------------------------------

function [s,C]=stretching(Vstretching,theta_s,theta_b,hc,N,kgrid,report)

%
% STRETCHING:  Compute ROMS vertical coordinate stretching function
%
% [s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report)
%
% Given vertical terrain-following vertical stretching parameters, this
% this routine computes the vertical stretching function used used in
% ROMS vertical coordinate transformation. Check the following link for
% details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% On Input:
%
%    Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS)
%                    Vstretching = 3,  R. Geyer BBL refinement
%    theta_s       S-coordinate surface control parameter (scalar)
%    theta_b       S-coordinate bottom control parameter (scalar)
%    hc            Width (m) of surface or bottom boundary layer in which
%                    higher vertical resolution is required during
%                    stretching (scalar)
%    N             Number of vertical levels (scalar)
%    kgrid         Depth grid type logical switch:
%                    kgrid = 0,        function at vertical RHO-points
%                    kgrid = 1,        function at vertical W-points
%    report        Flag to report detailed information (OPTIONAL):
%                    report = 0,       do not report
%                    report = 1,       report information
%
% On Output:
%
%    s             S-coordinate independent variable, [-1 <= s <= 0] at
%                    vertical RHO- or W-points (vector)
%    C             Nondimensional, monotonic, vertical stretching function,
%                    C(s), 1D array, [-1 <= C(s) <= 0]
%

% svn $Id: roms_get_grid.m 391 2010-03-11 15:04:48Z wilkin $
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

s=[];
C=[];

%----------------------------------------------------------------------------
%  Set several parameters.
%----------------------------------------------------------------------------

if (nargin < 7)
  report=0;
end

%- Np=N+1;

%----------------------------------------------------------------------------
% Compute ROMS S-coordinates vertical stretching function
%----------------------------------------------------------------------------

switch Vstretching
  
  case 1
    % Original vertical stretching function (Song and Haidvogel, 1994)
    
    ds=1.0/N;
    if (kgrid == 1)
      %-    Nlev=Np;
      lev=0:N;
      s=(lev-N).*ds;
    else
      %-    Nlev=N;
      lev=[1:N]-0.5;
      s=(lev-N).*ds;
    end
    if (theta_s > 0)
      Ptheta=sinh(theta_s.*s)./sinh(theta_s);
      Rtheta=tanh(theta_s.*(s+0.5))./(2.0*tanh(0.5*theta_s))-0.5;
      C=(1.0-theta_b).*Ptheta+theta_b.*Rtheta;
    else
      C=s;
    end
    
  case 2
    % A. Shchepetkin (UCLA-ROMS) vertical stretching function.
    
    alfa=1.0;
    beta=1.0;
    ds=1.0/N;
    if (kgrid == 1)
      %-    Nlev=Np;
      lev=0:N;
      s=(lev-N).*ds;
    else
      %-    Nlev=N;
      lev=[1:N]-0.5;
      s=(lev-N).*ds;
    end
    if (theta_s > 0)
      Csur=(1.0-cosh(theta_s.*s))/(cosh(theta_s)-1.0);
      if (theta_b > 0)
        Cbot=-1.0+sinh(theta_b*(s+1.0))/sinh(theta_b);
        weight=(s+1.0).^alfa.*(1.0+(alfa/beta).*(1.0-(s+1.0).^beta));
        C=weight.*Csur+(1.0-weight).*Cbot;
      else
        C=Csur;
      end
    else
      C=s;
    end
    
  case 3
    %  R. Geyer BBL vertical stretching function.
    
    ds=1.0/N;
    if (kgrid == 1)
      %-    Nlev=Np;
      lev=0:N;
      s=(lev-N).*ds;
    else
      %-    Nlev=N;
      lev=[1:N]-0.5;
      s=(lev-N).*ds;
    end
    if (theta_s > 0)
      exp_s=theta_s;      %  surface stretching exponent
      exp_b=theta_b;      %  bottom  stretching exponent
      alpha=3;            %  scale factor for all hyperbolic functions
      Cbot=log(cosh(alpha*(s+1).^exp_b))/log(cosh(alpha))-1;
      Csur=-log(cosh(alpha*abs(s).^exp_s))/log(cosh(alpha));
      weight=(1-tanh( alpha*(s+.5)))/2;
      C=weight.*Cbot+(1-weight).*Csur;
    else
      C=s;
    end
    
  otherwise
    error(['STRETCHING - Illegal parameter Vstretching = ' ...
      num2str(Vstretching)])
    
end

% Report S-coordinate parameters.

if (report)
  warning('ROMS_GET_GRID:noreport',...
    'Function STRETCHING embedded in roms_get_grid does not report')
end

%--------------------------------------------------------------------------

function [z]=set_depth(Vtransform,Vstretching, ...
  theta_s,theta_b,hc,N,igrid,h,zeta,report)
%
% SET_DEPTH:  Compute ROMS grid depth from vertical stretched variables
%
% [z]=set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
%               igrid, h, zeta);
%
% Given a batymetry (h), free-surface (zeta) and terrain-following parameters,
% this function computes the 3D depths for the request C-grid location. If the
% free-surface is not provided, a zero value is assumef resulting in unperturb
% depths.  This function can be used when generating initial conditions or
% climatology data for an application. Check the following link for details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% On Input:
%
%    Vtransform    Vertical transformation equation:
%                    Vtransform = 1,   original transformation
%
%                      z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]
%
%                      Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)
%
%                    Vtransform = 2,   new transformation
%
%                      z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)
%
%                       Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]
%    Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS)
%                    Vstretching = 3,  R. Geyer BBL refinement
%    theta_s       S-coordinate surface control parameter (scalar)
%    theta_b       S-coordinate bottom control parameter (scalar)
%    hc            Width (m) of surface or bottom boundary layer in which
%                    higher vertical resolution is required during
%                    stretching (scalar)
%    N             Number of vertical levels (scalar)
%    igrid         Staggered grid C-type (integer):
%                    igrid=1  => density points
%                    igrid=2  => streamfunction points
%                    igrid=3  => u-velocity points
%                    igrid=4  => v-velocity points
%                    igrid=5  => w-velocity points
%    h             Bottom depth, 2D array at RHO-points (m, positive),
%                    h(1:Lp+1,1:Mp+1).
%    zeta          Free-surface, 2D array at RHO-points (m), OPTIONAL,
%                    zeta(1:Lp+1,1:Mp+1).
%    report        Flag to report detailed information (OPTIONAL):
%                    report = 0,       do not report
%                    report = 1,       report information
%
% On Output:
%
%    z             Depths (m, negative), 3D array.
%

% svn $Id: roms_get_grid.m 391 2010-03-11 15:04:48Z wilkin $
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

z=[];

%----------------------------------------------------------------------------
%  Set several parameters.
%----------------------------------------------------------------------------

if (Vtransform < 1 || Vtransform > 2)
  error(['SET_DEPTH Illegal param Vtransform = ' num2str(Vtransfrom)])
end

if (Vstretching < 1 || Vstretching > 3)
  error(['SET_DEPTH Illegal param Vstretching = ' num2str(Vstretching)])
end

if (hc > min(min(h)) && Vtransform == 1)
  disp(['Vtranform = ' num2str(Vtransform)])
  disp(['hc        = ' num2str(hc)])
  disp(['hmax      = ' num2str(min(h(:)))])
  disp('SET_DEPTH - critical depth exceeds minimum bathymetry value.')
end

if (nargin < 9)
  zeta=zeros(size(h));
end

if (nargin < 10)
  report=1;
end

Np=N+1;
[Lp Mp]=size(h);
L=Lp-1;
M=Mp-1;

%- hmin=min(min(h));
%- hmax=max(max(h));

%----------------------------------------------------------------------------
% Compute vertical stretching function, C(k):
%----------------------------------------------------------------------------

if (report)
  warning('ROMS_GET_GRID:noreport',...
    'Function SET_DEPTH embedded in roms_get_grid does not report')
end

if (igrid == 5)
  kgrid=1;
else
  kgrid=0;
end

[s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report);

%----------------------------------------------------------------------------
%  Average bathymetry and free-surface at requested C-grid type.
%----------------------------------------------------------------------------

switch ( igrid )
  case 1
    hr=h;
    zetar=zeta;
  case 2
    hp=0.25.*(h(1:L,1:M)+h(2:Lp,1:M)+h(1:L,2:Mp)+h(2:Lp,2:Mp));
    zetap=0.25.*(zeta(1:L,1:M)+zeta(2:Lp,1:M)+zeta(1:L,2:Mp)+zeta(2:Lp,2:Mp));
  case 3
    hu=0.5.*(h(1:L,1:Mp)+h(2:Lp,1:Mp));
    zetau=0.5.*(zeta(1:L,1:Mp)+zeta(2:Lp,1:Mp));
  case 4
    hv=0.5.*(h(1:Lp,1:M)+h(1:Lp,2:Mp));
    zetav=0.5.*(zeta(1:Lp,1:M)+zeta(1:Lp,2:Mp));
  case 5
    hr=h;
    zetar=zeta;
end

%----------------------------------------------------------------------------
% Compute depths (m) at requested C-grid location.
%----------------------------------------------------------------------------

if (Vtransform == 1)
  switch ( igrid )
    case 1
      for k=1:N
        z0=(s(k)-C(k))*hc + C(k).*hr;
        z(:,:,k)=z0 + zetar.*(1.0 + z0./hr);
      end
    case 2
      for k=1:N
        z0=(s(k)-C(k))*hc + C(k).*hp;
        z(:,:,k)=z0 + zetap.*(1.0 + z0./hp);
      end
    case 3
      for k=1:N
        z0=(s(k)-C(k))*hc + C(k).*hu;
        z(:,:,k)=z0 + zetau.*(1.0 + z0./hu);
      end,
    case 4
      for k=1:N
        z0=(s(k)-C(k))*hc + C(k).*hv;
        z(:,:,k)=z0 + zetav.*(1.0 + z0./hv);
      end,
    case 5
      z(:,:,1)=-hr;
      for k=2:Np
        z0=(s(k)-C(k))*hc + C(k).*hr;
        z(:,:,k)=z0 + zetar.*(1.0 + z0./hr);
      end
  end
elseif (Vtransform == 2)
  switch ( igrid )
    case 1
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hr)./(hc+hr);
        z(:,:,k)=zetar+(zeta+hr).*z0;
      end
    case 2
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hp)./(hc+hp);
        z(:,:,k)=zetap+(zetap+hp).*z0;
      end
    case 3
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hu)./(hc+hu);
        z(:,:,k)=zetau+(zetau+hu).*z0;
      end
    case 4
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hv)./(hc+hv);
        z(:,:,k)=zetav+(zetav+hv).*z0;
      end
    case 5
      for k=1:Np,
        z0=(hc.*s(k)+C(k).*hr)./(hc+hr);
        z(:,:,k)=zetar+(zetar+hr).*z0;
      end
  end
end