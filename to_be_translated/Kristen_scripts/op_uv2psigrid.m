function [u,v,coords] = op_uv2psigrid(u,v,coords)
% Move u and v from their respective grids to a psi grid, along with v's
% coords.

% N = ndims(u);

% if isscalar(tind) %&& ~sum(strcmp(type,{'profile';'point'}))% One time step
    u = op_resize(u,ndims(u)-1);
    v = op_resize(v,ndims(v));
    coords = op_resize(coords,ndims(coords.xm));
% else % Multiple time steps
%     if ~isvector(u)
%         u = op_resize(u,ndims(u)-1); 
%         v = op_resize(v,ndims(v));
%         coords = op_resize(coords,ndims(coords));
%     end
% end
