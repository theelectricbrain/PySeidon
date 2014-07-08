function [sbias,coords] = aimet_sbias(SD,tind,type,varargin)
% Metric sbias, the bias of the mean speed toward ebb or flood tide.
% Can call this function like sbias = met_sbias(u,s) if u and s are
% already read in. Alternatively, call like sbias =
% met_sbias(iind,'i',1,nout) in order to have slices read in in a loop over
% time to save memory.

load 'matfiles/sill1/u.mat'; u=data;
load 'matfiles/sill1/v.mat'; v=data;
% [u,~]=roms_extract(SD,'u',SD.nctime(tind),type,varargin{:});
% [v,coords]=roms_extract(SD,'v',SD.nctime(tind),type,varargin{:});
% %% Resize u and v to be on same grid for calculating s
% % I think we can ignore the non-zslice cases since they are already
% % interpolated
% if ~sum(strcmp(type,{'profile';'point'}))
%     [u,v,coords] = op_uv2psigrid(u,v,coords,tind);
% end
% load 'uvh10.mat' 
% load 'uvsiteprofiles.mat'
% u = u(608:658,:,:); v = v(608:658,:,:); % neap
% u = u(327:377,:,:); v = v(327:377,:,:); % spring

[pa, ~] = principal_axis(u,v);
% save 'pah10.mat' pa coords
% load 'pah10.mat'
% coords = op_elimtdim(coords);    
coords = op_elimtdim(coords);
% ss = aimet_ss(pa,u,v);
% save 'ssh10.mat' ss coords
% load 'ssthalweg.mat'
[ue,ve,uf,vf] = op_uv2ebbflood(u,v,pa);
% [flood,ebb] = met_sbias_op(ss);
sbias = squeeze(nansum(sqrt(ue.^2+ve.^2),1)./nansum(sqrt(uf.^2+vf.^2),1));
    
% Want sbias around 0. sbias = 0 is sflood=sebb.
% sbias < 0 is percentage toward ebb (if multiplied by 100)
% sbias > 0 is percentage toward flood (if multiplied by 100)
sbias(sbias<1) = 1-sbias(sbias<=1); % flood bias
sbias(sbias>1) = -(1-1./sbias(sbias>1)); % ebb bias
end

% function [flood,ebb] = met_sbias_op(ss)
%     flood = squeeze(nansum(ss>0,1));
%     ebb = squeeze(nansum(ss<0,1));
% end