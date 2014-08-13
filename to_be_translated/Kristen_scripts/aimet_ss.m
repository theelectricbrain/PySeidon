function ss=aimet_ss(pa,u,v)
% function ss=met_ss(pa,u,v,xlims,grid)
% Calculate signed speed from principal axis values.
% This gives the wrong direction for flows in eddy fields.

pa = repmat(pa,[1 1 length(u(:,1,1))]); pa=permute(pa,[3 1 2]);

% % Eliminate entries with speed < cut-in speed
%  s=sqrt(u.^2+v.^2);
%  u(s<.7)=nan;
%  v(s<.7)=nan;
% If pa is nan, eliminate that velocity point
u(isnan(pa))=nan; v(isnan(pa))=nan;

% Rotate from principal axis to normal axes to identify points as ebb or flood
% Convert back to cartesian coords from compass coords
up=u.*cosd(-(-pa+90))-v.*sind(-(-pa+90));
% vp=u.*sind(-(-pa+90))+v.*cosd(-(-pa+90));

ue=u; uf=u; ve=v; vf=v; clear u v;
% Fill in opposite spots with nan's.
% If pa<0, do one, if pa>0, do the other
ue(logical((up<0).*(pa<0)))=nan; uf(logical((up>0).*(pa<0)))=nan;
ve(logical((up<0).*(pa<0)))=nan; vf(logical((up>0).*(pa<0)))=nan;
ue(logical((up>0).*(pa>0)))=nan; uf(logical((up<0).*(pa>0)))=nan;
ve(logical((up>0).*(pa>0)))=nan; vf(logical((up<0).*(pa>0)))=nan;

% ss = cosd(-pa+90).*(u.*(up>0))+sind(-pa+90).*(v.*(up>0)) + ... 
%     cosd(-pa+90).*(u.*(up<0))+sind(-pa+90).*(v.*(up<0));
% Project onto principal axis
ssf = cosd(-pa+90).*uf+sind(-pa+90).*vf;
sse = cosd(-pa+90).*ue+sind(-pa+90).*ve;
if ndims(sse) == 1
    ss = nansum(cat(2,ssf,sse),2);
elseif ndims(sse) == 2
    ss = nansum(cat(3,ssf,sse),3);
elseif ndims(sse) == 3
    ss = nansum(cat(4,ssf,sse),4);
end

% The following is incorrect
% sse = -sqrt(ue.^2+ve.^2);
% ssf = sqrt(uf.^2+vf.^2);
% ss.var = sse+ssf;