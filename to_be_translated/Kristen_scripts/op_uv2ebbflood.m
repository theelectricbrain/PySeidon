function [ue,ve,uf,vf] = op_uv2ebbflood(u,v,pa)
% Takes in u, v and principal axis angle to split into uebb, vebb, uflood,
% and vflood.
% This is an approximation at telling direction of tides but may not do a
% great job.

% load 'flooddir.mat' %in compass coordinates
% As filler just make one up since not doing direction anymore
flooddir = ones(size(pa));
for j=1:size(pa,1)
    for i=1:size(pa,2)
        if ~isnan(pa(j,i))
            temp = op_closestAngles([flooddir(j,i) pa(j,i)]);
            pa(j,i) = temp(2);
        end
    end
end

%% compare here before pa changes size
if abs(pa-flooddir) < 90 % then pa points in flood dir and up rotates toward flood tide
    dir = 1; % pa is toward flood
else % pa points in ebb direction
    dir = 0; % pa is toward ebb
end
clear flooddir

pa = repmat(pa,[1 1 length(u(:,1,1))]); pa=permute(pa,[3 1 2]);

% % Eliminate entries with speed < cut-in speed
%  s=sqrt(u.^2+v.^2);
%  u(s<.7)=nan;
%  v(s<.7)=nan;
% If pa is nan, eliminate that velocity point
u(isnan(pa))=nan; v(isnan(pa))=nan;

save matfiles/vtemp.mat v
clear v
%% Rotate from principal axis to normal axes to identify points as ebb or flood
% Convert back to cartesian coords from compass coords
up1=u.*cosd(-(-pa+90));%-v.*sind(-(-pa+90));
save matfiles/utemp.mat u up1
clear u up1
load matfiles/vtemp.mat
v=-v.*sind(-(-pa+90));
% save matfiles/utemp2.mat up2
clear pa
load matfiles/utemp.mat up1
up = up1+v;
clear up1 v
load matfiles/vtemp.mat
% up=u.*cosd(-(-pa+90))-v.*sind(-(-pa+90));
% vp=u.*sind(-(-pa+90))+v.*cosd(-(-pa+90));
ve=v; vf=v; 
clear v
load matfiles/utemp.mat
ue=u; uf=u; 
clear u;

%% Compare pa angle with predetermined angle in flooddir.mat that points
% approximately in flood direction. If difference is within 90 degrees,
% after running angles through op_closestAngles to be sure they are on the
% same ring, then pa is already pointing in flood dir. Otherwise pa is
% pointing in ebb direction. Then we nan out entries in ue, ve, uf, and vf
% accordingly.
if dir % then pa points in flood dir and up rotates toward flood tide
    ue(up>0) = nan;
    ve(up>0) = nan;
    uf(up<0) = nan;
    vf(up<0) = nan;
else % pa points in ebb direction
    ue(up<0) = nan;
    ve(up<0) = nan;
    uf(up>0) = nan;
    vf(up>0) = nan;
end

% % Fill in opposite spots with nan's.
% % If pa<0, do one, if pa>0, do the other. First change pa back to compass
% % coordinates for this check since atan2 gives answers for compass coords
% % between -180 and +180
% % ue(logical((up<0).*(-pa+90<0)))=nan; uf(logical((up>0).*(-pa+90<0)))=nan;
% % ve(logical((up<0).*(-pa+90<0)))=nan; vf(logical((up>0).*(-pa+90<0)))=nan;
% ue(logical((up>0).*(-pa+90>0)))=nan; uf(logical((up<0).*(-pa+90>0)))=nan;
% ve(logical((up>0).*(-pa+90>0)))=nan; vf(logical((up<0).*(-pa+90>0)))=nan;