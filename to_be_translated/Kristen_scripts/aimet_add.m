function [a,std]=aimet_add(u,v,pa)%,coords)
% Calculates asymmetry parameter a and standard deviation of angles on 
% each direction
% It assumes u is time x j x i or time x k x i or time x space x space
% Call as [a,stdf,stde]=met_add(u,v,pa) if have already read in/found u, v,
% and pa.
% Call as [a,stdf,stde]=met_add(iind,'i',N,nout) if need to save memory and
% read in in chunks of time.
% Meane and meanf are in compass coords
% clear coords
[ue,ve,uf,vf] = op_uv2ebbflood(u,v,pa);
save matfiles/uvef.mat ue ve uf vf
clear u v pa
% pa = repmat(pa,[1 1 length(u(:,1,1))]); pa=permute(pa,[3 1 2]);
% 
% % % Eliminate entries with speed < cut-in speed
% %  s=sqrt(u.^2+v.^2);
% %  u(s<.7)=nan;
% %  v(s<.7)=nan;
% % If pa is nan, eliminate that velocity point
% u(isnan(pa))=nan; v(isnan(pa))=nan;
% 
% % Rotate from principal axis to normal axes to identify points as ebb or flood
% % Convert back to cartesian coords from compass coords
% up=u.*cosd(-(-pa+90))-v.*sind(-(-pa+90));
% % vp=u.*sind(-(-pa+90))+v.*cosd(-(-pa+90));
% 
% ue=u; uf=u; ve=v; vf=v; clear u v;
% % Fill in opposite spots with nan's.
% % If pa<0, do one, if pa>0, do the other
% ue(logical((up<0).*(pa<0)))=nan; uf(logical((up>0).*(pa<0)))=nan;
% ve(logical((up<0).*(pa<0)))=nan; vf(logical((up>0).*(pa<0)))=nan;
% ue(logical((up>0).*(pa>0)))=nan; uf(logical((up<0).*(pa>0)))=nan;
% ve(logical((up>0).*(pa>0)))=nan; vf(logical((up<0).*(pa>0)))=nan;
% clear pa

% Get angles from vectors
% save uvefh10sc07.mat ue ve uf vf coords
phie = op_anglesFromVectors(ue,ve);
phif = op_anglesFromVectors(uf,vf);
% clear ue ve vf uf
save matfiles/phi.mat phie phif
% Calculate a, stde, stdf, std using cubic weighting with speed
s = sqrt(ue.^2+ve.^2); clear ue ve
w = abs(s.^3)./repmat(nansum(abs(s.^3)),[size(s,1) 1 1]); % Cubic weighting of speed
meane = squeeze(nansum(phie.*w)); % Weighted average with w normalized
stde = squeeze(sqrt(nansum(w.*(phie-repmat(reshape(meane,[1 size(meane)]),[size(phie,1) 1 1])).^2))); % Weighted standard deviation
s = sqrt(uf.^2+vf.^2); clear uf vf
w = abs(s.^3)./repmat(nansum(abs(s.^3)),[size(s,1) 1 1]); % Linear weighting of speed
meanf = squeeze(nansum(phif.*w)); % in compass coords
stdf = squeeze(sqrt(nansum(w.*(phif-repmat(reshape(meanf,[1 size(meanf)]),[size(phif,1) 1 1])).^2))); % Weighted standard deviation
a = abs(meane-meanf-180);
a(a>180)=abs(a(a>180)-360);
std=max(stde,stdf);
load matfiles/uvz50resized.mat 'coords'
% load matfiles/add.mat
% load matfiles/uvz50resized.mat 'coords'
ind = (a==180);
a(ind) = nan;
ind = (std==0);
std(ind) = nan;
save 'matfiles/add.mat' a std coords

% save 'matfiles/meanadcpvel.mat' meane meanf
% save 'matfiles/stdadcpvel.mat' stde stdf
% 
% save phih10sc07.mat phie phif coords
% save meanh10sc0w.mat meane meanf coords
% save stdh10sc0w.mat stde stdf coords
% save ah10sc0w.mat a coords
%     
%% Plot
% x = op_slicedist(coords.xm*100/0.00134529147982221,coords.ym*100/0.000901162790697185);
% x = squeeze(x(1,:,:)); y = squeeze(coords.zm(1,:,:));
color = lbmap(9,'Turqoise');
x = squeeze(coords.xm(1,:,:))/1000; y = squeeze(coords.ym(1,:,:))/1000;
aimet_plot(x,y,a,color,'fname','add/az50','clims',[0 40],...
    'xlabel','km','ylabel','km','zlabel','Bidirectionality (deg)')
aimet_plot(x,y,a,color,'fname','add/az50zoom1','clims',[0 40],...
    'xlabel','km','ylabel','km','zlabel','Bidirectionality (deg)','xlims',[10 30])
% std=max(stde,stdf);
aimet_plot(x,y,std,color,'fname','add/stdz50','clims',[0 40],...
    'xlabel','km','ylabel','km','zlabel','Directional Deviation (deg)')
aimet_plot(x,y,std,color,'fname','add/stdz50zoom1','clims',[0 40],...
    'xlabel','km','ylabel','km','zlabel','Directional Deviation (deg)','xlims',[10 30])
% % % x = op_slicedist(coords.xm*100/0.00134529147982221,coords.ym*100/0.000901162790697185);
% % % z = squeeze(coords.zm(1,:,:)); x = squeeze(x(1,:,:));
% % % % x = squeeze(coords.xm(1,:,:)); y = squeeze(coords.ym(1,:,:));
% % % fnamea = 'add/aslice1spring';
% % % aimet_plot(x,z,a,color,'fname',fnamea), hold on
% % % fnamestd = 'add/stdslice1spring';
% % % aimet_plot(x,z,std,color,'fname',fnamestd), hold on
% 
% end