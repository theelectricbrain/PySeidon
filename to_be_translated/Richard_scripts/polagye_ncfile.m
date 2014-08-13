clear

%% User input

%File names
%path=['/array/data1/073208o/workspace_matlab/runs/dg_2014/output/'];
path=['/EcoII/july_2012/output/'];
savepath=['/array/data1/rkarsten/'];
%    '13hr_run/output/'];
file='dngrid';

%Hydrodynamics metrics parameters
min_speed=0.2;%Sets the divide between slack and non-slack

%User sets coordinates of E-W and N-S boundaries of 
%rectangular region of interest 
%Particular code below is unnecessarily complicated;
%just set west, south, east and north

%CLA coordinates
%GP	tight	
		region=[-66.355 -66.32 44.245 44.2925];
filesave='GP_data.csv';
filesavemat='GP_data.mat';
%PP
%region=[-66.225 -66.195 44.37 44.41];
%filesave='PP_data.csv';
%DG
% region=[-65.79 -65.73 44.65 44.7];
% filesave='DG_data.csv';

west=region(1);
south=region(3);
east=region(2);
north=region(4);

%Expand CLA region to desired rectangular region
%by adding X metres to each side
% [westx,southy]=gmt_proj(west,south,'forward',...
%     'utm20t',[]);
% [eastx,northy]=gmt_proj(east,north,'forward',...
%     'utm20t',[]);
% %westx=westx-2000;%-CLAexpand;
% %eastx=eastx+4000;%+CLAexpand+1800;
% %southy=southy-1500;%-CLAexpand;
% %northy=northy+1250;%+CLAexpand;
% %westx=westx-0;%7500;%-CLAexpand;
% %eastx=eastx+6000;%+CLAexpand+1800;
% %southy=southy-5000;%-CLAexpand;
% %northy=northy+1500;%+CLAexpand;
% [west,south]=gmt_proj(westx,southy,'inverse',...
%     'utm20t',[]);
% [east,north]=gmt_proj(eastx,northy,'inverse',...
%     'utm20t',[]);

%Broader view of Minas Passage
%west=-64.707662;
%south=45.315916;
%east=-64.324625;
%north=45.401058;
%
delta=0;%0.0025;
ax=[west-delta east+delta south-delta north+delta];


%% Load nc variables 

ncfile=[path,file,'_0001_02.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
lontri = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonc'));
lattri = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latc'));
xc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xc'));
yc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'yc'));
x = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x'));
y = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y'));
t = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
trinodes = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nv'));
%u = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ua'));
%v = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'va'));
%el = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'));
h = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'));

lon=double(lon);
lat=double(lat);
%u=double(u);
%v=double(v);
h=double(h);
%el=double(el);
t=double(t);

lonlat(:,1)=lon;
lonlat(:,2)=lat;

%x, y and el interpolated to element centres
% lontri=(lon(trinodes(:,1))+lon(trinodes(:,2))+...
%         lon(trinodes(:,3)))/3;
% lattri=(lat(trinodes(:,1))+lat(trinodes(:,2))+...
%         lat(trinodes(:,3)))/3;
htri=(h(trinodes(:,1))+h(trinodes(:,2))+...
        h(trinodes(:,3)))/3; 
%eltri=(el(trinodes(:,1),:)+el(trinodes(:,2),:)+...
%        el(trinodes(:,3),:))/3;  


%% Use only those points in rectangular region of interest
%Note ua's spatial extent is limited to rectangular region 
%of interest before loading (full) va 
%%{
FORCE=(lontri>west & lontri<east & lattri>south & lattri<north);
lontri=lontri(FORCE);
lattri=lattri(FORCE);
xc=xc(FORCE);
yc=yc(FORCE);
htri=htri(FORCE);
%u = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ua'));
u = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 2 0],[length(trinodes) 2 length(t)]);
u=double(squeeze(mean(u,2))');
u=u(:,FORCE);
%v = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'va'));
v = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 2 0],[length(trinodes) 2 length(t)]);
v=double(squeeze(mean(v,2))');
v=v(:,FORCE);
%}

%switch u and v for flow that is almost north south

%% Computer hydrodynamics metrics

%Compute signed speed and direction
magvel = (u.^2 + v.^2).^(1/2);
%meanmagvel=(mean(magvel,1))';

dirvel = atan2(v, u); 
dirvel = (pi/2)-dirvel;
%neg_dirvel = dirvel<0;
%nonneg_dirvel = dirvel>=0;
%dirvel = (2*pi+dirvel).*(dirvel<0) + dirvel.*(dirvel>=0);
d = dirvel*180/pi;

clear *dirvel*

s = sign_speed_ncfile(u, v, magvel, d); 
%

%%Velocity characteristics
% - mean
% - max velocity
% - asymmetry

[s_fld, s_ebb, s_all] = s_mean_ncfile(min_speed, s);
[smax_fld, smax_ebb, smax_all] = s_max_ncfile(s);
smax_fld=smax_fld';smax_ebb=smax_ebb';smax_all=smax_all';
s_asym = abs(s_ebb./s_fld);

%%Directional characteristics

% - principal axes
% - standard deviation
% - directional asymmetry

[PA_fld, PA_ebb, PA_more , d_PA] = ...
    dir_PrincipalAxis_ncfile(d, s, min_speed, u, v);
[std_fld, std_ebb, std_all] = ...
    dir_std_ncfile(s, min_speed, d_PA);
dir_asym=abs(abs((PA_fld-PA_ebb))-180);

 std_fld(isnan(std_fld))=0;
 std_ebb(isnan(std_ebb))=0;
% 
 dir_asym(isnan(dir_asym))=0;

%Power characteristics
min_power=0;
[P_fld, P_ebb, P_all] = p_mean_ncfile(min_power, s);
P_asym = abs(P_ebb./P_fld);
P_asym(P_asym>10)=0;
P_asym(isnan(P_asym))=0;

%% Write to text file and terminate code ('stop')

%%{
combined_metrics=[lontri, lattri,smax_all, P_all, P_asym,...
    std_fld, std_ebb, dir_asym]; 
dlmwrite([savepath,filesave],combined_metrics,...
    'delimiter',',','precision','%.5f');
index=FORCE;
save([savepath,filesavemat],'lontri', 'lattri','smax_all', 'P_all', 'P_fld','P_ebb','P_asym','std_fld', 'std_ebb', 'dir_asym','lonlat','trinodes','h','index','x','y','xc','yc')
%fid=fopen('ambient_flow_metrics.csv','w+');
%fprintf(fid,'%2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f',...
%    combined_metrics);

%stop
%}

%% For visualizing region of interest (uncomment 'Write to..'
% above if want to visualize

for ii=3:8
fullgrid=zeros(size(FORCE,1),1);
fullgrid(FORCE)=combined_metrics(:,ii);

%Write output into existing NetCDF file

%ncid = netcdf.open([path,file,'_0001.nc'], 'NC_WRITE');           
%netcdf.reDef(ncid);
%dimid = netcdf.inqDimID(ncid,'nele');
%dir_asymmetry = netcdf.defVar(ncid,'dir_asymmetry','float',dimid);
%netcdf.putAtt(ncid,dir_asymmetry,'long_name','dir_asymmetry');
%netcdf.endDef(ncid);
%netcdf.putVar(ncid,dir_asymmetry,dir_asym);
%netcdf.close(ncid);

figure
%hold on
%plot([west west], [south north], 'Color', [0.7 0.7 0.7],...
%        'LineWidth', 2);
patch('Vertices',lonlat,'Faces',trinodes,'Cdata',fullgrid)
shading flat
%caxis([0 max(fullgrid)])
%caxis([0 1])
ch=colorbar;
set(ch,'fontsize',20)
axis(ax)
axis off
end
