function [UI,VI]=interp_vel(XI,YI,u,v,trinodes,nbe,uvnodell,nodell,a1u,a2u,varargin)

if length(varargin)==1
	host=varargin{1};

else	
	for i=1:length(trinodes)
		if inpolygon(XI,YI,nodell(trinodes(i,:),1),nodell(trinodes(i,:),2))==1
			host=i;
			break;
		end
	end
end
















			i=host;
			if (i~=0)
				e1=nbe(i,1);
				e2=nbe(i,2);
			 	e3=nbe(i,3);
				e1(e1==0)=length(trinodes)+1;
				e2(e2==0)=length(trinodes)+1;
				e3(e3==0)=length(trinodes)+1;
			 	x0c=XI-uvnodell(i,1);
			 	y0c=YI-uvnodell(i,2);
				dudx= a1u(i,1)*u(i,:)+a1u(i,2)*u(e1,:)+a1u(i,3)*u(e2,:)+a1u(i,4)*u(e3,:);
				dudy= a2u(i,1)*u(i,:)+a2u(i,2)*u(e1,:)+a2u(i,3)*u(e2,:)+a2u(i,4)*u(e3,:);
				dvdx= a1u(i,1)*v(i,:)+a1u(i,2)*v(e1,:)+a1u(i,3)*v(e2,:)+a1u(i,4)*v(e3,:);
				dvdy= a2u(i,1)*v(i,:)+a2u(i,2)*v(e1,:)+a2u(i,3)*v(e2,:)+a2u(i,4)*v(e3,:);
				UI= u(i,:) + dudx*x0c + dudy*y0c;
				VI= v(i,:) + dvdx*x0c + dvdy*y0c;
			else
				UI=NaN;
   				VI=NaN; 
			end
