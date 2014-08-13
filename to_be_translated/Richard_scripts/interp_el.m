function [HI]=interp_el(XI,YI,h,trinodes,uvnodell,nodell,aw0,awx,awy,varargin)

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
					n1  = trinodes(i,1);
					n2  = trinodes(i,2);
					n3  = trinodes(i,3);	
					x0c=XI-uvnodell(i,1);
			 		y0c=YI-uvnodell(i,2);				
					h0 = aw0(i,1)*h(n1,:)+aw0(i,2)*h(n2,:)+aw0(i,3)*h(n3,:);
					hx = awx(i,1)*h(n1,:)+awx(i,2)*h(n2,:)+awx(i,3)*h(n3,:);
					hy = awy(i,1)*h(n1,:)+awy(i,2)*h(n2,:)+awy(i,3)*h(n3,:);
					HI = h0 + hx*x0c + hy*y0c;
				else
					HI=NaN;   				
				end

