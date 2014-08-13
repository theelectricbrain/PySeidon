function [cidx,varargout]=closest_element(gpsin,uvnodell,varargin)


	d(:,1)=uvnodell(:,1)-gpsin(1);
	d(:,2)=uvnodell(:,2)-gpsin(2);
	d2=sqrt( d(:,1).^2 + d(:,2).^2);
	[closest idx]=sort(d2);	


if length(varargin)==0
	cidx=idx(1);
elseif length(varargin)==1
	number=varargin{1};
	cidx=idx(1:number);
elseif length(varargin)==3
	number=varargin{1};
	u=varargin{2};
	v=varargin{3};
	UI= u(:,idx(1:number));
	VI= v(:,idx(1:number));
	cidx=idx(1:number);
	varargout{1}=UI;
	varargout{2}=VI;	
end	

	

end

