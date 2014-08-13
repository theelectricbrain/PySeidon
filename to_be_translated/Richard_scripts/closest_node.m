function [cidx,varargout]=closest_node(gpsin,nodell,varargin)


	d(:,1)=nodell(:,1)-gpsin(1);
	d(:,2)=nodell(:,2)-gpsin(2);
	d2=sqrt( d(:,1).^2 + d(:,2).^2);
	[closest idx]=sort(d2);	




if length(varargin)==0
	cidx=idx(1);
elseif length(varargin)==1
	number=varargin{1};
	cidx=idx(1:number);
elseif length(varargin)==2
	number=varargin{1};
	h=varargin{2};
	HI= h(:,idx(1:number));
	cidx=idx(1:number);
	varargout{1}=HI;
end
				

	
end

