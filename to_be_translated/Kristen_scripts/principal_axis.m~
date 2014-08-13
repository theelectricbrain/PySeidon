%Brian Polagye
%May 7, 2010
%[PA, varxp_PA] = principal_axis(u,v)
%Description: revised version of principal axis code
%
%Input:
%   - x-velocity component (u)
%   - y-velocity component (v)
%
%Output:
%   - PA: principal axis direction in compass coordinates (degrees)
%   - varxp_PA: variance explained by principal axis

function [PA, varxp_PA] = principal_axis(u,v)

%     s=sqrt(u.^2+v.^2);
%     u(s<.7)=nan;
%     v(s<.7)=nan;

debug = 0;%1;

% keep code as brian had it if u is a vector
if isvector(u)
    [PA, varxp_PA] = pa(u,v);
elseif ndims(u) == 2 %&& ~isvector(u)
    PA = nan(length(u(1,:)),1);
    varxp_PA = nan(length(u(1,:)),1);
    for i = 1:length(u(1,:))
        if sum(~isnan((squeeze(u(:,i)))))>4 && sum(~isnan((squeeze(v(:,i)))))>4
            [PA(i), varxp_PA(i)] = pa(squeeze(u(:,i)),squeeze(v(:,i)));
        end
    end
% if u is a matrix, do differently to speed up, assume time x j x i
% Assume that first dimension is the time dimension
elseif ndims(u) == 3 %&& ~isvector(u)
    PA = nan(length(u(1,:,1)),length(u(1,1,:)));
    varxp_PA = nan(length(u(1,:,1)),length(u(1,1,:)));
    for j = 1:length(u(1,:,1))
        for i = 1:length(u(1,1,:))
            if sum(~isnan((squeeze(u(:,j,i)))))>4 && sum(~isnan((squeeze(v(:,j,i)))))>4
                [PA(j,i), varxp_PA(j,i)] = pa(squeeze(u(:,j,i)),squeeze(v(:,j,i)));
            end
        end
    end
% time x k x j x i    
elseif ndims(u) == 4    
    PA = zeros(length(u(1,:,1,1)),length(u(1,1,:,1)),length(u(1,1,1,:)));
    varxp_PA = zeros(length(u(1,:,1,1)),length(u(1,1,:,1)),length(u(1,1,1,:)));
    for k = 1:length(u(1,:,1,1))
        for j = 1:length(u(1,1,:,1))
            for i = 1:length(u(1,1,1,:))
                [PA(k,j,i), varxp_PA(k,j,i)] = pa(squeeze(u(:,k,j,i)),squeeze(v(:,k,j,i)));
            end
        end
    end
    
end


end

function [PA, varxp_PA] = pa(u,v)
    U = [u, v];                             %create velocity matrix
    ind = logical(~isnan(U(:,1)).*~isnan(U(:,2)));
    U = U(ind,:);                %eliminate NaN values
    if size(U(:,1))==1 %only one point
        PA=nan;
        varxp_PA=nan;
        return
    else
    U = U - repmat(mean(U,1),length(U),1);  %convert matrix to deviate form
    R = U'*U/(length(U)-1);                 %compute covariance matrix (alternatively - cov(U))
    [V,lambda]=eig(R);                      %calculate eigenvalues and eigenvectors for covariance matrix

    %sort eignvalues in descending order so that major axis is given by first eigenvector
    [lambda, ilambda]=sort(diag(lambda),'descend');     %sort in descending order with indices
    lambda=diag(lambda);                                %reconstruct the eigenvalue matrix
    V=V(:,ilambda);                                     %reorder the eigenvectors

    ra = atan2(V(2,1),V(2,2));   %rotation angle of major axis in radians relative to cartesian coordiantes
% disp('not currently in compass coordinates')
% PA = ra*180/pi;
    PA = -ra*180/pi+90;                         %express principal axis in compass coordinates
    varxp_PA = diag(lambda(1))/trace(lambda);   %variance captured by principal axis
    end
end