%Description: routine to generate signed speed (flood positive) from time
%   series
% flood_heading is the angle (clockwise from North) of the flood tide. 
% Outputs the signed speed and the principal axis at each vertical level. 
% Brian Polagye
% June 14, 2010
%
% Modifications:
% July 20th -- modified to output PA for each vertical level and to 
%               convert the flood_heading to a vector within the code. 
function [s_signed_all PA_all] = sign_speed(u_all, v_all, s_all, dir_all, flood_heading)

if length(flood_heading)==1
    flood_heading = flood_heading + [-90, +90];
end

%initialize signed magnitude
s_signed_all = NaN(size(s_all));

%generate signed speed
for i = 1:size(s_all,2)

    %determine signed speed
    u = u_all(:,i);     %u-component
    v = v_all(:,i);     %v-component
    dir = dir_all(:,i); %direction
    s = s_all(:,i);     %speed

    %determine principal axes - potentially a problem if axes are very kinked
    %   since this would misclassify part of ebb and flood
    [PA ~] = principal_axis(u, v);
    PA_all(i) = PA;

    %sign speed - eliminating wrap-around
    dir_PA = dir - PA;
  

    dir_PA(dir_PA<-90) = dir_PA(dir_PA<-90) + 360;
    dir_PA(dir_PA>270) = dir_PA(dir_PA>270) - 360;

    %general direction of flood passed as input argument
    if PA >=flood_heading(1) && PA <=flood_heading(2)
        ind_fld = find(dir_PA >= -90 & dir_PA<90);
        s_signed = -s;
        s_signed(ind_fld) = s(ind_fld);
    else
        ind_ebb = find(dir_PA >= -90 & dir_PA<90);
        s_signed = s;
        s_signed(ind_ebb) = -s(ind_ebb);
    end
    
    s_signed_all(:,i) = s_signed;

end

end

%% principal_axis: determine principal axis from directional scatter

function [PA, varxp_PA] = principal_axis(u,v)

U = [u, v];                             %create velocity matrix
U = U(~isnan(U(:,1)),:);                %eliminate NaN values
U = U - repmat(mean(U,1),length(U),1);  %convert matrix to deviate form
R = U'*U/(length(U)-1);                 %compute covariance matrix (alternatively - cov(U))
[V,lambda]=eig(R);                      %calculate eigenvalues and eigenvectors for covariance matrix

%sort eignvalues in descending order so that major axis is given by first eigenvector
[lambda, ilambda]=sort(diag(lambda),'descend');     %sort in descending order with indices
lambda=diag(lambda);                                %reconstruct the eigenvalue matrix
V=V(:,ilambda);                                     %reorder the eigenvectors

ra = atan2(V(2,1),V(2,2));   %rotation angle of major axis in radians relative to cartesian coordiantes

PA = -ra*180/pi+90;                         %express principal axis in compass coordinates
varxp_PA = diag(lambda(1))/trace(lambda);   %variance captured by principal axis

end