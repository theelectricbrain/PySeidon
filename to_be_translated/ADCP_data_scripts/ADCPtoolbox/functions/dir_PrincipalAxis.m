%Brian Polagye
%November 11, 2011

%Description: principal axis decomposition by vertical bin
%   0 degrees = N, CW positive
%   The tricky aspect is to account for is the quadrants represented by
%   flood and ebb with a 0/360 wrap around.

function [PA_fld, PA_ebb, PA_all, d_PA] = dir_PrincipalAxis(d, s, min_speed, u, v)

%inputs
% - z: vertical bins to process
% - d: current direction (t,z)
% - s: current speed (t,z)
% - min_speed: minimum speed to include for analysis
% - u: x-component of horizontal velocity (t,z)
% - v: y-component of horizontal velocity (t,z)

%Outputs
% - PA_fld, PA_ebb, PA_all: principal axes for flood, ebb, and composite tide
% - d_PA: current direction mapped into principal axis coordinates (no 0/360 wrap)

PA_fld = zeros(size(s,2),1);
PA_ebb = zeros(size(s,2),1);
PA_all = zeros(size(s,2),1);
d_PA = zeros(size(d));


for i = 1:size(s,2)     %loop through all depth bins
    
    ind_real = find(abs(s(:,i))>=min_speed);  %work only with speeds greater than min

    %determine principal axes - potentially problematic if axes are very kinked
    U = [u(ind_real,i), v(ind_real,i)];     %create velocity matrix
    U = U(~isnan(U(:,1)),:);                %eliminate NaN values
    U = U - repmat(mean(U,1),length(U),1);  %convert matrix to deviate form
    R = U'*U/(length(U)-1);                 %compute covariance matrix (alternatively - cov(U))
    [V,lambda]=eig(R);                      %calculate eigenvalues and eigenvectors for covariance matrix

    %sort eignvalues in descending order so that major axis is given by first eigenvector
    [~, ilambda]=sort(diag(lambda),'descend');     %sort in descending order with indices
    V=V(:,ilambda);                                %reorder the eigenvectors

    ra = atan2(V(2,1),V(2,2));      %rotation angle of major axis in radians relative to cartesian coordiantes
    PA_all(i) = -ra*180/pi+90;      %express principal axis in compass coordinates

    %create matrix of signed speeds - eliminate wrap-around 360 degrees
    d_PA(:,i) = d(:,i) - PA_all(i);

    d_PA(d_PA(:,i)<-90,i) = d_PA(d_PA(:,i)<-90,i) + 360;
    d_PA(d_PA(:,i)>270,i) = d_PA(d_PA(:,i)>270,i) - 360;

    %determine principal axis for flood and ebb, accounting for wrapping
    %   Here, the principle axes are treated as the mean direction for ebb
    %   or flood tides. This allows for directional assymetry. Ebb and
    %   flood directions have been previously assigned using similar
    %   algorithm.
    PA_fld(i) = mean(d_PA(s(:,i)>=min_speed,i)) + PA_all(i);
    PA_ebb(i) = mean(d_PA(s(:,i)<=-min_speed,i)) + PA_all(i);
end

end
