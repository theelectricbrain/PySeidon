function phi = op_anglesFromVectors(u,v)
% This function takes in vectors in the form (u,v) and compares them in
% order to find the angles of the vectors without any wrap-around issues.
% This is accomplished by finding the smallest difference between angles
% compared at different wrap-around values.
% This appears to work correctly.
% Angles are reported in compass coordinates.
% This same function is done when input angles in op_closestAngles.

%% Find initial angle values for all vectors in compass coordinates
phi = -atan2(v,u)*180/pi+90;

if ndims(u)==2
    for j = 1:size(u,2) % go through each column separately
        % Find first non-nan to use as comparison angle
        ind = find(~isnan(phi(:,j)),1);
        % If all entries in this column are nan's, move to next column
        if isempty(ind)
            continue
        end
        for i = ind+1:size(u,1) % Start loop at index after first non-nan
            % If phi(i,j) is nan, skip rest of loop to save time
            if isnan(phi(i,j))
                continue
            end
            % Test-difference 1: initial difference between angles
            diff1 = abs(phi(ind,j)-phi(i,j));
            % Test-difference 2: difference between angles when phi(i) is moved
            % down a ring
            diff2 = abs(phi(ind,j)-(phi(i,j)-360));
            % Test-difference 3: difference between angles with phi(i) is moved up
            % a ring
            diff3 = abs(phi(ind,j)-(phi(i,j)+360));
            if diff1 < diff2 && diff1 < diff3 % phi(i) stays the same
                phi(i,j) = phi(i,j);
            elseif diff2 < diff1 && diff2 < diff3 % phi(i) moves down a ring
                phi(i,j) = phi(i,j)-360;
            elseif diff3 < diff1 && diff3 < diff2 % phi(i) moves up a ring
                phi(i,j) = phi(i,j)+360;
            end
        end
    end
elseif ndims(u)==3
    for k = 1:size(u,3)
        for j = 1:size(u,2) % go through each column separately
            % Find first non-nan to use as comparison angle
            ind = find(~isnan(phi(:,j,k)),1);
            % If all entries in this column are nan's, move to next column
            if isempty(ind)
                continue
            end
            for i = ind+1:size(u,1) % Start loop at index after first non-nan
                % If phi(i,j) is nan, skip rest of loop to save time
                if isnan(phi(i,j,k))
                    continue
                end
                % Test-difference 1: initial difference between angles
                diff1 = abs(phi(ind,j,k)-phi(i,j,k));
                % Test-difference 2: difference between angles when phi(i) is moved
                % down a ring
                diff2 = abs(phi(ind,j,k)-(phi(i,j,k)-360));
                % Test-difference 3: difference between angles with phi(i) is moved up
                % a ring
                diff3 = abs(phi(ind,j,k)-(phi(i,j,k)+360));
                if diff1 < diff2 && diff1 < diff3 % phi(i) stays the same
                    phi(i,j,k) = phi(i,j,k);
                elseif diff2 < diff1 && diff2 < diff3 % phi(i) moves down a ring
                    phi(i,j,k) = phi(i,j,k)-360;
                elseif diff3 < diff1 && diff3 < diff2 % phi(i) moves up a ring
                    phi(i,j,k) = phi(i,j,k)+360;
                end
            end
        end
    end
end