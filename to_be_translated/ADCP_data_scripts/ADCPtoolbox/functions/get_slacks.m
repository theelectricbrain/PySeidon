%Brian Polagye
%April 3, 2010

%Description: identify slack waters in time series

%Change Log:
%   -6/9/10: modified to compute slacks at multiple depths

function [slack_ind, slack_t] = get_slacks(s, t, T_thresh)

%Inputs:
%   s: signed speed series (m/s)
%   t: times (matlab serial format) associated with each speed
%   T_thresh: the minimum actual tidal period (in hours) - oscillations around slack
%       water can lead to spuriously short periods

%Outputs:
%   slack_ind: indices associated with sign change
%   slack_t: interpolated times of slack water

for j = 1:size(s, 2)
    
    bin(j).slack_ind = find(sign(s(1:end-1,j))-sign(s(2:end,j))~=0);
    
    %eliminate false slack identifications caused by NaN values
    NaN_error = intersect(bin(j).slack_ind, find(isnan(s(:,j))));
    if ~isempty(NaN_error)
        NaN_error = [NaN_error; NaN_error-1];
        %error trap for NaN value at start of series
        if NaN_error(1) == 1;
            NaN_error = setxor(NaN_error,0);
        end
    end
    
    bin(j).slack_ind = setxor(bin(j).slack_ind, NaN_error);   %prune NaN errors from slack indices
    
    %prune artificial slack at t=1
    if bin(j).slack_ind(1) == 1, bin(j).slack_ind = bin(j).slack_ind(2:end); end
    %prune artificial slack at t=end
    if bin(j).slack_ind(end) == size(t,1), bin(j).slack_ind = bin(j).slack_ind(1:end-1); end     
    
    bin(j).slack_t = NaN(size(bin(j).slack_ind));
    
    for i = 1:length(bin(j).slack_ind)  
        %if any interpolation values are NaN or are exactly repeated, do not try to interpolate
        %either of these cases will throw an error from interp1
        if isnan(s(bin(j).slack_ind(i)-1,j)) || isnan(s(bin(j).slack_ind(i),j)) || ...
                isnan(s(bin(j).slack_ind(i)+1,j)) || length(unique([s(bin(j).slack_ind(i)-1,j), ...
                s(bin(j).slack_ind(i),j),s(bin(j).slack_ind(i)+1,j)])) < 3
            bin(j).slack_t(i) = t(bin(j).slack_ind(i));
        %standardard three point interpolation to identify time of slack
        else
            bin(j).slack_t(i) = interp1([s(bin(j).slack_ind(i)-1,j),s(bin(j).slack_ind(i),j),s(bin(j).slack_ind(i)+1,j)],...
                [t(bin(j).slack_ind(i)-1),t(bin(j).slack_ind(i)),t(bin(j).slack_ind(i)+1)],0);
        end
    end
    
    %eliminate "slacks" which are nearly concurrent and result from oscillations around slack water
    %note: this does truncate the time series by one slack
    T = (bin(j).slack_t(2:end)-bin(j).slack_t(1:end-1))*24;   %duration of slack water
    bin(j).slack_ind = bin(j).slack_ind(T>T_thresh);
    bin(j).slack_t = bin(j).slack_t(T>T_thresh);
    
end

%create a uniform slack structure for all bins
for i = 1:length(bin)
   slack_count(i) = length(bin(i).slack_ind);
end
slack_count = min(slack_count);

slack_ind = zeros(slack_count, length(bin));
slack_t = zeros(slack_count, length(bin));

for i = 1:length(bin)
    slack_ind(:,i) = bin(i).slack_ind(1:slack_count);
    slack_t(:,i) = bin(i).slack_t(1:slack_count);
end