function data = op_taverage(data,SD,tind)
% Do average over time period

% if not starting set of read-ins from 1, need to shift tind
if tind(1) ~= 1
    tind = 1:length(tind);
end

% N = length(tind);
% 
% dt = SD.nctime(2)-SD.nctime(1); % in days
% Perform average over input time indices
% data = squeeze(nansum(data,1)*dt/(SD.nctime(tind(end))-SD.nctime(tind(1))));
% Some nan's come in at points that are in the water and so need to not
% divide by them for average.
% data = squeeze(nansum(data,1)./sum(~isnan(data),1));
% In case we are taking a subset of a full already-read-in set
if ndims(data==3)
    data = squeeze(nansum(data(tind,:,:),1)./sum(~isnan(data(tind,:,:)),1));
elseif ndims(data==4)
    data = squeeze(nansum(data(tind,:,:,:),1)./sum(~isnan(data(tind,:,:,:)),1));
end
% data = squeeze(nansum(data,1)/N);