function coords = op_elimtdim(coords)
% Since coords has extra dimensions for time series, we can use this
% function to eliminate extra entries for time of coords fields. This
% assumes they are all the same. Assumes coords is a structure.

names = fieldnames(coords);

for i=1:length(names)
    field = getfield(coords,names{i});
    coords = setfield(coords,names{i},squeeze(field(1,:,:)));
end