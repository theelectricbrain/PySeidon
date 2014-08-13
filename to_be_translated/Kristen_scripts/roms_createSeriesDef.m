function SD = roms_createSeriesDef(dirname, basename)

% SD = roms_createSeriesDef(dirname, basename);
%
% makes a simple structure containing the input arguments and the timebase
% of the corresponding .nc files
%
% neil banas feb 2009

if dirname(end) ~= '/', dirname = [dirname '/']; end % make sure there's a trailing slash
SD.dirname = dirname;
SD.basename = basename;

[nctime,ncn] = roms_outputdir2timebase(dirname,basename);
SD.ncn = ncn;
SD.nctime = nctime;

