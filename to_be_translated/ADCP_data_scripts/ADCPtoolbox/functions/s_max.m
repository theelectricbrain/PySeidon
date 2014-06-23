% Brian Polagye
% November 11, 2011

function [s_max_fld, s_max_ebb, s_max_all] = s_max(s)

%Description: find maximum speed

%Inputs:
%   - s: horizontal velocity

%Outputs:
%   - s_max_fld: maximum flood velocity
%   - s_max_ebb: maximum ebb velocity
%   - s_max_all: maximum velocity over all tides

s_max_ebb = min(s);
s_max_fld = max(s);
s_max_all = max(abs(s));

end
