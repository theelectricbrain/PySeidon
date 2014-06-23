%Brian Polagye
%September 14, 2011

%Description: routine to calculate best fit to a power law for a vertical
%   profile of the form u(z) = u_0(z/z_0)^(1/alpha)

function out = fitpowerlaw(guess,u_data,z_data)

    %note: perhaps would also be interesting to search based on index which
    %minimizes, but unclear if that would ever end due to discretization
    %error

    %initialize guess - this will change during minimization
    a = guess(1);
    b = guess(2);
    
    %minimization target - difference between observation and profile model
    out = sum((u_data - a*(z_data).^b).^2); 
    
end
