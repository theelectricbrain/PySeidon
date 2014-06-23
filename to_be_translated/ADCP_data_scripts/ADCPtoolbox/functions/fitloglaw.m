%Brian Polagye
%November 7, 2011

%Description: routine to calculate best fit to a log law for a vertical
%   profile of the form u(z) = u_shear/K*ln(z/z0)

function out = fitloglaw(guess,u_data,z_data)

    K = 0.4;    %Von Karman's constant

    %initialize guess - this will change during minimization
    a = guess(1);   %shear velocity
    b = guess(2);   %bottom roughness
    
    %minimization target - difference between observation and profile model
    out = sum((u_data - a/K*log(z_data/b)).^2); 
    
end
