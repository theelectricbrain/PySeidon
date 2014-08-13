%Brian Polagye
%April 14, 2010

%Change log:
%   - 6/29/11 - modified to calculate number of elements internally
%   - 6/29/11 - modified to calculate ensemble along a specified dimension

%
%Description: function to calculate ensemble average of measurement data
%   can process nx1 and nxm dimension matrices

%Input:
%   x: time or spatial series to be ensembled
%   ens: elements per ensemble
%   ens_dim: dimension to take ensemble average along

function [x_ens] = calc_ensemble(x, ens, ens_dim)
    %initialize input
    if ens_dim == 1
        ens_size = [floor(size(x,1)/ens), size(x,2)];
    else
        ens_size = [size(x,1), floor(size(x,2)/ens)];
    end
    
    x_ens = NaN([ens_size,ens]);
    
    %create matrix of ensemble values
    for j = 1:ens
        if ens_dim == 1
            ind_ens = j:ens:size(x,1)-(ens-j);
            x_ens(:,:,j) = x(ind_ens,:);
        else
            ind_ens = j:ens:size(x,2)-(ens-j);
            x_ens(:,:,j) = x(:,ind_ens);
        end
    end
    
    %calculate ensemble
    x_ens = nanmean(x_ens,3);   %average along third dimension
    
end