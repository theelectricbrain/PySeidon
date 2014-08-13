function newdata = op_resize(data,dim)
% Resize data down in dimension dim by averaging. If input is a structure,
% perform resizing on each field of structure.

if isstruct(data)
    names = fieldnames(data);
    newdata = [];
   
    for i=1:length(names)
        field = getfield(data,names{i});
        N = ndims(field);

        if N == 1
            newdata = setfield(newdata,names{i},.5*(field(1:end-1)+field(2:end)));
        elseif N == 2
            if dim == 1
                newdata = setfield(newdata,names{i},.5*(field(1:end-1,:)+field(2:end,:)));
            elseif dim == 2
                newdata = setfield(newdata,names{i},.5*(field(:,1:end-1)+field(:,2:end)));
            end
        elseif N == 3
            if dim == 1
                newdata = setfield(newdata,names{i},.5*(field(1:end-1,:,:)+field(2:end,:,:)));
            elseif dim == 2
                newdata = setfield(newdata,names{i},.5*(field(:,1:end-1,:)+field(:,2:end,:)));
            elseif dim == 3
                newdata = setfield(newdata,names{i},.5*(field(:,:,1:end-1)+field(:,:,2:end)));
            end
        elseif N == 4
            if dim == 1
                newdata = setfield(newdata,names{i},.5*(field(1:end-1,:,:,:)+field(2:end,:,:,:)));
            elseif dim == 2
                newdata = setfield(newdata,names{i},.5*(field(:,1:end-1,:,:)+field(:,2:end,:,:)));
            elseif dim == 3
                newdata = setfield(newdata,names{i},.5*(field(:,:,1:end-1,:)+field(:,:,2:end,:)));
            elseif dim == 4
                newdata = setfield(newdata,names{i},.5*(field(:,:,:,1:end-1)+field(:,:,:,2:end)));
            end
        end
    end
else
    N = ndims(data);

    if N == 1
        newdata = .5*(data(1:end-1)+data(2:end));
    elseif N == 2
        if dim == 1
            newdata = .5*(data(1:end-1,:)+data(2:end,:));
        elseif dim == 2
            newdata = .5*(data(:,1:end-1)+data(:,2:end));
        end
    elseif N == 3
        if dim == 1
            newdata = .5*(data(1:end-1,:,:)+data(2:end,:,:));
        elseif dim == 2
            newdata = .5*(data(:,1:end-1,:)+data(:,2:end,:));
        elseif dim == 3
            newdata = .5*(data(:,:,1:end-1)+data(:,:,2:end));
        end
    elseif N == 4
        if dim == 1
            newdata = .5*(data(1:end-1,:,:,:)+data(2:end,:,:,:));
        elseif dim == 2
            newdata = .5*(data(:,1:end-1,:,:)+data(:,2:end,:,:));
        elseif dim == 3
            newdata = .5*(data(:,:,1:end-1,:)+data(:,:,2:end,:));
        elseif dim == 4
            newdata = .5*(data(:,:,:,1:end-1)+data(:,:,:,2:end));
        end
    end
end