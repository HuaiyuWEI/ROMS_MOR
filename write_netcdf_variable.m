function write_netcdf_variable(filename, variable_name, data, verify_readback)
%WRITE_NETCDF_VARIABLE Write a variable with NCWRITE and read it back.
%   WRITE_NETCDF_VARIABLE(FILENAME, VARIABLE_NAME, DATA) writes DATA, reads
%   the variable back and compares it with DATA cast to the type NCREAD
%   returns, so a double written to an NC_FLOAT variable passes. An error
%   is raised when the two differ.
%   WRITE_NETCDF_VARIABLE(..., false) skips the read-back.

    if nargin < 4 || isempty(verify_readback)
        verify_readback = true;
    end
    validateattributes(verify_readback, {'logical', 'numeric'}, ...
        {'scalar'}, mfilename, 'verify_readback');

    filename = char(filename);
    variable_name = char(variable_name);
    ncwrite(filename, variable_name, data);
    if ~logical(verify_readback)
        return
    end

    variable_info = ncinfo(filename, variable_name);
    max_chunk_elements = 5e6;
    if prod(variable_info.Size) <= max_chunk_elements
        verify_values(ncread(filename, variable_name), data, ...
            filename, variable_name);
        return
    end

    slab_elements = prod(variable_info.Size(1:end-1));
    slab_count = max(1, floor(max_chunk_elements / slab_elements));
    for first = 1:slab_count:variable_info.Size(end)
        start = ones(1, numel(variable_info.Size));
        start(end) = first;
        count = variable_info.Size;
        count(end) = min(slab_count, variable_info.Size(end)-first+1);
        actual = ncread(filename, variable_name, start, count);

        indices = repmat({':'}, 1, numel(variable_info.Size));
        indices{end} = first:first+count(end)-1;
        expected = data(indices{:});
        verify_values(actual, expected, filename, variable_name);
    end
end

function verify_values(actual, expected, filename, variable_name)
    if isnumeric(actual) || islogical(actual)
        expected = cast(expected, 'like', actual);
    end
    if ~isequal(size(actual), size(expected))
        error('write_netcdf_variable:SizeMismatch', ...
            'Variable "%s" in "%s" has readback size %s; expected %s.', ...
            variable_name, filename, mat2str(size(actual)), ...
            mat2str(size(expected)));
    end

    if ~isequaln(actual, expected)
        error('write_netcdf_variable:ReadbackMismatch', ...
            'Readback verification failed for variable "%s" in "%s".', ...
            variable_name, filename);
    end
end
