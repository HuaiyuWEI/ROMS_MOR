function assert_netcdf_no_nan(ncfile)
%ASSERT_NETCDF_NO_NAN Fail when any numeric NetCDF variable contains NaNs.

    ncfile = char(ncfile);
    if exist(ncfile, 'file') ~= 2
        error('assert_netcdf_no_nan:FileNotFound', ...
            'NetCDF file not found: "%s".', ncfile);
    end

    file_info = ncinfo(ncfile);
    variables_with_nan = cell(1, 0);
    max_chunk_elements = 5e6;

    for ivar = 1:numel(file_info.Variables)
        variable = file_info.Variables(ivar);
        if strcmp(variable.Datatype, 'char')
            continue
        end

        if prod(variable.Size) == 1
            contains_nan = any(isnan(ncread(ncfile, variable.Name)), 'all');
        else
            slab_elements = prod(variable.Size(1:end-1));
            slab_count = max(1, floor(max_chunk_elements / slab_elements));
            contains_nan = false;
            for first = 1:slab_count:variable.Size(end)
                start = ones(1, numel(variable.Size));
                start(end) = first;
                count = variable.Size;
                count(end) = min(slab_count, variable.Size(end)-first+1);
                data = ncread(ncfile, variable.Name, start, count);
                if any(isnan(data(:)))
                    contains_nan = true;
                    break
                end
            end
        end

        if contains_nan
            variables_with_nan{end+1} = variable.Name; %#ok<AGROW>
        end
    end

    if ~isempty(variables_with_nan)
        error('assert_netcdf_no_nan:NaNDetected', ...
            'NaN values found in "%s": %s.', ncfile, ...
            strjoin(variables_with_nan, ', '));
    end

    fprintf('No NetCDF variables in "%s" contain NaN values.\n', ncfile);
end
