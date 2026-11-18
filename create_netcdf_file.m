function [ncid, cleaner] = create_netcdf_file(filename, creation_mode)
%CREATE_NETCDF_FILE Create a NetCDF file with automatic cleanup on errors.
%   [NCID, CLEANER] = CREATE_NETCDF_FILE(FILENAME, CREATION_MODE) removes an
%   existing file, creates its parent directory when needed, and returns an
%   onCleanup object that closes NCID. Clear CLEANER after NETCDF.ENDDEF to
%   close the file before using high-level functions such as NCWRITE.

    if nargin < 2 || isempty(creation_mode)
        creation_mode = 'CLOBBER';
    end

    filename = char(filename);
    creation_mode = char(creation_mode);
    [parent_dir, ~, ~] = fileparts(filename);
    if ~isempty(parent_dir) && exist(parent_dir, 'dir') ~= 7
        [ok, message] = mkdir(parent_dir);
        if ~ok
            error('create_netcdf_file:CreateDirectoryFailed', ...
                'Could not create NetCDF parent directory "%s": %s', ...
                parent_dir, message);
        end
    end

    if exist(filename, 'file') == 2
        delete(filename);
    end

    ncid = netcdf.create(filename, creation_mode);
    cleaner = onCleanup(@() netcdf.close(ncid));
end
