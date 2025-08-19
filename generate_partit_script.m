function generate_partit_script(core_num_x, core_num_y, filename)
    % generate_partit_script(core_num_x, core_num_y, [filename])
    % Writes ./Neptune_input/<filename> that runs:
    %   partit core_num_x core_num_y neptune_grid.nc neptune_init.nc neptune_frc.nc

    if nargin < 3 || isempty(filename)
        filename = 'do_partit.sh';
    end

    % Ensure output dir exists
    outdir = fullfile('.', 'Neptune_input');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Full path for the script
    filepath = fullfile(outdir, filename);

    % Write script
    fid = fopen(filepath, 'w');
    if fid == -1
        error('Cannot open "%s" for writing.', filepath);
    end
    fprintf(fid, '#!/bin/bash\n\n');
    fprintf(fid, 'partit %d %d neptune_grid.nc neptune_init.nc neptune_frc.nc\n', core_num_x, core_num_y);
    fclose(fid);

    % Make executable on Unix/macOS only (Windows doesn't support +x via fileattrib)
    if isunix || ismac
        fileattrib(filepath, '+x', 'a');
    end

    fprintf('Script "%s" created at "%s" with core_num_x=%d, core_num_y=%d\n', filename, outdir, core_num_x, core_num_y);
end
