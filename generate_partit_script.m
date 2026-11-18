function generate_partit_script(core_num_x, core_num_y, filename)
    % generate_partit_script(core_num_x, core_num_y, [filename])
    % Writes ./Neptune_input/<filename>. If neptune_init.nc is present, it is
    % partitioned with grid/forcing; restart-only directories can partition
    % just grid/forcing.

    if nargin < 3 || isempty(filename)
        filename = 'do_partit.sh';
    end
    validateattributes(core_num_x, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'core_num_x');
    validateattributes(core_num_y, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'core_num_y');

    % Ensure output dir exists
    outdir = fullfile('.', 'Neptune_input');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Full path for the script
    filepath = fullfile(outdir, filename);

    % Write script
    [fid, cleaner] = open_text_file_for_write(filepath);
    fprintf(fid, '#!/bin/bash\n\n');
    fprintf(fid, 'set -euo pipefail\n\n');
    fprintf(fid, 'files="neptune_grid.nc"\n');
    fprintf(fid, 'if [[ -f neptune_init.nc ]]; then files="$files neptune_init.nc"; fi\n');
    fprintf(fid, 'files="$files neptune_frc.nc"\n');
    fprintf(fid, '# separate zero-flux forcing file (ridge experiments)\n');
    fprintf(fid, 'if [[ -f neptune_flux.nc ]]; then files="$files neptune_flux.nc"; fi\n\n');
    fprintf(fid, '# Remove partitions from an older processor layout before repartitioning\n');
    fprintf(fid, '# (any rank-suffix width: .00 for 20 ranks, .000 for 128, ...).\n');
    fprintf(fid, 'rm -f neptune_grid.[0-9]*.nc neptune_init.[0-9]*.nc \\\n');
    fprintf(fid, '      neptune_frc.[0-9]*.nc neptune_flux.[0-9]*.nc\n\n');
    fprintf(fid, 'partit %d %d $files\n', core_num_x, core_num_y);
    clear cleaner;

    if isunix || ismac
        fileattrib(filepath, '+x', 'a');
    end

    fprintf('Script "%s" created at "%s" with core_num_x=%d, core_num_y=%d\n', filename, outdir, core_num_x, core_num_y);
end
