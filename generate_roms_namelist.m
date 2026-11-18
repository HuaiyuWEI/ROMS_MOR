function generate_roms_namelist(...
    NTIMES, dt, THETA_S, THETA_B, hc, ...
    rho0, Zob, Tcoef, T0, Scoef, S0, ...
    NDTFAST, frc_files, filename)

    if nargin < 11
        error('Eleven inputs are required before the optional NDTFAST, frc_files, filename.');
    end
    if nargin < 12 || isempty(NDTFAST)
        NDTFAST = 20;   % barotropic steps per baroclinic step
    end
    validateattributes(NDTFAST, {'numeric'}, ...
        {'scalar', 'integer', 'positive'}, mfilename, 'NDTFAST');
    if nargin < 13 || isempty(frc_files)
        frc_files = {'Neptune_input/neptune_frc.nc'};
    elseif ischar(frc_files) || isstring(frc_files)
        frc_files = {char(frc_files)};
    end
    if nargin < 14 || isempty(filename)
        filename = 'neptune.in';
    end

    [fid, cleaner] = open_text_file_for_write(filename);

    % Title
    fprintf(fid, 'title:\n');
    fprintf(fid, '   Idealized simulation of eddy-driven topostrophy, using input files to specify grid, initial, and forcing\n\n');

    % Time stepping
    fprintf(fid, 'time_stepping: NTIMES   dt[sec]  NDTFAST  NINFO\n');
    fprintf(fid, '               %d      %d       %d    144\n\n', NTIMES, dt, NDTFAST);

    % S-coordinates
    fprintf(fid, 'S-coord: THETA_S,   THETA_B,    hc (m)\n');
    fprintf(fid, '          %.1fD0        %.1fD0    %d\n\n', THETA_S, THETA_B, hc);

    % Grid
    fprintf(fid, 'grid:  filename\n');
    fprintf(fid, '     Neptune_input/neptune_grid.nc\n\n');

    % Forcing (one filename per line; the model searches them in order)
    fprintf(fid, 'forcing: filename\n');
    for k = 1:numel(frc_files)
        fprintf(fid, '     %s\n', frc_files{k});
    end
    fprintf(fid, '\n');

    % Initial
    fprintf(fid, 'initial: NRREC  filename\n');
    fprintf(fid, '          1\n');
    fprintf(fid, '     Neptune_input/neptune_init.nc\n\n');

    % Output name
    fprintf(fid, 'output_root_name:\n');
    fprintf(fid, '      neptune\n\n');

    % Lateral viscosity
    fprintf(fid, 'lateral_visc:   VISC2,    VISC4    [m^2/sec for all]\n');
    fprintf(fid, '                 0.       0.\n\n');

    % Reference density
    fprintf(fid, 'rho0:\n');
    fprintf(fid, '      %.6g\n\n', rho0);

    % Sponge viscosity: read_inp requires the v_sponge keyword whether or
    % not SPONGE is defined; 0 keeps the sponge inactive.
    fprintf(fid, 'v_sponge:  V_SPONGE [m^2/sec]\n');
    fprintf(fid, '           0.\n\n');

    % Tracer diffusivity
    fprintf(fid, 'tracer_diff2: TNU2(1:NT)           [m^2/sec for all]\n');
    fprintf(fid, ' 0. 0. 0. 0. 0. 0. 0. 0. 0.\n\n');

    % Bottom drag
    fprintf(fid, 'bottom_drag:     RDRG [m/s],  RDRG2,  Zob [m],  Cdb_min, Cdb_max\n');
    fprintf(fid, '                  0.E-4       0.0E-3   %.1E     1.E-4    1.E-1\n\n', Zob);

    % EOS
    fprintf(fid, 'lin_rho_eos:  Tcoef    T0    Scoef   S0\n');
    fprintf(fid, '              %.2f   %.1f   %.1f  %.1f\n', Tcoef, T0, Scoef, S0);

    clear cleaner;

    fprintf('ROMS namelist file written to "%s"\n', filename);
end
