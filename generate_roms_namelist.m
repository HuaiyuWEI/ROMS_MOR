function generate_roms_namelist(...
    NTIMES, dt, THETA_S, THETA_B, hc, ...
    rho0,Zob, Tcoef, T0, Scoef, S0, ...
    filename)

    if nargin < 11
        error('All 11 arguments are required.');
    end
    if nargin < 12
        filename = 'neptune.in';
    end

    fid = fopen(filename, 'w');

    % Title
    fprintf(fid, 'title:\n');
    fprintf(fid, '   Idealized simulation of eddy-driven topostrophy, using input files to specify grid, initial, and forcing\n\n');

    % Time stepping
    fprintf(fid, 'time_stepping: NTIMES   dt[sec]  NDTFAST  NINFO\n');
    fprintf(fid, '               %d      %d       20    144\n\n', NTIMES, dt);

    % S-coordinates
    fprintf(fid, 'S-coord: THETA_S,   THETA_B,    hc (m)\n');
    fprintf(fid, '          %.1fD0        %.1fD0    %d\n\n', THETA_S, THETA_B, hc);

    % Grid
    fprintf(fid, 'grid:  filename\n');
    fprintf(fid, '     Neptune_input/neptune_grid.nc\n\n');

    % Forcing
    fprintf(fid, 'forcing: filename\n');
    fprintf(fid, '     Neptune_input/neptune_frc.nc\n\n');

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

    % Reference density (customizable)
    fprintf(fid, 'rho0:\n');
    fprintf(fid, '      %.6g\n\n', rho0);

    % Tracer diffusivity
    fprintf(fid, 'tracer_diff2: TNU2(1:NT)           [m^2/sec for all]\n');
    fprintf(fid, ' 0. 0. 0. 0. 0. 0. 0. 0. 0.\n\n');

    % Bottom drag
    fprintf(fid, 'bottom_drag:     RDRG [m/s],  RDRG2,  Zob [m],  Cdb_min, Cdb_max\n');
    fprintf(fid, '                  0.E-4       0.0E-3   %.1E     1.E-4    1.E-1\n\n', Zob);

    % gamma2
    fprintf(fid, 'gamma2:\n');
    fprintf(fid, '                  1.D0\n\n');

    % EOS
    fprintf(fid, 'lin_rho_eos:  Tcoef    T0    Scoef   S0\n');
    fprintf(fid, '              %.2f   %.1f   %.1f  %.1f\n', Tcoef, T0, Scoef, S0);

    fclose(fid);

    fprintf('ROMS namelist file written to "%s"\n', filename);
end
