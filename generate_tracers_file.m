function generate_tracers_file(passive_names, passive_long_names, passive_units, filename)
%GENERATE_TRACERS_FILE Write tracer metadata for UCLA-ROMS passive tracers.

    if nargin < 1
        error('passive_names is required.');
    end
    if nargin < 2 || isempty(passive_long_names)
        passive_long_names = passive_names;
    end
    if nargin < 3 || isempty(passive_units)
        passive_units = repmat({'nondimensional'}, size(passive_names));
    end
    if nargin < 4 || isempty(filename)
        filename = 'tracers.opt';
    end

    passive_names = cellstr(passive_names);
    passive_long_names = cellstr(passive_long_names);
    passive_units = cellstr(passive_units);

    if numel(passive_long_names) ~= numel(passive_names) || numel(passive_units) ~= numel(passive_names)
        error('passive_names, passive_long_names, and passive_units must have the same length.');
    end

    validate_fortran_strings(passive_names, 42, 'passive_names');
    validate_fortran_strings(passive_units, 42, 'passive_units');
    validate_fortran_strings(passive_long_names, 60, 'passive_long_names');

    [fid, cleaner] = open_text_file_for_write(filename);

    fprintf(fid, '! User notes:\n');
    fprintf(fid, '!   1) Passive tracer variable names must match variables in neptune_init.nc.\n');
    fprintf(fid, '!   2) t_ana_frc=2 uses source-code zero flux for no tracer source/sink.\n');
    fprintf(fid, '!   3) nt and nt_passive are set in param.opt.\n\n');

    fprintf(fid, '! USER INPUT 1:\n');
    fprintf(fid, '      interp_t = 0 ! Factor 2 scaling only for now\n\n');

    fprintf(fid, '! - TEMP & SALT:\n\n');
    fprintf(fid, '      wrt_t     (itemp)=.True.;    wrt_t_avg(itemp)=.True.\n');
    fprintf(fid, '      wrt_t_dia (itemp)=.True.\n\n');
    fprintf(fid, '#ifdef SALINITY\n');
    fprintf(fid, '      wrt_t     (isalt)=.True.;    wrt_t_avg(isalt)=.True.\n');
    fprintf(fid, '      wrt_t_dia (isalt)=.True.\n');
    fprintf(fid, '#endif\n\n');

    fprintf(fid, '! - PASSIVE TRACERS:\n\n');
    for itrc = 1:numel(passive_names)
        name = fortran_string(passive_names{itrc});
        units = fortran_string(passive_units{itrc});
        lname = fortran_string(passive_long_names{itrc});

        fprintf(fid, '      itot=itot+1\n');
        fprintf(fid, '      wrt_t     (itot)=.True.;    wrt_t_avg(itot)=.True.\n');
        fprintf(fid, '      wrt_t_dia (itot)=.True.\n');
        fprintf(fid, '      t_vname(itot)=''%s'';        t_units(itot)=''%s''\n', name, units);
        fprintf(fid, '      t_tname(itot)=''frc_time'';  t_ana_frc(itot)=2\n');
        fprintf(fid, '      t_lname(itot)=''%s''\n\n', lname);
    end

    fprintf(fid, '! - BGC TRACERS: (must come after passive tracers)\n');

    clear cleaner;
    fprintf('Tracer config file written to "%s"\n', filename);
end

function validate_fortran_strings(values, max_len, label)
    for i = 1:numel(values)
        if length(char(values{i})) > max_len
            error('%s{%d} exceeds the %d-character ROMS metadata limit.', label, i, max_len);
        end
    end
end

function value = fortran_string(value)
    value = char(value);
    value = strrep(value, '''', '''''');
end
