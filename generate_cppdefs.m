function generate_cppdefs(enable_tnudg, enable_uvnudg, enable_convec, filename)
    if nargin < 1 || isempty(enable_tnudg)
        enable_tnudg = false;
    end
    if nargin < 2 || isempty(enable_uvnudg)
        enable_uvnudg = false;
    end
    if nargin < 3 || isempty(enable_convec)
        enable_convec = true;   % LMD_CONVEC on unless the experiment says otherwise
    end

    if ischar(enable_tnudg) || isstring(enable_tnudg)
        filename = enable_tnudg;
        enable_tnudg = false;
        enable_uvnudg = false;
        enable_convec = true;
    elseif ischar(enable_uvnudg) || isstring(enable_uvnudg)
        filename = enable_uvnudg;
        enable_uvnudg = false;
        enable_convec = true;
    elseif ischar(enable_convec) || isstring(enable_convec)
        filename = enable_convec;
        enable_convec = true;
    elseif nargin < 4 || isempty(filename)
        filename = 'cppdefs.opt';
    end

    [fid, cleaner] = open_text_file_for_write(filename);

    % Header
    fprintf(fid, '/* This is "cppdefs.opt": MODEL CONFIGURATION FILE\n');
    fprintf(fid, '   ==== == ============ ===== ============= ==== */\n\n');
    fprintf(fid, '/*\n');
    fprintf(fid, ' * CHOOSE ONLY ONE PRIMARY FLAG FOR SWITCH LIST BELOW\n');
    fprintf(fid, ' */\n\n');

    % Primary flag
    fprintf(fid, '#define NEPTUNE_IDEAL\n\n');

    % Start block
    fprintf(fid, '#if defined NEPTUNE_IDEAL\n\n');

    % Basics
    fprintf(fid, '        /* Basics */\n');
    fprintf(fid, '# define SOLVE3D\n');
    fprintf(fid, '# define UV_ADV\n');
    fprintf(fid, '# define UV_COR\n');
    fprintf(fid, '# define ADV_ISONEUTRAL\n');
    fprintf(fid, '# define DIAGNOSTICS\n\n');

    if enable_tnudg
        fprintf(fid, '# define TNUDG_Huaiyu\n\n');
    else
        fprintf(fid, '# undef TNUDG_Huaiyu\n\n');
    end

    if enable_uvnudg
        fprintf(fid, '# define UVNUDG_Huaiyu\n\n');
    else
        fprintf(fid, '# undef UVNUDG_Huaiyu\n\n');
    end
    fprintf(fid, '/* # undef SPONGE */\n\n');

    % Equation of State (commented out)
    fprintf(fid, '        /* Equation of State */\n');
    fprintf(fid, '/* # define NONLIN_EOS */\n');
    fprintf(fid, '/* # define SPLIT_EOS  */\n');
    fprintf(fid, '/* # define SALINITY   */\n\n');

    % Mixing
    fprintf(fid, '        /* Mixing */\n');
    fprintf(fid, '        /*        - lateral */\n');
    fprintf(fid, '# define UV_VIS2\n');
    fprintf(fid, '# define TS_DIF2\n\n');

    fprintf(fid, '        /*        - vertical */\n');
    fprintf(fid, '# define LMD_MIXING\n');
    fprintf(fid, '# define LMD_KPP\n');
    fprintf(fid, '# define LMD_NONLOCAL\n');
    fprintf(fid, '# define LMD_RIMIX\n');
    if enable_convec
        fprintf(fid, '# define LMD_CONVEC\n');
    else
        fprintf(fid, '# undef LMD_CONVEC\n');
    end
    fprintf(fid, '# define LMD_BKPP\n\n');

    % Grid configuration
    fprintf(fid, '        /* Grid Configuration */\n');
    fprintf(fid, '/* # define CURVGRID  */\n');
    fprintf(fid, '/* # define SPHERICAL */\n');
    fprintf(fid, '# define MASKING\n\n');

    % Boundaries
    fprintf(fid, '        /* Boundaries */\n');
    fprintf(fid, '# define NS_PERIODIC\n');
    fprintf(fid, '# define EW_PERIODIC\n\n');

    % Restart
    fprintf(fid, '        /* Restart */\n');
    fprintf(fid, '# define EXACT_RESTART\n\n');

    % Alternate case
    fprintf(fid, '#elif defined DUMMY_CASE\n\n');
    fprintf(fid, '# define AVERAG\n\n');
    fprintf(fid, '#endif\n\n');

    % Include
    fprintf(fid, '#include "set_global_definitions.h"\n');

    clear cleaner;
    fprintf('cppdefs.opt written to "%s"\n', filename);
end
