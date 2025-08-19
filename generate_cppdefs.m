function generate_cppdefs(filename)
    if nargin < 1 || isempty(filename)
        filename = 'cppdefs.opt';
    end

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open "%s" for writing.', filename);
    end

    % ===== Header you requested =====
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
    fprintf(fid, '# define LMD_CONVEC\n');
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

    % Output & Restart
    fprintf(fid, '        /* Output Options */\n');
    fprintf(fid, '# define MASK_LAND_DATA\n\n');
    fprintf(fid, '        /* Restart */\n');
    fprintf(fid, '# define EXACT_RESTART\n\n');

    % Alternate case
    fprintf(fid, '#elif defined DUMMY_CASE\n\n');
    fprintf(fid, '# define AVERAG\n\n');
    fprintf(fid, '#endif\n\n');

    % Include
    fprintf(fid, '#include "set_global_definitions.h"\n');

    fclose(fid);
    fprintf('cppdefs.opt written to "%s"\n', filename);
end
