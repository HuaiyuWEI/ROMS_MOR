function generate_MoreDiag_file( ...
    diag_avg, output_period, nrpf, ...
    diag_uv, diag_trc, diag_pflx, ...
    diag_norm, filename)

    if nargin < 6
        error('Six inputs are required before the optional diag_norm, filename.');
    end
    if nargin < 8 || isempty(filename)
        filename = 'diagnostics.opt';
    end
    if nargin >= 7 && (ischar(diag_norm) || isstring(diag_norm))
        filename = char(diag_norm);   % the file name may also be given as 7th argument
        diag_norm = [];
    end
    if nargin < 7 || isempty(diag_norm)
        % false = dz-integrated momentum terms (m^2/s^2); true = the same
        % terms divided by dz (m/s^2). The energy_budget module requires
        % diag_norm = true.
        diag_norm = false;
    end

    % Ensure logicals
    diag_avg = logical(diag_avg);
    diag_uv  = logical(diag_uv);
    diag_trc = logical(diag_trc);
    diag_pflx= logical(diag_pflx);
    diag_norm= logical(diag_norm);

    % Map logical -> Fortran logical string
    tf = {'.false.','.true.'};
    tfstr = @(x) tf{1 + (x~=0)};   % returns '.false.' or '.true.'

    validateattributes(output_period, {'numeric'}, {'scalar', 'positive'}, mfilename, 'output_period');
    validateattributes(nrpf, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'nrpf');

    [fid, cleaner] = open_text_file_for_write(filename);

    fprintf(fid, '      ! ***************************************************************\n');
    fprintf(fid, '      ! User inputs:\n');
    fprintf(fid, '      !\n');
    fprintf(fid, '      !   Momentum terms:\n');
    fprintf(fid, '      !     Pressure gradient; coriolis; adv.; dissipation from adv.; h. mix.; and v. mix.\n');
    fprintf(fid, '      !\n');
    fprintf(fid, '      !   Tracer terms:\n');
    fprintf(fid, '      !     H. adv.; dissipation from h. adv.; v. adv.; dissipation from v. adv.; h. mix.; and v. mix.\n');
    fprintf(fid, '      !\n');
    fprintf(fid, '      !   Notes:\n');
    fprintf(fid, '      !     1) need DIAGNOSTICS flag in cppdefs.opt\n');
    fprintf(fid, '      !     2) units are du/dt*dz - m^2/s^2 (or v/w/tracer) - i.e. vertically integrated in cell\n');
    fprintf(fid, '      !     3) if tracer diagnostics required set diag_trc==T,\n');
    fprintf(fid, '      !        select which tracers you want in wrt_t_dia array in tracers.opt\n');
    fprintf(fid, '      !     4) for an example try Examples/Diagnostics/\n');
    fprintf(fid, '      !     5) averaging frequency = rec_rate * time_step\n');
    fprintf(fid, '      !        averaging is expensive and will likely slow your simulation by 30-40%%\n');
    fprintf(fid, '      !        with u, v, temp and salinity diagnostics.\n');
    fprintf(fid, '      !     6) history (snap-shot) diagnostics is relatively inexpensive.\n\n');

    fprintf(fid, '      ! normalization should typically be on for 3d budgets and off for vertically integrated versions\n');
    fprintf(fid, '      logical, parameter, public :: diag_norm     = %s   ! Normalization of diagnostic terms\n\n', tfstr(diag_norm));
    fprintf(fid, '      logical, parameter         :: diag_avg      = %s    ! compute history (=F) or averages (=T)\n', tfstr(diag_avg));
    fprintf(fid, '      integer, parameter         :: output_period =  %d       ! output period\n', output_period);
    fprintf(fid, '      integer, parameter         :: nrpf          =  %d        ! total recs per file\n\n', nrpf);

    fprintf(fid, '      logical, parameter, public :: diag_uv       = %s   ! Momentum diagnostics\n', tfstr(diag_uv));
    fprintf(fid, '      logical, parameter, public :: diag_trc      = %s   ! Selected tracers diagnostics\n', tfstr(diag_trc));
    fprintf(fid, '      logical, parameter, public :: diag_pflx     = %s   ! Baroclinic pressure fluxes\n\n', tfstr(diag_pflx));

    fprintf(fid, '      real,    parameter         :: timescale     = %d ! timescale for filtering (in seconds)\n\n', 24*3600);

    % nf90_double is written as a bare identifier, not as a quoted string
    fprintf(fid, '      integer, parameter         :: diag_prec     = %s ! Precision of output variables (nf90_float/nf90_double)\n\n', 'nf90_double');

    fprintf(fid, '      !  If diag_trc=.true., select for which tracers diagnostics are desired:\n');
    fprintf(fid, '      !  Set wrt_t_dia(...) in tracers.opt.\n\n');

    fprintf(fid, '      ! End of user inputs\n');
    fprintf(fid, '      ! ***************************************************************\n');

    clear cleaner;
    fprintf('Diagnostics config file written to "%s"\n', filename);
end
