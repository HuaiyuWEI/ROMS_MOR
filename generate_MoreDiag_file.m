function generate_MoreDiag_file( ...
    diag_avg, output_period, nrpf, ...
    diag_uv, diag_trc, diag_pflx, ...
    filename)

    if nargin < 7 || isempty(filename)
        filename = 'diagnostics.opt';
    end

    % Ensure logicals
    diag_avg = logical(diag_avg);
    diag_uv  = logical(diag_uv);
    diag_trc = logical(diag_trc);
    diag_pflx= logical(diag_pflx);

    % Map logical -> Fortran logical string
    tf = {'.false.','.true.'};
    tfstr = @(x) tf{1 + (x~=0)};   % returns '.false.' or '.true.'

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open "%s" for writing.', filename);
    end

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
    fprintf(fid, '      !        select which tracers you want in wrt_t_diag array in tracers.opt\n');
    fprintf(fid, '      !     4) for an example try Examples/Diagnostics/\n');
    fprintf(fid, '      !     5) averaging frequency = rec_rate * time_step\n');
    fprintf(fid, '      !        averaging is expensive and will likely slow your simulation by 30-40%%\n');
    fprintf(fid, '      !        with u, v, temp and salinity diagnostics.\n');
    fprintf(fid, '      !     6) history (snap-shot) diagnostics is relatively inexpensive.\n\n');

    fprintf(fid, '      logical, parameter         :: diag_avg      = %s    ! compute history (=F) or averages (=T)\n', tfstr(diag_avg));
    fprintf(fid, '      integer, parameter         :: output_period =  %d       ! output period\n', output_period);
    fprintf(fid, '      integer, parameter         :: nrpf          =  %d        ! total recs per file\n\n', nrpf);

    fprintf(fid, '      logical, parameter, public :: diag_uv       = %s   ! Momentum diagnostics\n', tfstr(diag_uv));
    fprintf(fid, '      logical, parameter, public :: diag_trc      = %s   ! Selected tracers diagnostics\n', tfstr(diag_trc));
    fprintf(fid, '      logical, parameter, public :: diag_pflx     = %s   ! Baroclinic pressure fluxes\n\n', tfstr(diag_pflx));

    fprintf(fid, '      real,    parameter         :: timescale     = %d ! timescale for filtering (in seconds)\n\n', 24*3600);

    % If you truly want an identifier, not a string, keep %s with bare text:
    fprintf(fid, '      integer, parameter         :: diag_prec     = %s ! Precision of output variables (nf90_float/nf90_double)\n\n', 'nf90_double');

    fprintf(fid, '      ! End of user inputs\n');
    fprintf(fid, '      ! ***************************************************************\n');

    fclose(fid);
    fprintf('Diagnostics config file written to "%s"\n', filename);
end
