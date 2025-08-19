function generate_param_file(LLm, MMm, N, NP_XI, NP_ETA, NSUB_X, NSUB_E, nt, nt_passive, filename)
    if nargin < 10
        filename = 'param.opt';  % default filename
    end

    fid = fopen(filename, 'w');

    fprintf(fid, '! Parameter options file\n\n');

    fprintf(fid, '! Dimensions of Physical Grid and array dimensions:\n');
    fprintf(fid, '!----------- -- -------- ---- --- ----- -----------\n');
    fprintf(fid, '! LLm   Number of the internal points of the PHYSICAL grid in XI-\n');
    fprintf(fid, '! MMm   and ETA-directions, excluding physical side boundary points,\n');
    fprintf(fid, '!       peroodic ghost points, and MPI-margins (if any).\n\n');

    fprintf(fid, '! Domain subdivision parameters:\n');
    fprintf(fid, '!------- ----------- -----------\n');
    fprintf(fid, '! NP_XI,  NP_ETA     number of MPI subdomains in XI-, ETA-directions;\n');
    fprintf(fid, '! NSUB_X, NSUB_E     number of shared memory subdomains (tiles) in XI- and ETA-directions;\n\n');

    fprintf(fid, '! Number of tracers\n');
    fprintf(fid, '!------- -----------\n');
    fprintf(fid, '! nt  must be 2 or more if Salinity is defined\n\n');

    fprintf(fid, '      integer, parameter :: LLm=%d, MMm=%d, N=%d\n\n', LLm, MMm, N);
    fprintf(fid, '      integer, parameter :: NP_XI  = %d, NP_ETA = %d\n', NP_XI, NP_ETA);
    fprintf(fid, '      integer, parameter :: NSUB_X = %d, NSUB_E = %d\n\n', NSUB_X, NSUB_E);
    fprintf(fid, '      integer, parameter :: nt = %d\n', nt);
    fprintf(fid, '      integer, parameter :: nt_passive = %d\n\n', nt_passive);


    fclose(fid);
    fprintf('Parameter file written to "%s"\n', filename);
end
