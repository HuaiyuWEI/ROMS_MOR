function generate_BasicDiag_file(...
    wrt_rst, per_rst, nrpf_rst, ...
    wrt_his, per_his, nrpf_his, ...
    wrt_avg, per_avg, nrpf_avg, ...
    filename)

    if nargin < 12
        filename = 'ocean_vars.opt';
    end

    fid = fopen(filename, 'w');
    fprintf(fid, '      ! ****************************************************************\n');
    fprintf(fid, '      ! user inputs\n');

    % Restart File Output
    fprintf(fid, '      logical,parameter :: wrt_file_rst      = .%s.\n', tfstr(wrt_rst));
    fprintf(fid, '      real,parameter    :: output_period_rst = %d\n', per_rst);
    fprintf(fid, '      integer,parameter :: nrpf_rst          = %d\n\n', nrpf_rst);

    % History File Output
    fprintf(fid, '      logical,parameter :: wrt_file_his      = .%s.\n', tfstr(wrt_his));
    fprintf(fid, '      real,parameter    :: output_period_his =  %d\n', per_his);
    fprintf(fid, '      integer,parameter :: nrpf_his          =  %d\n', nrpf_his);

    % History Output Variables
    fprintf(fid, '      logical,parameter :: wrt_Z =.true.,\n');
    fprintf(fid, '     &                     wrt_Ub=.true.,\n');
    fprintf(fid, '     &                     wrt_Vb=.true.,\n');
    fprintf(fid, '     &                     wrt_U=.true.,\n');
    fprintf(fid, '     &                     wrt_V=.true.,\n');
    fprintf(fid, '     &                     wrt_R=.false.,\n');
    fprintf(fid, '     &                     wrt_O=.true.,\n');
    fprintf(fid, '     &                     wrt_W=.true.,\n');
    fprintf(fid, '     &                     wrt_Akv=.true.,\n');
    fprintf(fid, '     &                     wrt_Akt=.true.,\n');
    fprintf(fid, '     &                     wrt_Aks=.false.,\n');
    fprintf(fid, '     &                     wrt_Hbls=.true.,\n');
    fprintf(fid, '     &                     wrt_Hbbl=.true.\n\n');

    % Averaged Output
    fprintf(fid, '      logical,parameter :: wrt_file_avg      = .%s.\n', tfstr(wrt_avg));
    fprintf(fid, '      real,parameter    :: output_period_avg = %d\n', per_avg);
    fprintf(fid, '      integer,parameter :: nrpf_avg          = %d\n', nrpf_avg);

    % Averaged Output Variables
    fprintf(fid, '      logical,parameter :: wrt_avg_Z =.true.,\n');
    fprintf(fid, '     &                     wrt_avg_Ub=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_Vb=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_U=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_V=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_R=.false.,\n');
    fprintf(fid, '     &                     wrt_avg_O=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_W=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_Akv=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_Akt=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_Aks=.false.,\n');
    fprintf(fid, '     &                     wrt_avg_Hbls=.true.,\n');
    fprintf(fid, '     &                     wrt_avg_Hbbl=.true.\n\n');

    % Final code check flag
    fprintf(fid, '      logical :: code_check = .false.\n');
    fprintf(fid, '      ! end user inputs\n');
    fprintf(fid, '      ! ****************************************************************\n');

    fclose(fid);
    fprintf('Output configuration file written to "%s"\n', filename);
end

function s = tfstr(tf)
    % Convert logical to Fortran .true./.false. string
    if tf
        s = 'true';
    else
        s = 'false';
    end
end
