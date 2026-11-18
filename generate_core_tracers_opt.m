function generate_core_tracers_opt(wrt_salt, filename)
%GENERATE_CORE_TRACERS_OPT Write tracers.opt.
% temp and salt are both mandatory core tracers (nt >= 2). This file
% controls their output flags: with wrt_salt = false the salt field
% (all zeros in temperature-only setups) is computed but kept out of
% the his/avg/dia files.

    if nargin < 1 || isempty(wrt_salt)
        wrt_salt = false;
    end
    if nargin < 2 || isempty(filename)
        filename = 'tracers.opt';
    end

    tf = {'.false.','.true.'};
    tfstr = @(x) tf{1 + (x~=0)};

    [fid, cleaner] = open_text_file_for_write(filename);

    fprintf(fid, '! tracers.opt\n');
    fprintf(fid, '! temp and salt are mandatory core tracers; this file controls\n');
    fprintf(fid, '! their output (and would declare additional passive tracers).\n\n');
    fprintf(fid, '      interp_t = 0 ! interpolate forcing from a coarser input grid (=1) or not (=0)\n\n');
    fprintf(fid, '      wrt_t(itemp)     = .true.\n');
    fprintf(fid, '      wrt_t_avg(itemp) = .true.\n');
    fprintf(fid, '      wrt_t_dia(itemp) = .true.\n\n');
    fprintf(fid, '      wrt_t(isalt)     = %s\n', tfstr(wrt_salt));
    fprintf(fid, '      wrt_t_avg(isalt) = %s\n', tfstr(wrt_salt));
    fprintf(fid, '      wrt_t_dia(isalt) = %s\n', tfstr(wrt_salt));

    clear cleaner;
    fprintf('Core tracers config written to "%s"\n', filename);
end
