function generate_energy_budget_file( ...
    do_energy_budget, output_avg, output_period, nrpf, filename)
%GENERATE_ENERGY_BUDGET_FILE Write energy_budget.opt.
% Kinetic-energy budget output module:
%   do_energy_budget : master switch
%   output_avg       : averaged (=T) or snapshot (=F) output
%   output_period    : output period in seconds
%   nrpf             : records per file
% NOTE: when enabled, the module requires diag_uv = .true. and
% diag_norm = .true. in diagnostics.opt (the model stops otherwise).

    if nargin < 4
        error('Four inputs are required before the optional filename.');
    end
    if nargin < 5 || isempty(filename)
        filename = 'energy_budget.opt';
    end

    do_energy_budget = logical(do_energy_budget);
    output_avg = logical(output_avg);
    validateattributes(output_period, {'numeric'}, {'scalar', 'positive'}, mfilename, 'output_period');
    validateattributes(nrpf, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'nrpf');

    tf = {'.false.','.true.'};
    tfstr = @(x) tf{1 + (x~=0)};

    [fid, cleaner] = open_text_file_for_write(filename);

    fprintf(fid, '\n');
    fprintf(fid, '      !-----------------------------------------------------------\n');
    fprintf(fid, '      ! User inputs for energy_budget output module\n');
    fprintf(fid, '      !\n');
    fprintf(fid, '      logical,public,parameter :: do_energy_budget = %s\n', tfstr(do_energy_budget));
    fprintf(fid, '      logical,public,parameter :: output_avg = %s\n', tfstr(output_avg));
    fprintf(fid, '      real   ,parameter :: output_period = %g    ! in seconds\n', output_period);
    fprintf(fid, '      integer,parameter :: nrpf   = %d             ! number of frames per file\n', nrpf);
    fprintf(fid, '      !\n');
    fprintf(fid, '      !-----------------------------------------------------------\n');

    clear cleaner;
    fprintf('Energy budget config file written to "%s"\n', filename);
end
