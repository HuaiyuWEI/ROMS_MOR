function expDir = NewExp_Neptune_ridge_2D(parentName, opts)
% NEWEXP_NEPTUNE_RIDGE_2D  Generate a 2-D (along-ridge-uniform) ridge run
% driven by its parent 3-D run's time- and along-ridge-mean vertical mixing.
%
%   NewExp_Neptune_ridge_2D()                   default parent, default options
%   NewExp_Neptune_ridge_2D(parentName)
%   NewExp_Neptune_ridge_2D(parentName, opts)   opts fields listed below
%
% Everything except the along-ridge extent, the forcing and the vertical
% mixing is taken from the PARENT's own files rather than regenerated from
% the generator defaults, so edits made to the parent after it was
% generated (dt, restart period, ...) carry over, and the grid, vertical
% coordinate, initial stratification, EOS, bottom drag and nudging strips
% match the parent exactly:
%   grid / init / flux files -> the parent's, cut to opts.Ny eta rows
%   neptune_frc.nc           -> zero wind stress (no external forcing)
%   vertical mixing          -> KPP off; Akv/Akt prescribed (ANA_VMIX +
%                               VMIX_PRESCRIBED_Huaiyu) as the time- and
%                               y-mean of the parent's avg-file Akv/Akt
%   neptune.in, *.opt, *.F   -> the parent's, with the minimal 2-D edits
%
% The ridge varies in x (xi) and is uniform along y (eta), so the 2-D run
% keeps all LLm cross-ridge points and only Ny interior points in eta,
% doubly periodic like the parent.
%
% parentName is the name of the parent's folder in the output root. The
% folder must hold the parent's joined averages files
% neptune_avg.<stamp>.nc, from which Akv and Akt are read.
% NewExp_Neptune_ridge.m writes ocean_vars.opt with the averages file
% switched off (wrt_file_avg), so switch it on before running the parent.
%
% opts (all optional):
%   expname     default: parent name with 'Ridge_' -> 'Ridge2D_', then
%               '_NoSurfMix' when surfaceMixingDepth > 0, then '_Ny<Ny>'.
%               The SLURM job name is the last '_' segment of the
%               experiment name (e.g. "Ny4").
%   Ny          interior along-ridge points (default 4). The minimum is 2:
%               with NP_ETA = 1 the tile is its own periodic neighbour and
%               roms_mpi packs rows 1:bf and ny-bf+1:ny (bf = 2), so ny < 2
%               would send a halo row as data. A y-uniform state stays
%               exactly uniform, so Ny only sets cost.
%   avgStamps   parent avg stamps to average (default: the last nAvg)
%   nAvg        number of trailing joined avg files to average (default 5)
%   foldMixing  symmetrize Akv/Akt across the crest (default false)
%   npXi        MPI tasks along xi; NP_ETA = 1 (default 20)
%   targetYear  run_roms_auto TARGET_YEAR (default 2030)
%   walltime    SLURM wall time per segment, hours (default 12)
%   overwrite   replace an existing 2-D experiment folder (default false)
%   productionDiags  switch on the online momentum and temperature budgets
%               (dz-normalized), the averages file, daily snapshots and the
%               KE budget (default true). The values are set explicitly,
%               whatever the parent's .opt files hold.
%   surfaceMixingDepth  depth [m] above which the parent's surface-enhanced
%               mixing is replaced by the interior background bgAkv/bgAkt
%               (default 1000; 0 keeps the full parent profile). The focus is
%               the mixing-driven abyssal flow. Checked per parent: error if
%               the cutoff comes within 500 m of the shallowest seabed;
%               warning if surface mixing reaches below it.
%   bgAkv, bgAkt  interior background viscosity / diffusivity [m^2/s]
%               (default 1e-4 / 1e-5 = nuwm / nuws in lmd_vmix.F, the
%               internal-wave floor KPP applies away from shear instability)
%
% Caveat: the time-mean Akt overstates the parent's actual mixing flux,
% because Akt and N^2 are anti-correlated.
%
% Huaiyu Wei

if nargin < 1 || isempty(parentName)
    parentName = 'Ridge_KPPconv_hr2000_tau075_alphaT2em4_zob1em3_ref';
end
if nargin < 2
    opts = struct();
end
opts = fillDefaults(opts, struct('expname', '', 'Ny', 4, 'avgStamps', {{}}, ...
    'nAvg', 5, 'foldMixing', false, 'npXi', 20, 'targetYear', 2030, ...
    'walltime', 12, 'overwrite', false, 'productionDiags', true, ...
    'surfaceMixingDepth', 1000, 'bgAkv', 1e-4, 'bgAkt', 1e-5));
validateattributes(opts.surfaceMixingDepth, {'numeric'}, {'scalar', 'nonnegative'}, mfilename, 'opts.surfaceMixingDepth');
validateattributes(opts.Ny, {'numeric'}, {'scalar', 'integer', '>=', 2}, mfilename, 'opts.Ny');
validateattributes(opts.npXi, {'numeric'}, {'scalar', 'integer', 'positive', '<', 128}, mfilename, 'opts.npXi');

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
config = newexp_defaults(fullfile(script_dir, 'NewExp_Neptune_ridge.m'));

parentDir = fullfile(config.output_root, parentName);
parentIn  = fullfile(parentDir, 'Neptune_input');
mustExist(parentDir, 'dir');
for f = {'neptune.in', 'cppdefs.opt', 'param.opt', 'Makefile'}
    mustExist(fullfile(parentDir, f{1}), 'file');
end
for f = {'neptune_grid.nc', 'neptune_init.nc', 'neptune_frc.nc', 'neptune_flux.nc'}
    mustExist(fullfile(parentIn, f{1}), 'file');
end

if isempty(opts.expname)
    opts.expname = regexprep(parentName, '^Ridge_', 'Ridge2D_');
    if strcmp(opts.expname, parentName)
        opts.expname = ['Ridge2D_' parentName];
    end
    if opts.surfaceMixingDepth > 0
        opts.expname = [opts.expname '_NoSurfMix'];
    end
    opts.expname = sprintf('%s_Ny%d', opts.expname, opts.Ny);
end
expDir = fullfile(config.output_root, opts.expname);
inDir  = fullfile(expDir, 'Neptune_input');

fprintf('\n2-D ridge experiment\n  parent : %s\n  new    : %s\n', parentDir, expDir);
if exist(expDir, 'dir') && ~opts.overwrite   % fail before any heavy reading
    error('NewExp_Neptune_ridge_2D:exists', ...
        '%s already exists; pass opts.overwrite = true to replace it.', expDir);
end


%% Parent dimensions

paramTxt = fileread(fullfile(parentDir, 'param.opt'));
LLm = parseInt(paramTxt, 'LLm');
MMm = parseInt(paramTxt, 'MMm');
N   = parseInt(paramTxt, 'N');
Ny  = opts.Ny;
fprintf('  parent grid LLm=%d MMm=%d N=%d  ->  2-D grid LLm=%d MMm=%d N=%d\n', ...
    LLm, MMm, N, LLm, Ny, N);


%% The parent must be along-ridge uniform, or cutting rows is meaningless

parentGrid = fullfile(parentIn, 'neptune_grid.nc');
checkVars = {parentGrid, 'h'; parentGrid, 'f'; parentGrid, 'pm'; parentGrid, 'pn'; ...
    parentGrid, 'mask_rho'};
gridInfo = ncinfo(parentGrid);
for nm = {'t_targ', 'nudgc', 'uv_nudgc'}
    if any(strcmp({gridInfo.Variables.Name}, nm{1}))
        checkVars(end+1, :) = {parentGrid, nm{1}}; %#ok<AGROW>
    end
end
checkVars(end+1, :) = {fullfile(parentIn, 'neptune_init.nc'), 'temp'};
for k = 1:size(checkVars, 1)
    A = ncread(checkVars{k, 1}, checkVars{k, 2});
    dev = max(abs(A - A(:, 1, :)), [], 'all');
    scale = max(abs(A), [], 'all');
    if dev > 1e-10 * max(scale, eps)
        error('NewExp_Neptune_ridge_2D:notUniform', ...
            ['Parent field %s varies along eta (max dev %.3g of %.3g); a 2-D cut ', ...
             'of a non-uniform parent would not represent it.'], checkVars{k, 2}, dev, scale);
    end
end
fprintf('  parent along-ridge uniformity verified (%d fields)\n', size(checkVars, 1));


%% Time- and along-ridge-mean mixing from the parent's avg files

avgAll = dir(fullfile(parentDir, 'neptune_avg.*.nc'));
avgAll = sort({avgAll(~cellfun(@isempty, ...
    regexp({avgAll.name}, '^neptune_avg\.\d{14}\.nc$', 'once'))).name});
if isempty(opts.avgStamps)
    if numel(avgAll) < opts.nAvg
        error('Only %d joined avg files in %s (need %d).', numel(avgAll), parentDir, opts.nAvg);
    end
    avgFiles = avgAll(end-opts.nAvg+1:end);
else
    avgFiles = cellfun(@(s) sprintf('neptune_avg.%s.nc', s), opts.avgStamps, 'UniformOutput', false);
    for f = avgFiles
        mustExist(fullfile(parentDir, f{1}), 'file');
    end
end
fprintf('  mixing mean over %d parent avg files: %s ... %s\n', ...
    numel(avgFiles), avgFiles{1}, avgFiles{end});

[AkvMean, AktMean] = deal(0);
for i = 1:numel(avgFiles)
    fn = fullfile(parentDir, avgFiles{i});
    AkvMean = AkvMean + yMean(ncread(fn, 'Akv')) / numel(avgFiles);
    AktMean = AktMean + yMean(ncread(fn, 'Akt')) / numel(avgFiles);
end
% The avg files leave Akt unset (zero) in the periodic ghost columns
% (xi_rho index 1 and LLm+2). Fill both fields from their periodic
% counterparts so the input does not rely on the model's halo exchange.
if ~isempty(regexp(fileread(fullfile(parentDir, 'cppdefs.opt')), ...
        '(?m)^\s*#\s*define\s+EW_PERIODIC', 'once'))
    AkvMean([1 end], :, :) = AkvMean([end-1 2], :, :);
    AktMean([1 end], :, :) = AktMean([end-1 2], :, :);
end
if opts.foldMixing   % scalar diffusivities are even about the crest
    AkvMean = 0.5 * (AkvMean + flip(AkvMean, 1));
    AktMean = 0.5 * (AktMean + flip(AktMean, 1));
end

% Surface-enhanced mixing -> interior background (abyssal focus)
if opts.surfaceMixingDepth > 0
    [AkvMean, AktMean] = removeSurfaceMixing(AkvMean, AktMean, parentGrid, opts);
end

if any(~isfinite([AkvMean(:); AktMean(:)])) || any([AkvMean(:); AktMean(:)] < 0)
    error('Parent mean Akv/Akt contains NaN/Inf or negative values.');
end
% Interior faces only: ROMS takes the seabed and surface fluxes from the
% drag/flux boundary conditions, and KPP leaves the surface face at zero.
inner = @(A) A(:, :, 2:end-1);
if any(inner(AkvMean) <= 0, 'all') || any(inner(AktMean) <= 0, 'all')
    error('Parent mean Akv/Akt is zero on an interior w-level.');
end
fprintf('  mean Akv on interior faces [%.3g, %.3g], Akt [%.3g, %.3g] m^2/s (%s)\n', ...
    min(inner(AkvMean), [], 'all'), max(inner(AkvMean), [], 'all'), ...
    min(inner(AktMean), [], 'all'), max(inner(AktMean), [], 'all'), ...
    ternary(opts.foldMixing, 'folded across the crest', 'not folded'));


%% Experiment folder

if exist(expDir, 'dir')
    rmdir(expDir, 's');   % opts.overwrite, checked above
end
mkdir(inDir);


%% Input files: the parent's, cut to Ny along-ridge rows

newLen = containers.Map({'eta_rho', 'eta_u', 'eta_v', 'eta_psi'}, ...
    {Ny + 2, Ny + 2, Ny + 1, Ny + 1});

gridFile = fullfile(inDir, 'neptune_grid.nc');
sliceNc(parentGrid, gridFile, newLen);
dy = 1 / mean(ncread(gridFile, 'pn'), 'all');
ncwrite(gridFile, 'el', Ny * dy);

nccreate(gridFile, 'Akv_pres', 'Dimensions', {'xi_rho', LLm + 2, 'eta_rho', Ny + 2, 's_w', N + 1});
ncwriteatt(gridFile, 'Akv_pres', 'long_name', 'prescribed vertical viscosity at w-points');
ncwriteatt(gridFile, 'Akv_pres', 'units', 'meter2 second-1');
ncwrite(gridFile, 'Akv_pres', repmat(AkvMean, [1, Ny + 2, 1]));
nccreate(gridFile, 'Akt_pres', 'Dimensions', {'xi_rho', LLm + 2, 'eta_rho', Ny + 2, 's_w', N + 1});
ncwriteatt(gridFile, 'Akt_pres', 'long_name', 'prescribed tracer vertical diffusivity at w-points');
ncwriteatt(gridFile, 'Akt_pres', 'units', 'meter2 second-1');
ncwrite(gridFile, 'Akt_pres', repmat(AktMean, [1, Ny + 2, 1]));
ncwriteatt(gridFile, '/', 'prescribed_mixing_parent', parentName);
ncwriteatt(gridFile, '/', 'prescribed_mixing_avg_files', strjoin(avgFiles, ' '));
ncwriteatt(gridFile, '/', 'prescribed_mixing_folded', double(opts.foldMixing));
ncwriteatt(gridFile, '/', 'prescribed_mixing_surface_cutoff_m', opts.surfaceMixingDepth);
ncwriteatt(gridFile, '/', 'prescribed_mixing_background_Akv_Akt', [opts.bgAkv opts.bgAkt]);

sliceNc(fullfile(parentIn, 'neptune_init.nc'), fullfile(inDir, 'neptune_init.nc'), newLen);
sliceNc(fullfile(parentIn, 'neptune_flux.nc'), fullfile(inDir, 'neptune_flux.nc'), newLen);

% No external forcing: zero wind stress on the parent's frc schema, held
% constant by two records (the same layout as the zero-flux file).
frcInfo = ncinfo(fullfile(parentIn, 'neptune_frc.nc'));
frcLen = [newLen.keys, {'frc_time'}];
frcVal = [newLen.values, {2}];
writeZeroForcing(frcInfo, fullfile(inDir, 'neptune_frc.nc'), containers.Map(frcLen, frcVal));
fprintf('  input files written (grid, init, zero-stress frc, flux)\n');


%% Configuration: the parent's files with the minimal 2-D edits

copyfile(fullfile(parentDir, 'Makefile'), expDir);
optFiles = dir(fullfile(parentDir, '*.opt'));
for f = {optFiles.name}
    if ~ismember(f{1}, {'cppdefs.opt', 'param.opt'})
        copyfile(fullfile(parentDir, f{1}), expDir);
    end
end

% Production diagnostics. diag_norm must be true for the KE budget
% (energy_budget.F stops otherwise), so momentum terms come out
% dz-normalized.
if opts.productionDiags
    fprintf('  production diagnostics:\n');
    prod = { ...
        'diagnostics.opt',   'diag_norm',         '.true.'; ...
        'diagnostics.opt',   'diag_avg',          '.true.'; ...
        'diagnostics.opt',   'diag_uv',           '.true.'; ...
        'diagnostics.opt',   'diag_trc',          '.true.'; ...
        'ocean_vars.opt',    'wrt_file_avg',      '.true.'; ...
        'ocean_vars.opt',    'output_period_his', '86400'; ...
        'ocean_vars.opt',    'nrpf_his',          '370'; ...
        'energy_budget.opt', 'do_energy_budget',  '.true.'; ...
        'energy_budget.opt', 'output_period',     '3.16224e+07'};
    for k = 1:size(prod, 1)
        setOptValue(fullfile(expDir, prod{k, 1}), prod{k, 2}, prod{k, 3});
    end
end

% param.opt: Ny interior along-ridge points, 1-D MPI decomposition in xi
paramNew = regexprep(paramTxt, '(MMm\s*=\s*)\d+', sprintf('$1%d', Ny));
paramNew = regexprep(paramNew, '(NP_XI\s*=\s*)\d+', sprintf('$1%d', opts.npXi));
paramNew = regexprep(paramNew, '(NP_ETA\s*=\s*)\d+', '$11');
assert(parseInt(paramNew, 'MMm') == Ny && parseInt(paramNew, 'NP_XI') == opts.npXi ...
    && parseInt(paramNew, 'NP_ETA') == 1 && parseInt(paramNew, 'LLm') == LLm, ...
    'param.opt edit failed');
writeText(fullfile(expDir, 'param.opt'), paramNew);

% cppdefs.opt: KPP family off, prescribed mixing on
cpp = splitlines(fileread(fullfile(parentDir, 'cppdefs.opt')));
isLmd = ~cellfun(@isempty, regexp(cpp, '^\s*#\s*(define|undef)\s+LMD_\w+', 'once'));
cpp = regexprep(cpp, '^(\s*#\s*)define(\s+LMD_\w+)', '$1undef$2');
insertAt = find(isLmd, 1, 'last');
if isempty(insertAt)
    error('No LMD_* flags found in the parent cppdefs.opt.');
end
cpp = [cpp(1:insertAt); ...
    {'        /*        - 2-D run: KPP off; Akv/Akt prescribed from the parent'; ...
     '                    3-D run''s mean mixing (ana_vmix.h, grid.F) */'; ...
     '# define ANA_VMIX'; '# define VMIX_PRESCRIBED_Huaiyu'}; cpp(insertAt+1:end)];
cppTxt = strjoin(cpp, newline);
assert(isempty(regexp(cppTxt, '(?m)^\s*#\s*define\s+LMD_', 'once')), 'LMD_* still defined');
writeText(fullfile(expDir, 'cppdefs.opt'), cppTxt);

% neptune.in: new title, and the vertical_mixing keyword that read_inp
% requires once LMD_MIXING is undefined (ana_vmix overwrites the values)
nin = splitlines(fileread(fullfile(parentDir, 'neptune.in')));
iTitle = find(startsWith(strtrim(nin), 'title:'), 1);
nin{iTitle + 1} = sprintf('   2-D (along-ridge-uniform) run, KPP off, Akv/Akt prescribed from %s', parentName);
if ~any(startsWith(strtrim(nin), 'vertical_mixing:'))
    vm = {'vertical_mixing:  Akv_bak  Akt_bak  [m^2/sec]'; '                  0.       0.'; ''};
    iAt = find(startsWith(strtrim(nin), 'v_sponge:'), 1);
    if isempty(iAt), iAt = numel(nin) + 1; end
    nin = [nin(1:iAt-1); vm; nin(iAt:end)];
end
writeText(fullfile(expDir, 'neptune.in'), strjoin(nin, newline));

% Fortran: the parent's overrides, with grid.F replaced by the generator's
% VMIX-capable version (identical to the parent's apart from the guarded
% VMIX_PRESCRIBED_Huaiyu blocks, checked here), plus the ana_vmix.h hook
genGrid = fileread(fullfile(script_dir, 'grid.F'));
parGridF = fullfile(parentDir, 'grid.F');
if exist(parGridF, 'file') && ~strcmp(stripVmixBlocks(genGrid), stripVmixBlocks(fileread(parGridF)))
    warning('NewExp_Neptune_ridge_2D:gridDiffers', ...
        ['The generator grid.F differs from the parent''s beyond the VMIX blocks; ', ...
         'the 2-D run uses the generator version.']);
end
parentF = dir(fullfile(parentDir, '*.F'));
for f = {parentF.name}
    if ~strcmp(f{1}, 'grid.F')
        copyfile(fullfile(parentDir, f{1}), expDir);
    end
end
copyfile(fullfile(script_dir, 'grid.F'), expDir);
copyfile(fullfile(script_dir, 'ana_vmix.h'), expDir);


%% Run, partition and join scripts (shared partition: npXi cores only)

oldDir = pwd;
restoreDir = onCleanup(@() cd(oldDir));
cd(expDir);
generate_run_script(1, config.slurm_account, config.slurm_email, 'run_roms_auto', ...
    opts.targetYear, opts.expname, opts.walltime, config.roms_root, opts.npXi);
generate_partit_script(opts.npXi, 1);
generate_ncjoin_script(config.slurm_account, config.slurm_email);
clear restoreDir


%% Provenance and final checks

writeText(fullfile(expDir, 'README_2D.txt'), sprintf([ ...
    '2-D ridge run generated %s by NewExp_Neptune_ridge_2D.m\n', ...
    'parent           : %s\n', ...
    'along-ridge Ny   : %d interior points (MMm), NP_XI = %d, NP_ETA = 1\n', ...
    'forcing          : zero wind stress; zero surface heat/salt flux (parent''s)\n', ...
    'vertical mixing  : KPP off; Akv/Akt = time- and y-mean of the parent''s\n', ...
    '                   avg files %s ... %s (%d files), folded = %d\n', ...
    'surface mixing   : above %g m replaced by background Akv %g / Akt %g (0 m = kept)\n', ...
    'diagnostics      : production set (budgets, avg file, daily his, KE budget) = %d\n', ...
    'unchanged        : grid, vertical coordinate, initial state, EOS, bottom\n', ...
    '                   drag, nudging strips, dt (the parent''s files)\n'], ...
    char(datetime('now', 'Format', 'yyyy-MM-dd')), parentName, Ny, opts.npXi, ...
    avgFiles{1}, avgFiles{end}, numel(avgFiles), opts.foldMixing, ...
    opts.surfaceMixingDepth, opts.bgAkv, opts.bgAkt, opts.productionDiags));

for f = {'neptune_grid.nc', 'neptune_init.nc', 'neptune_frc.nc', 'neptune_flux.nc'}
    assert_netcdf_no_nan(fullfile(inDir, f{1}));
end
fprintf('  done: %s\n  next: upload the folder, then on Expanse: sbatch run_roms_auto\n\n', expDir);

end


%% ------------------------------------------------------------------------

function [Akv, Akt] = removeSurfaceMixing(Akv, Akt, gridFile, opts)
% Replace Akv/Akt on w-levels shallower than opts.surfaceMixingDepth by the
% interior background. Depths use the ROMS set_depth transform at zeta = 0
% (sea level varies by centimetres, irrelevant for a ~1000 m cutoff):
%   z_w(k) = h (hc s_k + Cs_w(k) h) / (h + hc),  s_k = (k - N)/N.
D = opts.surfaceMixingDepth;
h = ncread(gridFile, 'h');   h = h(:, 1);
hc = ncread(gridFile, 'hc'); Cs = ncread(gridFile, 'Cs_w');
N = numel(Cs) - 1;
s = ((0:N) - N) / N;
depth = -(h .* (hc * s + h .* Cs(:)') ./ (h + hc));   % (xi, s_w), positive down
replace = depth < D;

% The cutoff must stay well clear of the bottom-enhanced layer
if min(h) - D < 500
    error('NewExp_Neptune_ridge_2D:surfaceCutoff', ...
        ['surfaceMixingDepth = %g m is within %g m of the shallowest seabed ', ...
         '(%.0f m); it would remove abyssal mixing. Use a shallower cutoff.'], ...
        D, min(h) - D, min(h));
end

% The parent's background, checked against the requested values
V = squeeze(Akv); T = squeeze(Akt);
vMin = min(V(:, 2:end-1), [], 'all'); tMin = min(T(:, 2:end-1), [], 'all');
if abs(vMin / opts.bgAkv - 1) > 0.01 || abs(tMin / opts.bgAkt - 1) > 0.01
    warning('NewExp_Neptune_ridge_2D:background', ...
        ['Parent interior minimum Akv %.3g / Akt %.3g differs from the requested ', ...
         'background %.3g / %.3g.'], vMin, tMin, opts.bgAkv, opts.bgAkt);
end

% Surface mixing should be fully contained above the cutoff: the first
% level below it should already be background, or the profile gets a step
kBelow = arrayfun(@(i) find(~replace(i, :), 1, 'last'), (1:size(V, 1))');
ratioT = T(sub2ind(size(T), (1:size(T, 1))', kBelow)) / opts.bgAkt;
ratioV = V(sub2ind(size(V), (1:size(V, 1))', kBelow)) / opts.bgAkv;
if max([ratioT; ratioV]) > 1.1
    warning('NewExp_Neptune_ridge_2D:surfaceStep', ...
        ['Surface-enhanced mixing reaches below %g m (just below the cutoff: ', ...
         'Akt up to %.2fx, Akv up to %.2fx background); the prescribed profile ', ...
         'steps there. Consider a deeper surfaceMixingDepth.'], D, max(ratioT), max(ratioV));
end

removedT = max(T(replace)); removedV = max(V(replace));
V(replace) = opts.bgAkv;
T(replace) = opts.bgAkt;
Akv = reshape(V, size(Akv));
Akt = reshape(T, size(Akt));
nRep = sum(replace, 2);
fprintf(['  surface mixing above %g m -> background Akv %.0e / Akt %.0e ', ...
    '(%d-%d top w-levels per column; removed max Akv %.3g, Akt %.3g; ', ...
    'just below: Akt <= %.2fx, Akv <= %.2fx background)\n'], ...
    D, opts.bgAkv, opts.bgAkt, min(nRep), max(nRep), removedV, removedT, max(ratioT), max(ratioV));
end

function opts = fillDefaults(opts, def)
unknown = setdiff(fieldnames(opts), fieldnames(def));
if ~isempty(unknown)
    error('NewExp_Neptune_ridge_2D:badOption', 'Unknown option(s): %s', strjoin(unknown, ', '));
end
for f = fieldnames(def)'
    if ~isfield(opts, f{1})
        opts.(f{1}) = def.(f{1});
    end
end
end

function mustExist(p, kind)
if exist(p, kind) == 0
    error('NewExp_Neptune_ridge_2D:missing', 'Missing %s: %s', kind, p);
end
end

function v = parseInt(txt, name)
tok = regexp(txt, ['(?<![A-Za-z_])' name '\s*=\s*(\d+)'], 'tokens', 'once');
if isempty(tok)
    error('Cannot find %s in param.opt.', name);
end
v = str2double(tok{1});
end

function M = yMean(A)
% Along-ridge mean over the interior eta rows (rows 1 and end are the
% periodic copies), keeping (xi, 1, s_w). The avg files store single
% precision; accumulate in double.
A = double(A);
A(abs(A) > 1e30) = NaN;
M = mean(A(:, 2:end-1, :, :), [2 4], 'omitnan');
end

function sliceNc(src, dst, newLen)
% Copy src to dst keeping only the first newLen(dim) entries of the
% listed dimensions (the parent is uniform along them).
I = resizeSchema(ncinfo(src), newLen);
ncwriteschema(dst, I);
for v = 1:numel(I.Variables)
    V = I.Variables(v);
    if isempty(V.Dimensions)
        ncwrite(dst, V.Name, ncread(src, V.Name));
    else
        ncwrite(dst, V.Name, ncread(src, V.Name, ones(1, numel(V.Size)), V.Size));
    end
end
end

function writeZeroForcing(I, dst, newLen)
% Parent forcing schema at the new sizes, two constant zero records.
I = resizeSchema(I, newLen);
for v = 1:numel(I.Variables)
    keep = ~strcmp({I.Variables(v).Attributes.Name}, 'cycle_length');
    I.Variables(v).Attributes = I.Variables(v).Attributes(keep);
end
ncwriteschema(dst, I);
for v = 1:numel(I.Variables)
    V = I.Variables(v);
    if strcmp(V.Name, 'frc_time')
        ncwrite(dst, V.Name, [0; 99999]);
    else
        ncwrite(dst, V.Name, zeros(V.Size));
    end
end
end

function I = resizeSchema(I, newLen)
for d = 1:numel(I.Dimensions)
    if isKey(newLen, I.Dimensions(d).Name)
        I.Dimensions(d).Length = newLen(I.Dimensions(d).Name);
    end
end
for v = 1:numel(I.Variables)
    for k = 1:numel(I.Variables(v).Dimensions)
        nm = I.Variables(v).Dimensions(k).Name;
        if isKey(newLen, nm)
            I.Variables(v).Dimensions(k).Length = newLen(nm);
            I.Variables(v).Size(k) = newLen(nm);
        end
    end
    I.Variables(v).ChunkSize = [];   % the parent's chunks may exceed the new sizes
end
end

function s = stripVmixBlocks(s)
s = regexprep(s, '#ifdef VMIX_PRESCRIBED_Huaiyu.*?#endif\r?\n', '');
end

function setOptValue(path, key, value)
% Set "<decl> :: key = value" in a Fortran .opt file, keeping the
% declaration and any trailing comment. The key must occur exactly once.
txt = fileread(path);
pat = ['(?m)^[^!\n]*::\s*' key '\s*=\s*([^!\s]+)'];
ext = regexp(txt, pat, 'tokenExtents');
if numel(ext) ~= 1
    error('NewExp_Neptune_ridge_2D:optKey', '%s: found "%s" %d times (expected 1).', ...
        path, key, numel(ext));
end
v = ext{1};
old = txt(v(1):v(2));
writeText(path, [txt(1:v(1)-1) value txt(v(2)+1:end)]);
[~, name, ext2] = fileparts(path);
fprintf('    %-18s %-18s %s -> %s\n', [name ext2], key, old, value);
end

function writeText(path, txt)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot write %s', path);
end
fprintf(fid, '%s', txt);
if ~endsWith(txt, newline)
    fprintf(fid, '\n');
end
fclose(fid);
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
