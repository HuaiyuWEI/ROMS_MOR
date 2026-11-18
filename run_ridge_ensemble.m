function run_ridge_ensemble(mode, subset, skip_existing)
% RUN_RIDGE_ENSEMBLE  Generate the ridge perturbation-experiment set.
%
%   run_ridge_ensemble()                   list the planned experiments (no generation)
%   run_ridge_ensemble('go')               generate all of them
%   run_ridge_ensemble('go', [2 5])        generate only these rows of the listing
%   run_ridge_ensemble('go', [], true)     skip experiments whose folder already
%                                          exists (e.g. to resume after a failure)
%
% One-at-a-time perturbations about the reference values below: RF_F0_rot
% (weak/strong forcing), h_ridge (short/tall ridge), zob (low/high Cd),
% alpha0 (low/high N2 via alphaT).
% Every experiment name is built from its perturbed values, so names and
% parameters cannot disagree. All four perturbable parameters are passed
% to NewExp_Neptune_ridge explicitly for every run, so the runs never
% depend on the defaults inside the generator.
%
% By default existing experiment folders are REGENERATED (input files are
% overwritten in place). If a regenerated experiment also exists on the
% cluster, delete its Neptune_input/*.<node>.nc partition files there so
% run_roms_auto repartitions the new inputs instead of reusing stale ones.
% Everything else about a run (forcing type, nudging, convection,
% resolution, time step, ...) comes from NewExp_Neptune_ridge.

if nargin < 1 || isempty(mode)
    mode = 'list';
end
if nargin < 2
    subset = [];
end
if nargin < 3 || isempty(skip_existing)
    skip_existing = false;
end
if ~any(strcmpi(mode, {'list', 'go'}))
    error('run_ridge_ensemble:badMode', ...
        'mode must be ''list'' or ''go'' (got ''%s'')', mode);
end
validateattributes(skip_existing, {'logical', 'numeric'}, {'scalar', 'binary'}, ...
    mfilename, 'skip_existing');
skip_existing = logical(skip_existing);
doGenerate = strcmpi(mode, 'go');

%% Reference values of the perturbable parameters

ref = struct();
ref.h_ridge   = 2000;    % ridge height (m)
ref.RF_F0_rot = 0.75;    % random-forcing amplitude
ref.alpha0    = 2e-4;    % thermal expansion coefficient (= 0.2/rho0, rho0 = 1000)
ref.zob       = 1e-3;    % bottom roughness (m)

% Optional qualifier appended to every experiment name (e.g. '_v2').
name_suffix = '';

%% Perturbation table: {tag, parameter, value}; empty parameter = reference run

perturbations = { ...
    'ref',           '',           [];    ...
    'weakforcing',   'RF_F0_rot',  0.25;  ...
    'strongforcing', 'RF_F0_rot',  1.25;  ...
    'shortridge',    'h_ridge',    1000;  ...
    'tallridge',     'h_ridge',    3000;  ...
    'highCd',        'zob',        1e-2;  ...
    'lowCd',         'zob',        1e-4;  ...
    'lowN2',         'alpha0',     1e-4;  ...
    'highN2',        'alpha0',     3e-4;  ...
    };

%% Build the run list

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
config = newexp_defaults(fullfile(script_dir, 'NewExp_Neptune_ridge.m'));

nRuns = size(perturbations, 1);
runs = struct('expname', cell(nRuns, 1), 'tag', [], 'overrides', [], 'exists', []);
for k = 1:nRuns
    [tag, param, value] = perturbations{k, :};
    p = ref;
    if ~isempty(param)
        if ~isfield(ref, param)
            error('run_ridge_ensemble:badParam', ...
                'Unknown parameter ''%s'' in perturbation row %d.', param, k);
        end
        p.(param) = value;
    end
    ovr = p;
    ovr.expname = make_expname(p, [tag name_suffix]);
    runs(k).tag = tag;
    runs(k).overrides = ovr;
    runs(k).expname = ovr.expname;
    runs(k).exists = exist(fullfile(config.output_root, ovr.expname), 'dir') == 7;
end

[~, iUnique] = unique({runs.expname});
if numel(iUnique) < nRuns
    dup = setdiff(1:nRuns, iUnique);
    error('run_ridge_ensemble:duplicateName', ...
        'Perturbation rows produce a duplicate experiment name: %s', ...
        runs(dup(1)).expname);
end

if isempty(subset)
    subset = 1:nRuns;
end
validateattributes(subset, {'numeric'}, ...
    {'vector', 'integer', 'positive', '<=', nRuns}, mfilename, 'subset');

%% Report the plan

fprintf('\nOutput root: %s\n\n', config.output_root);
fprintf('%-3s %-14s %-52s %s\n', '#', 'tag', 'expname', 'status');
for k = 1:nRuns
    if ~ismember(k, subset)
        status = 'not selected';
    elseif runs(k).exists && skip_existing
        status = 'EXISTS (will skip)';
    elseif runs(k).exists
        status = 'EXISTS (will REGENERATE in place)';
    else
        status = 'planned';
    end
    fprintf('%-3d %-14s %-52s %s\n', k, runs(k).tag, runs(k).expname, status);
end

if ~doGenerate
    fprintf('\nDry run only -- call run_ridge_ensemble(''go'') to generate.\n');
    return
end

%% Generate

start_dir = pwd;
restore_dir = onCleanup(@() cd(start_dir));   % generator cds into each experiment

nDone = 0;
nSkip = 0;
nRegen = 0;
for k = subset(:)'
    if runs(k).exists && skip_existing
        fprintf('\n[%d/%d] %s already exists -- skipping.\n', k, nRuns, runs(k).expname);
        nSkip = nSkip + 1;
        continue
    end
    if runs(k).exists
        fprintf('\n===== [%d/%d] REGENERATING %s (folder exists; input files overwritten) =====\n', ...
            k, nRuns, runs(k).expname);
        nRegen = nRegen + 1;
    else
        fprintf('\n===== [%d/%d] generating %s =====\n', k, nRuns, runs(k).expname);
    end
    t0 = tic;
    try
        NewExp_Neptune_ridge(runs(k).overrides);
    catch err
        cd(start_dir);
        fprintf(2, ['FAILED on %s after %.1f min ', ...
            '(%d generated, %d skipped before the failure).\n'], ...
            runs(k).expname, toc(t0)/60, nDone, nSkip);
        rethrow(err);
    end
    cd(start_dir);
    close all
    nDone = nDone + 1;
    fprintf('===== done: %s (%.1f min) =====\n', runs(k).expname, toc(t0)/60);
end

fprintf('\nEnsemble complete: %d generated (%d of them regenerated over an existing folder), %d skipped.\n', ...
    nDone, nRegen, nSkip);
fprintf('Output root: %s\n', config.output_root);
fprintf('Next: upload the experiments to Expanse and submit each one''s run_roms_auto.\n');
if nRegen > 0
    fprintf(['NOTE: for regenerated experiments that already exist on the cluster, delete\n' ...
        'their Neptune_input/*.<node>.nc partition files there -- run_roms_auto skips\n' ...
        'repartitioning when partitions are present and would run on the STALE inputs.\n']);
end

end


function name = make_expname(p, tag)
% Ridge_KPPconv_hr{H}_tau{100*RF}_alphaT{m}em4_zob{m}em{e}_{tag}
name = sprintf('Ridge_KPPconv_hr%d_tau%03d_alphaT%sem4_zob%s_%s', ...
    round(p.h_ridge), round(100*p.RF_F0_rot), ...
    num_token(p.alpha0*1e4), zob_token(p.zob), tag);
end


function tok = zob_token(zob)
% zob as normalized mantissa/exponent: 1e-3 -> '1em3', 3e-3 -> '3em3'
e = ceil(-log10(zob) - 1e-9);
m = round(zob * 10^e, 6);
if m >= 10
    m = m / 10;
    e = e - 1;
end
tok = sprintf('%sem%d', num_token(m), e);
end


function tok = num_token(x)
% Compact numeric token: integers as-is, decimal point as 'p' (0.5 -> '0p5')
x = round(x, 6);
if abs(x - round(x)) < 1e-9
    tok = sprintf('%d', round(x));
else
    tok = strrep(sprintf('%g', x), '.', 'p');
end
end
