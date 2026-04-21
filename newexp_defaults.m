function config = newexp_defaults(script_path)
%NEWEXP_DEFAULTS Shared defaults for the public NewExp_Expanse scripts.
%   CONFIG = NEWEXP_DEFAULTS(SCRIPT_PATH) returns a struct with portable
%   path defaults and optional cluster metadata used by the experiment
%   generators in this folder.

    if nargin < 1 || isempty(script_path)
        script_path = mfilename('fullpath');
    end

    script_dir = fileparts(script_path);
    project_root = fileparts(script_dir);

    config = struct();
    config.script_dir = script_dir;
    config.project_root = project_root;

    % Prefer the original local experiment root when it is available, and
    % fall back to a package-local folder for portable public use.
    config.preferred_output_root = 'E:\Data_ROMS\Neptune';
    config.fallback_output_root = fullfile(script_dir, 'generated_experiments');

    % Optional helper/toolbox paths. Leave roms_tools_dir empty if those
    % utilities are already on your MATLAB path or are not needed.
    config.roms_tools_dir = '';
    config.postprocess_dir = fullfile(project_root, 'PostProcess_2026');

    % Optional SLURM metadata. Leave blank to omit those directives in the
    % generated submission scripts and edit them manually on your cluster.
    config.slurm_account = 'cla327';
    config.slurm_email = 'hwei1';

    % Resolve the actual output root now so callers do not need a separate
    % helper script for this decision.
    config.output_root = resolve_output_root( ...
        config.preferred_output_root, config.fallback_output_root);
end

function output_root = resolve_output_root(preferred_output_root, fallback_output_root)
    candidates = {preferred_output_root, fallback_output_root};

    for i = 1:numel(candidates)
        candidate = candidates{i};
        if isempty(candidate)
            continue
        end

        if exist(candidate, 'dir') ~= 7
            [ok, msg] = mkdir(candidate);
            if ~ok
                warning('Could not create output root "%s": %s', candidate, msg);
                continue
            end
        end

        if ~is_directory_writable(candidate)
            warning('Output root "%s" exists but is not writable.', candidate);
            continue
        end

        output_root = candidate;
        if i > 1
            warning('Preferred output root "%s" is unavailable; using "%s" instead.', ...
                preferred_output_root, output_root);
        end
        return
    end

    error('Could not create or access any configured experiment output root.');
end

function tf = is_directory_writable(directory)
    test_file = tempname(directory);
    fid = fopen(test_file, 'w');
    if fid == -1
        tf = false;
        return
    end

    fclose(fid);
    delete(test_file);
    tf = true;
end
