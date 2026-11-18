function config = newexp_defaults(script_path)
%NEWEXP_DEFAULTS Settings shared by the experiment generators.
%   CONFIG = NEWEXP_DEFAULTS(SCRIPT_PATH) returns a struct with the output
%   roots, optional helper paths and cluster settings. SCRIPT_PATH is the
%   full path of the calling script (default: this file); the fallback
%   output root is placed in its folder. The values are the author's;
%   change them before use.

    if nargin < 1 || isempty(script_path)
        script_path = mfilename('fullpath');
    end

    script_dir = fileparts(script_path);
    project_root = fileparts(script_dir);

    config = struct();
    config.script_dir = script_dir;
    config.project_root = project_root;

    % Experiment directories go to preferred_output_root and, when that
    % cannot be created or is not writable, to generated_experiments/ in
    % this folder.
    config.preferred_output_root = 'E:\Data_ROMS\Neptune_ridge';
    config.fallback_output_root = fullfile(script_dir, 'generated_experiments');

    % Optional folders added to the MATLAB path; empty when not needed.
    config.roms_tools_dir = '';
    config.postprocess_dir = fullfile(project_root, 'analysis_ridge');

    % SLURM account and mail address; empty values leave those directives
    % out of the generated scripts.
    config.slurm_account = 'cla327';
    config.slurm_email = 'hwei1';

    % Path of the ucla-roms (development branch) clone on the cluster; the
    % generated run script exports it as ROMS_ROOT before compiling.
    config.roms_root = '~/work2026/ucla-roms-dev';

    % Take the first root that exists or can be created, and is writable.
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

        % On Unix, a Windows drive-letter path is a valid relative name;
        % skip it explicitly instead of creating a literal "E:\..." folder.
        if ~ispc && ~isempty(regexp(candidate, '^[A-Za-z]:[\\/]', 'once'))
            warning('Output root "%s" is not valid on this platform.', candidate);
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
