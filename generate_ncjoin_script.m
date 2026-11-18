function generate_ncjoin_script(account_name,user_name,filename)
%GENERATE_NCJOIN_SCRIPT Write the ncjoin SLURM script in the current directory.
% The script runs one ncjoin per output type (avg, dia, ekb) and one per
% calendar year of history files, up to MAX_JOIN_JOBS at a time.
% Default filename: 'do_ncjoin'

    if nargin < 1 || isempty(account_name)
        account_name = '';
    end

    if nargin < 2 || isempty(user_name)
        user_name = '';
    end

    if nargin < 3 || isempty(filename)
        filename = 'do_ncjoin';
    end

    [fid, cleaner] = open_text_file_for_write(filename);

    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#SBATCH --job-name="ncj_mpi50"\n\n');

    fprintf(fid, '#       %%j = job ID, %%N = short host name of the node\n');
    fprintf(fid, '#SBATCH --output="ncjoin_mpi50.%%j.%%N.log"\n\n');

    fprintf(fid, '#       If using less than 128 cores then the partition is ''shared''\n');
    fprintf(fid, '#       or ''large-shared'' if extra memory is needed (only use if ''shared'' failed due to memory error).\n');
    fprintf(fid, '#       Note, the job charge costs 4x more for ''large-shared'' than ''shared'' or ''compute'', use allocation wisely!\n');
    fprintf(fid, '#SBATCH --partition=shared\n');
    fprintf(fid, '#SBATCH --nodes=1\n\n');

    fprintf(fid, '#       Request number of cores: (Expanse has 128 cores per node)\n');
    fprintf(fid, '#       The ''compute'' partition charges 128 cores for every node requested,\n');
    fprintf(fid, '#       whatever --ntasks-per-node is.\n');
    fprintf(fid, '#       Up to this many ncjoin processes run at the same time (MAX_JOIN_JOBS).\n');
    fprintf(fid, '#SBATCH --ntasks-per-node=12\n\n');

    fprintf(fid, '#       Leave this at 1:\n');
    fprintf(fid, '#SBATCH --cpus-per-task=1\n\n');

    fprintf(fid, '#       Memory: default is 1GB on all nodes. However, you can request 2GB per core on\n');
    fprintf(fid, '#       shared/compute/debug or 15.5GB/core on large-shared at no extra cost.\n');
    fprintf(fid, '#       Charged on number of cores or fraction of total memory, whichever is greater.\n');
    fprintf(fid, '#       shared/compute/debug total mem=256G, large-shared total mem=2000G\n');
    fprintf(fid, '#            for total memory required use: #SBATCH --mem=256G\n');
    fprintf(fid, '#         or for memory per cores      use: #SBATCH --mem-per-cpu=2G\n');
    fprintf(fid, '#SBATCH --mem=64G\n');
    if ~isempty(account_name)
        fprintf(fid, '#SBATCH --account=%s\n', account_name);
    else
        fprintf(fid, '##SBATCH --account=<set-your-account-if-required>\n');
    end
    fprintf(fid, '#SBATCH --export=ALL\n');
    fprintf(fid, '#       Time duration requested for run:\n');
    fprintf(fid, '#SBATCH -t 24:00:00\n');
    if ~isempty(user_name)
        fprintf(fid, '#SBATCH --mail-user=%s\n', user_name);
        fprintf(fid, '#SBATCH --mail-type=ALL\n');
    else
        fprintf(fid, '# Email notifications disabled: no mail address is set (slurm_email in newexp_defaults.m).\n');
    end
    fprintf(fid, '\n');

    fprintf(fid, '#-----------------------------------------------------------------\n\n');
    fprintf(fid, 'set -euo pipefail\n\n');

    write_line(fid, 'MAX_JOIN_JOBS="${MAX_JOIN_JOBS:-${SLURM_NTASKS:-12}}"');
    write_line(fid, 'if ! [[ "${MAX_JOIN_JOBS}" =~ ^[1-9][0-9]*$ ]]; then');
    write_line(fid, '    echo "ERROR: MAX_JOIN_JOBS must be a positive integer, got ${MAX_JOIN_JOBS}." >&2');
    write_line(fid, '    exit 2');
    write_line(fid, 'fi');
    write_line(fid, 'join_failed=0');
    write_line(fid, 'join_pids=()');
    write_line(fid, 'join_labels=()');
    write_line(fid, '');

    fprintf(fid, '# Flags needed for mvapich2:\n');
    fprintf(fid, 'export MV2_USE_RDMA_CM=0\n');
    fprintf(fid, 'export MV2_IBA_HCA=mlx5_2\n');
    fprintf(fid, 'export MV2_DEFAULT_PORT=1\n\n');

    fprintf(fid, 'module purge\n');
    fprintf(fid, 'module load slurm\n');
    fprintf(fid, 'module load cpu/0.15.4  intel/19.1.1.217  mvapich2/2.3.4\n');
    fprintf(fid, 'module load netcdf-c/4.7.4\n');
    fprintf(fid, 'module load netcdf-fortran/4.5.3\n');
    fprintf(fid, '\n');
    fprintf(fid, 'if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then\n');
    fprintf(fid, '    cd "${SLURM_SUBMIT_DIR}"\n');
    fprintf(fid, 'fi\n');
    fprintf(fid, '\n');
    write_line(fid, '# SIM_RUNNING controls whether the newest output set of each type is left');
    write_line(fid, '# alone, because "ncjoin -d" deletes partial files and would kill the run');
    write_line(fid, '# if it touched the files ROMS is writing into.');
    write_line(fid, '#   SIM_RUNNING=auto : (default) detect via squeue -- the run counts as');
    write_line(fid, '#                      active if any other job of this user is running or');
    write_line(fid, '#                      queued with this directory as its working directory');
    write_line(fid, '#   SIM_RUNNING=1    : force "run is still active" (skip newest sets)');
    write_line(fid, '#   SIM_RUNNING=0    : force "run is finished" (join everything)');
    write_line(fid, '# Override at submit time:  sbatch --export=ALL,SIM_RUNNING=0 do_ncjoin');
    write_line(fid, 'SIM_RUNNING="${SIM_RUNNING:-auto}"');
    write_line(fid, '');
    write_line(fid, 'sim_is_running() {');
    write_line(fid, '    case "${SIM_RUNNING}" in');
    write_line(fid, '        1|yes|true)  return 0 ;;');
    write_line(fid, '        0|no|false)  return 1 ;;');
    write_line(fid, '    esac');
    write_line(fid, '');
    write_line(fid, '    if ! command -v squeue >/dev/null 2>&1; then');
    write_line(fid, '        echo "WARNING: squeue not found; assuming the simulation is still running." >&2');
    write_line(fid, '        return 0');
    write_line(fid, '    fi');
    write_line(fid, '');
    write_line(fid, '    local here jobid workdir');
    write_line(fid, '    here="$(pwd -P)"');
    write_line(fid, '    while read -r jobid workdir; do');
    write_line(fid, '        [[ -n "${jobid}" ]] || continue');
    write_line(fid, '        [[ "${jobid}" == "${SLURM_JOB_ID:-}" ]] && continue');
    write_line(fid, '        workdir="$(readlink -f "${workdir}" 2>/dev/null || printf ''%s'' "${workdir}")"');
    write_line(fid, '        if [[ "${workdir}" == "${here}" ]]; then');
    write_line(fid, '            echo "Detected active/queued job ${jobid} in this directory."');
    write_line(fid, '            return 0');
    write_line(fid, '        fi');
    write_line(fid, '    done < <(squeue -h -u "${USER:-$(whoami)}" -o ''%A %Z'' 2>/dev/null)');
    write_line(fid, '');
    write_line(fid, '    return 1');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, '# Timestamp of a partial file: second dot-field of neptune_<type>.<stamp>.<rank>.nc');
    write_line(fid, 'stamp_of() {');
    write_line(fid, '    local stamp="${1#*.}"');
    write_line(fid, '    printf ''%s\n'' "${stamp%%.*}"');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, '# Newest (largest) timestamp among the partial files of one output type.');
    write_line(fid, 'newest_stamp() {');
    write_line(fid, '    local prefix="$1" file stamp newest=""');
    write_line(fid, '    shopt -s nullglob');
    write_line(fid, '    for file in ${prefix}.*.*.nc; do');
    write_line(fid, '        stamp="$(stamp_of "${file}")"');
    write_line(fid, '        if [[ "${stamp}" > "${newest}" ]]; then');
    write_line(fid, '            newest="${stamp}"');
    write_line(fid, '        fi');
    write_line(fid, '    done');
    write_line(fid, '    shopt -u nullglob');
    write_line(fid, '    printf ''%s\n'' "${newest}"');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'wait_for_join_jobs() {');
    write_line(fid, '    local i pid label status');
    write_line(fid, '    for ((i=0; i<${#join_pids[@]}; i++)); do');
    write_line(fid, '        pid="${join_pids[$i]}"');
    write_line(fid, '        label="${join_labels[$i]}"');
    write_line(fid, '        if wait "${pid}"; then');
    write_line(fid, '            echo "Completed ncjoin for ${label}."');
    write_line(fid, '        else');
    write_line(fid, '            status=$?');
    write_line(fid, '            echo "ERROR: ncjoin failed for ${label} with status ${status}." >&2');
    write_line(fid, '            join_failed=1');
    write_line(fid, '        fi');
    write_line(fid, '    done');
    write_line(fid, '    join_pids=()');
    write_line(fid, '    join_labels=()');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'join_if_present() {');
    write_line(fid, '    local pattern="$1"');
    write_line(fid, '    local exclude_stamp="${2:-}"');
    write_line(fid, '    shopt -s nullglob');
    write_line(fid, '    local candidates=( $pattern )');
    write_line(fid, '    shopt -u nullglob');
    write_line(fid, '');
    write_line(fid, '    local files=() file');
    write_line(fid, '    for file in "${candidates[@]}"; do');
    write_line(fid, '        if [[ -n "${exclude_stamp}" && "$(stamp_of "${file}")" == "${exclude_stamp}" ]]; then');
    write_line(fid, '            continue');
    write_line(fid, '        fi');
    write_line(fid, '        files+=( "${file}" )');
    write_line(fid, '    done');
    write_line(fid, '');
    write_line(fid, '    if (( ${#files[@]} < ${#candidates[@]} )); then');
    write_line(fid, '        echo "Leaving in-progress set ${pattern%%.*}.${exclude_stamp}.*.nc alone (simulation still running)."');
    write_line(fid, '    fi');
    write_line(fid, '');
    write_line(fid, '    if (( ${#files[@]} == 0 )); then');
    write_line(fid, '        echo "Skipping ${pattern}: no joinable files."');
    write_line(fid, '        return 0');
    write_line(fid, '    fi');
    write_line(fid, '');
    write_line(fid, '    ncjoin -d "${files[@]}" &');
    write_line(fid, '    join_pids+=("$!")');
    write_line(fid, '    join_labels+=("${pattern}")');
    write_line(fid, '    if (( ${#join_pids[@]} >= MAX_JOIN_JOBS )); then');
    write_line(fid, '        wait_for_join_jobs');
    write_line(fid, '    fi');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'avg_exclude=""');
    write_line(fid, 'dia_exclude=""');
    write_line(fid, 'ekb_exclude=""');
    write_line(fid, 'his_exclude=""');
    write_line(fid, 'if sim_is_running; then');
    write_line(fid, '    avg_exclude="$(newest_stamp neptune_avg)"');
    write_line(fid, '    dia_exclude="$(newest_stamp neptune_dia)"');
    write_line(fid, '    ekb_exclude="$(newest_stamp neptune_ekb)"');
    write_line(fid, '    his_exclude="$(newest_stamp neptune_his)"');
    write_line(fid, '    echo "Simulation still active (SIM_RUNNING=${SIM_RUNNING}): the newest output set of each type will not be joined."');
    write_line(fid, 'else');
    write_line(fid, '    echo "No active simulation detected (SIM_RUNNING=${SIM_RUNNING}): joining all output sets."');
    write_line(fid, 'fi');
    write_line(fid, '');
    write_line(fid, 'join_if_present "neptune_avg.*.*.nc" "${avg_exclude}"');
    write_line(fid, 'join_if_present "neptune_dia.*.*.nc" "${dia_exclude}"');
    write_line(fid, 'join_if_present "neptune_ekb.*.*.nc" "${ekb_exclude}"');
    write_line(fid, 'shopt -s nullglob');
    write_line(fid, 'declare -A history_years=()');
    write_line(fid, 'for history_file in neptune_his.*.*.nc; do');
    write_line(fid, '    stamp="$(stamp_of "${history_file}")"');
    write_line(fid, '    year="${stamp:0:4}"');
    write_line(fid, '    if [[ "${year}" =~ ^[0-9]{4}$ ]]; then');
    write_line(fid, '        history_years["${year}"]=1');
    write_line(fid, '    fi');
    write_line(fid, 'done');
    write_line(fid, 'shopt -u nullglob');
    write_line(fid, '');
    write_line(fid, 'for year in "${!history_years[@]}"; do');
    write_line(fid, '    join_if_present "neptune_his.${year}*.*.nc" "${his_exclude}"');
    write_line(fid, 'done');
    write_line(fid, 'wait_for_join_jobs');
    write_line(fid, 'if (( join_failed != 0 )); then');
    write_line(fid, '    echo "ERROR: one or more ncjoin operations failed." >&2');
    write_line(fid, '    exit 1');
    write_line(fid, 'fi');
    clear cleaner;

    if isunix || ismac
        fileattrib(filename, '+x', 'a');
    end

    fprintf('Wrote "%s" in the current directory.\n', filename);
end

function write_line(fid, line)
    fprintf(fid, '%s\n', line);
end
