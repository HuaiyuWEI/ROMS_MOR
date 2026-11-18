function generate_run_script(num_nodes, account_name, user_name, filename, target_year, expname, walltime, roms_root, shared_ntasks)
    % Generate the auto-restarting Expanse SLURM batch script
    % (development-branch ROMS: exports ROMS_ROOT/FCOMP/NETCDF_ROOT before
    % compiling and pre-checks that the NHMG library has been built).
    %
    %   generate_run_script(num_nodes, account_name, user_name, filename, ...
    %       target_year, expname, walltime, roms_root, shared_ntasks)
    %
    % Only num_nodes is required. filename defaults to 'run_roms_auto', and
    % target_year, the calendar year at which the resubmission stops, to
    % 2014. roms_root is the path that the script exports as ROMS_ROOT.
    % shared_ntasks (fewer than 128) runs the job on the 'shared' partition
    % with that many MPI tasks instead of whole 'compute' nodes.
    %
    % When expname is supplied, the SLURM job name is taken from its trailing
    % '_'-delimited segment (e.g. 'Ridge_..._highN2' -> job-name="highN2").
    %
    % walltime is the requested SLURM wall time per job segment: either a
    % numeric hour count (36 -> '36:00:00') or a SLURM time string such as
    % '24:00:00' or '2-00:00:00'. Default '48:00:00' (the Expanse compute
    % partition maximum).

    if nargin < 4 || isempty(filename)
        filename = 'run_roms_auto';
    end

    if nargin < 7 || isempty(walltime)
        walltime = '48:00:00';
    end
    walltime_text = format_walltime(walltime);

    if nargin < 8 || isempty(roms_root)
        roms_root = '~/work2026/ucla-roms-dev';   % default: the author's clone (see newexp_defaults.m)
    end

    if nargin < 6 || isempty(expname)
        expname = '';
    end

    if nargin < 2 || isempty(account_name)
        account_name = '';
    end

    if nargin < 3 || isempty(user_name)
        user_name = '';
    end

    if nargin < 5 || isempty(target_year)
        target_year = 2014;
    end
    validateattributes(num_nodes, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'num_nodes');

    if isnumeric(target_year)
        target_year_text = sprintf('%04d', target_year);
    else
        target_year_text = char(target_year);
    end

    if isempty(regexp(target_year_text, '^[0-9]{4}$', 'once'))
        error('target_year must be a four-digit year.');
    end

    num_tasks = num_nodes * 128;

    % Optional 9th argument: run on the 'shared' partition with this many
    % MPI tasks (< 128) instead of whole 128-core 'compute' nodes. Used by
    % small jobs such as the 2-D ridge runs; omit it for the 3-D runs.
    if nargin < 9 || isempty(shared_ntasks)
        shared_ntasks = [];
    else
        validateattributes(shared_ntasks, {'numeric'}, ...
            {'scalar', 'integer', 'positive', '<', 128}, mfilename, 'shared_ntasks');
        num_tasks = shared_ntasks;
    end

    job_name = job_name_from_expname(expname);

    [fid, cleaner] = open_text_file_for_write(filename); %#ok<ASGLU>

    write_line(fid, '#!/bin/bash');
    fprintf(fid, '#SBATCH --job-name="%s"\n', job_name);
    write_line(fid, '');
    write_line(fid, '# %j = job ID, %N = short host name of the node');
    write_line(fid, '#SBATCH --output="log.%j.%N.log"');
    write_line(fid, '');
    if isempty(shared_ntasks)
        write_line(fid, '# Can only use a max of 32 nodes on ''compute'' partition:');
        write_line(fid, '#SBATCH --partition=compute');
        fprintf(fid, '#SBATCH --nodes=%d\n\n', num_nodes);
        write_line(fid, '# Request number of cores (Expanse has 128 cores per node)');
        write_line(fid, '#SBATCH --ntasks-per-node=128');
    else
        write_line(fid, '# Fewer than 128 cores: the ''shared'' partition bills only the cores');
        write_line(fid, '# requested (''compute'' would charge the whole 128-core node).');
        write_line(fid, '#SBATCH --partition=shared');
        write_line(fid, '#SBATCH --nodes=1');
        fprintf(fid, '#SBATCH --ntasks=%d\n', shared_ntasks);
        write_line(fid, '#SBATCH --mem-per-cpu=2G');
    end
    if ~isempty(account_name)
        fprintf(fid, '#SBATCH --account=%s\n', account_name);
    else
        write_line(fid, '##SBATCH --account=<set-your-account-if-required>');
    end
    write_line(fid, '#SBATCH --export=ALL');
    write_line(fid, '');
    write_line(fid, '# Requested wall time per job segment (Expanse compute allows at most 48 h).');
    write_line(fid, '# Ask Slurm for a signal before the limit so the next job can be queued even');
    write_line(fid, '# if this run is about to be stopped.');
    fprintf(fid, '#SBATCH -t %s\n', walltime_text);
    write_line(fid, '#SBATCH --signal=B:USR1@600');
    write_line(fid, '');
    write_line(fid, '# Email notifications:');
    if ~isempty(user_name)
        fprintf(fid, '#SBATCH --mail-user=%s\n', user_name);
        write_line(fid, '#SBATCH --mail-type=ALL');
    else
        write_line(fid, '# Email notifications disabled: no mail address is set (slurm_email in newexp_defaults.m).');
    end
    write_line(fid, '');
    write_line(fid, '#-----------------------------------------------------------------');
    write_line(fid, '');
    write_line(fid, 'set -euo pipefail');
    write_line(fid, 'export MV2_SMP_USE_CMA=0');
    write_line(fid, '');
    write_line(fid, '# Stop once the newest neptune_rst.YYYYMMDDhhmmss[.rank].nc restart set has');
    write_line(fid, '# reached this calendar year. Override at submit time, for example:');
    write_line(fid, '#   sbatch --export=ALL,TARGET_YEAR=2035 run_roms_auto');
    fprintf(fid, 'TARGET_YEAR="${TARGET_YEAR:-%s}"\n', target_year_text);
    write_line(fid, '');
    write_line(fid, 'INPUT_FILE="${INPUT_FILE:-neptune.in}"');
    write_line(fid, 'RESTART_PREFIX="${RESTART_PREFIX:-neptune_rst}"');
    fprintf(fid, 'NTASKS="${SLURM_NTASKS:-%d}"\n', num_tasks);
    write_line(fid, 'SCRIPT_NAME="${SCRIPT_NAME:-run_roms_auto}"');
    write_line(fid, 'BACKUP_DIR="${BACKUP_DIR:-backup}"');
    write_line(fid, 'NCJOIN_SCRIPT="${NCJOIN_SCRIPT:-do_ncjoin}"');
    write_line(fid, '');
    write_line(fid, 'if ! [[ "${TARGET_YEAR}" =~ ^[0-9]{4}$ ]]; then');
    write_line(fid, '  echo "ERROR: TARGET_YEAR must be a four-digit year, got ''${TARGET_YEAR}''." >&2');
    write_line(fid, '  exit 2');
    write_line(fid, 'fi');
    write_line(fid, '');
    write_line(fid, 'source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"');
    write_line(fid, 'experiment_dir=""');
    write_line(fid, '');
    write_line(fid, 'for candidate_dir in "${SLURM_SUBMIT_DIR:-}" "${PWD}" "${source_dir}"; do');
    write_line(fid, '  [[ -n "${candidate_dir}" ]] || continue');
    write_line(fid, '  if [[ -f "${candidate_dir}/${INPUT_FILE}" && -d "${candidate_dir}/Neptune_input" ]]; then');
    write_line(fid, '    experiment_dir="$(cd "${candidate_dir}" && pwd -P)"');
    write_line(fid, '    break');
    write_line(fid, '  fi');
    write_line(fid, 'done');
    write_line(fid, '');
    write_line(fid, 'if [[ -z "${experiment_dir}" ]]; then');
    write_line(fid, '  echo "ERROR: could not locate the experiment directory." >&2');
    write_line(fid, '  echo "       Expected to find ${INPUT_FILE} and Neptune_input/." >&2');
    write_line(fid, '  echo "       SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}" >&2');
    write_line(fid, '  echo "       PWD=${PWD}" >&2');
    write_line(fid, '  echo "       script source dir=${source_dir}" >&2');
    write_line(fid, '  exit 2');
    write_line(fid, 'fi');
    write_line(fid, '');
    write_line(fid, 'SCRIPT_DIR="${experiment_dir}"');
    write_line(fid, 'SCRIPT_PATH="${SCRIPT_DIR}/${SCRIPT_NAME}"');
    write_line(fid, '');
    write_line(fid, 'if [[ ! -f "${SCRIPT_PATH}" ]]; then');
    write_line(fid, '  echo "ERROR: cannot find ${SCRIPT_PATH} for automatic resubmission." >&2');
    write_line(fid, '  echo "       Submit this job from the experiment directory, or set SCRIPT_NAME." >&2');
    write_line(fid, '  exit 2');
    write_line(fid, 'fi');
    write_line(fid, '');
    write_line(fid, 'cd "${SCRIPT_DIR}"');
    write_line(fid, '');
    write_line(fid, '# A run that blew up leaves this sentinel; no job of the chain runs');
    write_line(fid, '# while the file exists.');
    write_line(fid, 'if [[ -f BLOWUP_DETECTED ]]; then');
    write_line(fid, '  echo "ERROR: BLOWUP_DETECTED sentinel found in ${SCRIPT_DIR}." >&2');
    write_line(fid, '  echo "       A previous run blew up. Fix the setup (e.g. dt), then" >&2');
    write_line(fid, '  echo "       re-enable runs with: rm BLOWUP_DETECTED" >&2');
    write_line(fid, '  exit 3');
    write_line(fid, 'fi');
    write_line(fid, '');
    write_line(fid, 'latest_restart_file() {');
    write_line(fid, '  local latest=""');
    write_line(fid, '  local latest_index=""');
    write_line(fid, '  local rst base rest index suffix rank');
    write_line(fid, '');
    write_line(fid, '  shopt -s nullglob');
    write_line(fid, '  for rst in ./"${RESTART_PREFIX}".*.nc; do');
    write_line(fid, '    base="${rst##*/}"');
    write_line(fid, '    rest="${base#${RESTART_PREFIX}.}"');
    write_line(fid, '    index="${rest%%.*}"');
    write_line(fid, '    suffix="${rest#${index}.}"');
    write_line(fid, '');
    write_line(fid, '    [[ "${index}" =~ ^[0-9]{14}$ ]] || continue');
    write_line(fid, '    if [[ "${suffix}" != "nc" ]]; then');
    write_line(fid, '      rank="${suffix%.nc}"');
    write_line(fid, '      # rank suffix of ANY width (ROMS uses the digit count of the rank count:');
    write_line(fid, '      # .00 for 20 ranks, .000 for 128). Never skip a restart set here: a');
    write_line(fid, '      # missed restart would silently cold-start from the initial file.');
    write_line(fid, '      [[ "${rank}" =~ ^[0-9]+$ ]] || continue');
    write_line(fid, '    fi');
    write_line(fid, '');
    write_line(fid, '    if [[ -z "${latest_index}" || "${index}" > "${latest_index}" ]]; then');
    write_line(fid, '      latest_index="${index}"');
    write_line(fid, '      latest="./${RESTART_PREFIX}.${index}.nc"');
    write_line(fid, '    fi');
    write_line(fid, '  done');
    write_line(fid, '  shopt -u nullglob');
    write_line(fid, '');
    write_line(fid, '  printf ''%s\n'' "${latest}"');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'restart_index() {');
    write_line(fid, '  local base index');
    write_line(fid, '  base="${1##*/}"');
    write_line(fid, '  index="${base#${RESTART_PREFIX}.}"');
    write_line(fid, '  index="${index%%.*}"');
    write_line(fid, '  printf ''%s\n'' "${index}"');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'target_reached() {');
    write_line(fid, '  local latest index year');
    write_line(fid, '  latest="$(latest_restart_file)"');
    write_line(fid, '');
    write_line(fid, '  [[ -n "${latest}" ]] || return 1');
    write_line(fid, '');
    write_line(fid, '  index="$(restart_index "${latest}")"');
    write_line(fid, '  year="${index:0:4}"');
    write_line(fid, '');
    write_line(fid, '  (( 10#${year} >= 10#${TARGET_YEAR} ))');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'backup_old_restart_files() {');
    write_line(fid, '  local latest latest_index rst base rest index suffix rank dest');
    write_line(fid, '');
    write_line(fid, '  latest="$(latest_restart_file)"');
    write_line(fid, '  if [[ -z "${latest}" ]]; then');
    write_line(fid, '    echo "No restart files found; skipping restart backup cleanup."');
    write_line(fid, '    return 0');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  latest_index="$(restart_index "${latest}")"');
    write_line(fid, '  mkdir -p "${BACKUP_DIR}"');
    write_line(fid, '');
    write_line(fid, '  shopt -s nullglob');
    write_line(fid, '  for rst in ./"${RESTART_PREFIX}".*.nc; do');
    write_line(fid, '    base="${rst##*/}"');
    write_line(fid, '    rest="${base#${RESTART_PREFIX}.}"');
    write_line(fid, '    index="${rest%%.*}"');
    write_line(fid, '    suffix="${rest#${index}.}"');
    write_line(fid, '');
    write_line(fid, '    [[ "${index}" =~ ^[0-9]{14}$ ]] || continue');
    write_line(fid, '    if [[ "${suffix}" != "nc" ]]; then');
    write_line(fid, '      rank="${suffix%.nc}"');
    write_line(fid, '      # rank suffix of any width, as in latest_restart_file');
    write_line(fid, '      [[ "${rank}" =~ ^[0-9]+$ ]] || continue');
    write_line(fid, '    fi');
    write_line(fid, '    [[ "${index}" != "${latest_index}" ]] || continue');
    write_line(fid, '');
    write_line(fid, '    dest="${BACKUP_DIR}/${base}"');
    write_line(fid, '    if [[ -e "${dest}" ]]; then');
    write_line(fid, '      echo "Backup already contains ${dest}; leaving ${rst} in place."');
    write_line(fid, '      continue');
    write_line(fid, '    fi');
    write_line(fid, '');
    write_line(fid, '    echo "Moving old restart ${rst} to ${BACKUP_DIR}/."');
    write_line(fid, '    mv "${rst}" "${BACKUP_DIR}/"');
    write_line(fid, '  done');
    write_line(fid, '  shopt -u nullglob');
    write_line(fid, '');
    write_line(fid, '  echo "Kept latest restart set in the run directory: ${latest}"');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'update_initial_file() {');
    write_line(fid, '  local input_file="$1"');
    write_line(fid, '  local restart_file="$2"');
    write_line(fid, '  local tmp_file="${input_file}.tmp.$$"');
    write_line(fid, '');
    write_line(fid, '  if [[ ! -f "${input_file}.initial" ]]; then');
    write_line(fid, '    cp -p "${input_file}" "${input_file}.initial"');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  if awk -v rst="${restart_file}" ''');
    write_line(fid, '    BEGIN {');
    write_line(fid, '      state = 0');
    write_line(fid, '      changed = 0');
    write_line(fid, '    }');
    write_line(fid, '    /^[[:space:]]*initial:[[:space:]]+NRREC[[:space:]]+filename[[:space:]]*$/ {');
    write_line(fid, '      state = 1');
    write_line(fid, '      print');
    write_line(fid, '      next');
    write_line(fid, '    }');
    write_line(fid, '    state == 1 {');
    write_line(fid, '      # EXACT_RESTART writes two records per restart event; a restart');
    write_line(fid, '      # must read the second one (cold starts keep NRREC=1 for init).');
    write_line(fid, '      print "          2"');
    write_line(fid, '      state = 2');
    write_line(fid, '      next');
    write_line(fid, '    }');
    write_line(fid, '    state == 2 {');
    write_line(fid, '      print "     " rst');
    write_line(fid, '      state = 0');
    write_line(fid, '      changed++');
    write_line(fid, '      next');
    write_line(fid, '    }');
    write_line(fid, '    {');
    write_line(fid, '      print');
    write_line(fid, '    }');
    write_line(fid, '    END {');
    write_line(fid, '      if (changed != 1) {');
    write_line(fid, '        exit 2');
    write_line(fid, '      }');
    write_line(fid, '    }');
    write_line(fid, '  '' "${input_file}" > "${tmp_file}"; then');
    write_line(fid, '    mv "${tmp_file}" "${input_file}"');
    write_line(fid, '  else');
    write_line(fid, '    rm -f "${tmp_file}"');
    write_line(fid, '    echo "ERROR: could not update the initial-condition filename in ${input_file}." >&2');
    write_line(fid, '    exit 2');
    write_line(fid, '  fi');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'partitioned_input_files_present() {');
    write_line(fid, '  local input_dir="Neptune_input"');
    write_line(fid, '  local expected_parts="${NTASKS}"');
    write_line(fid, '  local require_init="${1:-1}"');
    write_line(fid, '  local -a grid_parts frc_parts init_parts flux_parts');
    write_line(fid, '');
    write_line(fid, '  # A single rank reads the unpartitioned files directly.');
    write_line(fid, '  (( expected_parts == 1 )) && return 0');
    write_line(fid, '  # Rank-suffix width = digits of the rank COUNT, as in ROMS insert_nodes');
    write_line(fid, '  # (20 ranks -> .00-.19, 128 ranks -> .000-.127). partit writes one digit');
    write_line(fid, '  # fewer when the count is exactly 10, 100 or 1000.');
    write_line(fid, '  local width="${#expected_parts}"');
    write_line(fid, '  local rank_glob first_part');
    write_line(fid, '  rank_glob=$(printf ''[0-9]%.0s'' $(seq 1 "${width}"))');
    write_line(fid, '  first_part=$(printf "%0${width}d" 0)');
    write_line(fid, '');
    write_line(fid, '  shopt -s nullglob');
    write_line(fid, '  grid_parts=("${input_dir}"/neptune_grid.${rank_glob}.nc)');
    write_line(fid, '  frc_parts=("${input_dir}"/neptune_frc.${rank_glob}.nc)');
    write_line(fid, '  init_parts=("${input_dir}"/neptune_init.${rank_glob}.nc)');
    write_line(fid, '  flux_parts=("${input_dir}"/neptune_flux.${rank_glob}.nc)');
    write_line(fid, '  shopt -u nullglob');
    write_line(fid, '');
    write_line(fid, '  (( ${#grid_parts[@]} == expected_parts )) || return 1');
    write_line(fid, '  (( ${#frc_parts[@]} == expected_parts )) || return 1');
    write_line(fid, '  if [[ "${require_init}" == "1" && -f "${input_dir}/neptune_init.nc" ]]; then');
    write_line(fid, '    (( ${#init_parts[@]} == expected_parts )) || return 1');
    write_line(fid, '  fi');
    write_line(fid, '  if [[ -f "${input_dir}/neptune_flux.nc" ]]; then');
    write_line(fid, '    (( ${#flux_parts[@]} == expected_parts )) || return 1');
    write_line(fid, '  fi');
    write_line(fid, '  # Partitions made by a partit built without the fix in partit.F');
    write_line(fid, '  # lack all scalar variables (spherical, xl, ...); treat such sets');
    write_line(fid, '  # as absent so they are regenerated instead of silently reused.');
    write_line(fid, '  # (grep must read the whole stream: with pipefail, grep -q''s early');
    write_line(fid, '  # exit SIGPIPEs ncdump and turns a successful match into a failure)');
    write_line(fid, '  if command -v ncdump >/dev/null && [[ -f "${input_dir}/neptune_grid.${first_part}.nc" ]]; then');
    write_line(fid, '    if ! ncdump -h "${input_dir}/neptune_grid.${first_part}.nc" | grep spherical >/dev/null; then');
    write_line(fid, '      echo "Partitioned grid lacks scalar variables (partit built without the fix in partit.F); repartitioning."');
    write_line(fid, '      return 1');
    write_line(fid, '    fi');
    write_line(fid, '  fi');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'partition_input_files_if_needed() {');
    write_line(fid, '  local require_init="$1"');
    write_line(fid, '');
    write_line(fid, '  if partitioned_input_files_present "${require_init}"; then');
    write_line(fid, '    echo "Partitioned input files are already present; skipping partit."');
    write_line(fid, '    return 0');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  if [[ ! -f Neptune_input/do_partit.sh ]]; then');
    write_line(fid, '    echo "ERROR: Neptune_input/do_partit.sh is missing; cannot partition input files." >&2');
    write_line(fid, '    exit 2');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  echo "Partitioned input files not found; running Neptune_input/do_partit.sh."');
    write_line(fid, '  (cd Neptune_input && sh do_partit.sh)');
    write_line(fid, '');
    write_line(fid, '  if ! partitioned_input_files_present "${require_init}"; then');
    write_line(fid, '    echo "ERROR: partit finished, but expected partitioned input files are still missing." >&2');
    write_line(fid, '    exit 2');
    write_line(fid, '  fi');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'compile_roms_model() {');
    write_line(fid, '  echo "Compiling ROMS model."');
    write_line(fid, '  if [[ ! -d "${ROMS_ROOT}/src" ]]; then');
    write_line(fid, '    echo "ERROR: ROMS_ROOT=${ROMS_ROOT} does not contain src/." >&2');
    write_line(fid, '    exit 2');
    write_line(fid, '  fi');
    write_line(fid, '  # The experiment needs a development-branch tree, which has');
    write_line(fid, '  # src/energy_budget.F and no src/Makedefs.inc.');
    write_line(fid, '  if [[ -f "${ROMS_ROOT}/src/Makedefs.inc" || ! -f "${ROMS_ROOT}/src/energy_budget.F" ]]; then');
    write_line(fid, '    echo "ERROR: ROMS_ROOT=${ROMS_ROOT} is not a ucla-roms" >&2');
    write_line(fid, '    echo "       development-branch tree, which this experiment needs." >&2');
    write_line(fid, '    exit 2');
    write_line(fid, '  fi');
    write_line(fid, '  if [[ ! -f "${ROMS_ROOT}/NHMG/lib/libnhmg.a" ]]; then');
    write_line(fid, '    echo "ERROR: NHMG library missing; build it once with:" >&2');
    write_line(fid, '    echo "       cd ${ROMS_ROOT}/NHMG && make" >&2');
    write_line(fid, '    exit 2');
    write_line(fid, '  fi');
    write_line(fid, '  set +e');
    write_line(fid, '  make');
    write_line(fid, '  local make_status=$?');
    write_line(fid, '  set -e');
    write_line(fid, '');
    write_line(fid, '  if (( make_status == 0 )); then');
    write_line(fid, '    return 0');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  echo "Initial make failed with status ${make_status}." >&2');
    write_line(fid, '  if [[ ! -d Compile || ! -f Compile/Makefile ]]; then');
    write_line(fid, '    echo "Compile/Makefile is not available for a safe retry." >&2');
    write_line(fid, '    return "${make_status}"');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  echo "Prebuilding ROMS check tools serially, then retrying make once." >&2');
    write_line(fid, '  (cd Compile && rm -f strings.o strings.mod checkkwds.o checkkwds && make strings.o && make checkkwds)');
    write_line(fid, '  make');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, '# ROMS exits with status 0 even on BLOWUP (it shuts down cleanly after');
    write_line(fid, '# writing the emergency history record), so the exit code cannot be');
    write_line(fid, '# trusted alone: scan this job''s own SLURM log for the blow-up markers.');
    write_line(fid, 'blowup_detected() {');
    write_line(fid, '  local logfile');
    write_line(fid, '  for logfile in log.${SLURM_JOB_ID:-none}.*.log; do');
    write_line(fid, '    [[ -f "${logfile}" ]] || continue');
    write_line(fid, '    if grep -qE ''Abnormal termination: BLOWUP|setting may day flag'' "${logfile}"; then');
    write_line(fid, '      return 0');
    write_line(fid, '    fi');
    write_line(fid, '  done');
    write_line(fid, '  return 1');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'resubmit_requested=0');
    write_line(fid, '');
    write_line(fid, 'submit_next_job() {');
    write_line(fid, '  local reason="$1"');
    write_line(fid, '  local latest');
    write_line(fid, '  local -a dependency_arg=()');
    write_line(fid, '');
    write_line(fid, '  if [[ -f "${SCRIPT_DIR}/BLOWUP_DETECTED" ]]; then');
    write_line(fid, '    echo "BLOWUP_DETECTED is present; refusing to submit a continuation job (${reason})." >&2');
    write_line(fid, '    return 1');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  if (( resubmit_requested == 1 )); then');
    write_line(fid, '    echo "Continuation job already submitted; skipping duplicate request (${reason})."');
    write_line(fid, '    return 0');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  latest="$(latest_restart_file)"');
    write_line(fid, '  if [[ -z "${latest}" ]]; then');
    write_line(fid, '    echo "No restart file was found; cannot safely submit a continuation job." >&2');
    write_line(fid, '    return 1');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  if target_reached; then');
    write_line(fid, '    echo "Target year ${TARGET_YEAR} has been reached by ${latest}; not resubmitting."');
    write_line(fid, '    submit_ncjoin_job');
    write_line(fid, '    return 0');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  resubmit_requested=1');
    write_line(fid, '  if [[ -n "${SLURM_JOB_ID:-}" ]]; then');
    write_line(fid, '    dependency_arg=(--dependency="afterany:${SLURM_JOB_ID}")');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  echo "Submitting continuation job (${reason}) from latest restart ${latest}."');
    write_line(fid, '  sbatch "${dependency_arg[@]}" "${SCRIPT_PATH}"');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, '# Submit the ncjoin job once the simulation has reached the target year.');
    write_line(fid, '# The afterany dependency makes it start only after this job has left the');
    write_line(fid, '# queue, so do_ncjoin''s SIM_RUNNING=auto detection sees a finished run and');
    write_line(fid, '# joins every output set.');
    write_line(fid, 'submit_ncjoin_job() {');
    write_line(fid, '  local ncjoin_path="${SCRIPT_DIR}/${NCJOIN_SCRIPT}"');
    write_line(fid, '  local -a dependency_arg=()');
    write_line(fid, '');
    write_line(fid, '  if [[ ! -f "${ncjoin_path}" ]]; then');
    write_line(fid, '    echo "WARNING: ${ncjoin_path} not found; not submitting the ncjoin job." >&2');
    write_line(fid, '    return 0');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  if [[ -n "${SLURM_JOB_ID:-}" ]]; then');
    write_line(fid, '    dependency_arg=(--dependency="afterany:${SLURM_JOB_ID}")');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  echo "Simulation finished: submitting ${NCJOIN_SCRIPT} to join the output files."');
    write_line(fid, '  sbatch "${dependency_arg[@]}" "${ncjoin_path}"');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'handle_walltime_signal() {');
    write_line(fid, '  echo "Received Slurm pre-time-limit signal at $(date)."');
    write_line(fid, '  submit_next_job "pre-time-limit signal" || true');
    write_line(fid, '}');
    write_line(fid, '');
    write_line(fid, 'trap handle_walltime_signal USR1');
    write_line(fid, '');
    write_line(fid, 'echo "Target restart year: ${TARGET_YEAR}"');
    write_line(fid, 'latest_restart="$(latest_restart_file)"');
    write_line(fid, 'is_restart_run=0');
    write_line(fid, '');
    write_line(fid, 'if [[ -n "${latest_restart}" ]]; then');
    write_line(fid, '  is_restart_run=1');
    write_line(fid, '  latest_index="$(restart_index "${latest_restart}")"');
    write_line(fid, '  echo "Latest restart before run: ${latest_restart} (${latest_index})"');
    write_line(fid, '');
    write_line(fid, '  if target_reached; then');
    write_line(fid, '    echo "Target year ${TARGET_YEAR} is already reached; exiting."');
    write_line(fid, '    backup_old_restart_files');
    write_line(fid, '    submit_ncjoin_job');
    write_line(fid, '    exit 0');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  update_initial_file "${INPUT_FILE}" "${latest_restart}"');
    write_line(fid, '  echo "Updated ${INPUT_FILE} to restart from ${latest_restart}."');
    write_line(fid, '  backup_old_restart_files');
    write_line(fid, 'else');
    write_line(fid, '  echo "No restart file found; using the existing initial-condition entry in ${INPUT_FILE}."');
    write_line(fid, 'fi');
    write_line(fid, '');
    write_line(fid, '# MVAPICH2 flags:');
    write_line(fid, 'export MV2_USE_RDMA_CM=0');
    write_line(fid, 'export MV2_IBA_HCA=mlx5_2');
    write_line(fid, 'export MV2_DEFAULT_PORT=1');
    write_line(fid, '');
    write_line(fid, 'module purge');
    write_line(fid, 'module load slurm');
    write_line(fid, 'module load cpu/0.15.4  intel/19.1.1.217  mvapich2/2.3.4');
    write_line(fid, 'module load netcdf-c/4.7.4');
    write_line(fid, 'module load netcdf-fortran/4.5.3');
    write_line(fid, '');
    write_line(fid, '# Build environment for the development-branch ROMS.');
    write_line(fid, '# ROMS_ROOT is set here: with --export=ALL the job would otherwise');
    write_line(fid, '# inherit the login shell''s value, which may point to another tree.');
    fprintf(fid, 'export ROMS_ROOT=%s\n', roms_root);
    write_line(fid, 'export FCOMP=mpifort');
    write_line(fid, 'export NETCDF_ROOT="${NETCDF_FORTRANHOME}"');
    write_line(fid, 'export NETCDFC_ROOT="${NETCDF_CHOME}"');
    write_line(fid, '');
    write_line(fid, 'compile_roms_model');
    write_line(fid, '');
    write_line(fid, 'if (( is_restart_run == 1 )); then');
    write_line(fid, '  partition_input_files_if_needed 0');
    write_line(fid, 'else');
    write_line(fid, '  partition_input_files_if_needed 1');
    write_line(fid, 'fi');
    write_line(fid, '');
    write_line(fid, 'srun --mpi=pmi2 -n "${NTASKS}" roms "${INPUT_FILE}" &');
    write_line(fid, 'roms_pid=$!');
    write_line(fid, '');
    write_line(fid, 'while true; do');
    write_line(fid, '  set +e');
    write_line(fid, '  wait "${roms_pid}"');
    write_line(fid, '  roms_status=$?');
    write_line(fid, '  set -e');
    write_line(fid, '');
    write_line(fid, '  if kill -0 "${roms_pid}" 2>/dev/null; then');
    write_line(fid, '    continue');
    write_line(fid, '  fi');
    write_line(fid, '');
    write_line(fid, '  break');
    write_line(fid, 'done');
    write_line(fid, '');
    write_line(fid, 'if blowup_detected; then');
    write_line(fid, '  echo "The model log reports a blow-up; stopping the resubmission chain." >&2');
    write_line(fid, '  {');
    write_line(fid, '    echo "blow-up detected at $(date) in job ${SLURM_JOB_ID:-unknown}"');
    write_line(fid, '    echo "fix the run setup, then re-enable runs with: rm BLOWUP_DETECTED"');
    write_line(fid, '  } > BLOWUP_DETECTED');
    write_line(fid, '  exit 1');
    write_line(fid, 'elif (( roms_status == 0 )); then');
    write_line(fid, '  backup_old_restart_files');
    write_line(fid, '  submit_next_job "normal completion"');
    write_line(fid, 'else');
    write_line(fid, '  echo "ROMS exited with status ${roms_status}; not submitting a new job unless one was already queued by the signal handler." >&2');
    write_line(fid, '  exit "${roms_status}"');
    write_line(fid, 'fi');

    clear cleaner
    if isunix || ismac
        fileattrib(filename, '+x', 'a');
    end

    fprintf('Auto-restarting SLURM batch script written to "%s"\n', filename);
end

function write_line(fid, line)
    fprintf(fid, '%s\n', line);
end

function walltime_text = format_walltime(walltime)
    % Normalize the requested wall time: numeric hours -> 'H:00:00';
    % otherwise require a SLURM time string 'HH:MM:SS' or 'D-HH:MM:SS'.
    if isnumeric(walltime)
        validateattributes(walltime, {'numeric'}, ...
            {'scalar', 'positive', 'integer'}, mfilename, 'walltime');
        walltime_text = sprintf('%d:00:00', walltime);
    else
        walltime_text = char(walltime);
        if isempty(regexp(walltime_text, '^(\d+-)?\d{1,3}:\d{2}:\d{2}$', 'once'))
            error(['walltime must be numeric hours or a SLURM time string ', ...
                'like ''36:00:00'' or ''2-00:00:00'', got ''%s''.'], walltime_text);
        end
    end

    dash_split = strsplit(walltime_text, '-');
    hours_text = extractBefore([dash_split{end} ':'], ':');
    total_hours = str2double(hours_text);
    if numel(dash_split) == 2
        total_hours = total_hours + 24*str2double(dash_split{1});
    end
    if total_hours > 48
        warning('generate_run_script:walltime', ...
            ['Requested wall time %s exceeds the 48 h limit of the Expanse ', ...
            'compute partition; SLURM will reject the job there.'], walltime_text);
    end
end

function job_name = job_name_from_expname(expname)
    % Derive the SLURM job name from the trailing '_'-delimited segment of
    % expname. Falls back to 'roms_auto' when expname is empty or has no
    % usable trailing segment.
    job_name = 'roms_auto';

    expname = char(expname);
    if isempty(expname)
        return;
    end

    parts = strsplit(expname, '_');
    last = strtrim(parts{end});
    if ~isempty(last)
        job_name = last;
    end
end
