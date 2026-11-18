# Idealized ridge experiments for UCLA-ROMS

MATLAB scripts that generate complete experiment directories for idealized
mid-ocean-ridge simulations with UCLA-ROMS
([CESR-lab/ucla-roms](https://github.com/CESR-lab/ucla-roms), `development`
branch). An experiment directory holds:

- the ROMS grid, initial-condition and forcing NetCDF files
- the ROMS namelist (`neptune.in`) and option files (`*.opt`)
- the SLURM run and join scripts and the partition script, for the Expanse
  cluster
- a Makefile and the source overrides the experiment is compiled with

## Entry points

- `NewExp_Neptune_ridge.m`: generates one 3-D ridge experiment. Called with
  a struct, it overrides `expname`, `h_ridge`, `RF_F0_rot`, `alpha0` or
  `zob`.
- `run_ridge_ensemble.m`: generates a set of one-at-a-time perturbation
  experiments through the generator. Without arguments it only lists them.
- `NewExp_Neptune_ridge_2D.m`: generates a 2-D (along-ridge-uniform) run
  from an existing 3-D experiment, with zero wind stress and the vertical
  mixing prescribed from the 3-D run's time- and along-ridge mean. It reads
  `Akv` and `Akt` from the joined `neptune_avg.*.nc` files in the 3-D
  experiment directory, so the 3-D run must write the averages file; the
  3-D generator writes `ocean_vars.opt` with `wrt_file_avg` off.

The experiment-specific parameters (forcing, topography, thermal expansion
coefficient, resolution) are set near the top of `NewExp_Neptune_ridge.m`;
the vertical stretching and the initial temperature profile are set further
down, in the sections `generate and plot vertical grid` and `Initial
Stratification Setup`. The run-control parameters are near its end: `time_step` and `ndtfast`
(barotropic step = `time_step/ndtfast`), `runtime_request` (SLURM wall time
per job segment) and `restart_stop_year` (the year at which the automatic
resubmission stops).

`create_netcdf_file.m`, `write_netcdf_variable.m` and
`assert_netcdf_no_nan.m` create the NetCDF files, read variables back after
writing them, and check the finished files for NaNs.

## ROMS version and source overrides

The scripts are written for the `development` branch at commit `a0e2896`
(2026-11-09). The experiment Makefile copies `$(ROMS_ROOT)/src` into
`Compile/`, copies the `*.F`, `*.h` and `*.opt` files of the experiment
directory over it, and compiles there. The generators put these source
files into the experiment directory:

| File | Change relative to `src` | Copied |
|---|---|---|
| `rho_eos.F`, `alfabeta.F` | `SALINITY` guards in the linear equation of state, needed to compile without `SALINITY` | always |
| `grid.F` | reads the nudging fields `t_targ`, `nudgc` and `uv_nudgc` from the grid file; with `VMIX_PRESCRIBED_Huaiyu`, also `Akv_pres` and `Akt_pres` | with `restoring` or `use_uvnudg`; always for 2-D runs |
| `step3d_t_ISO.F` | nudges temperature toward `t_targ` in the boundary strips (`TNUDG_Huaiyu`) | with `restoring` or `use_uvnudg` |
| `step3d_uv2.F` | damps `u`, `v`, `ubar` and `vbar` in the same strips (`UVNUDG_Huaiyu`) | with `restoring` or `use_uvnudg` |
| `ana_vmix.h` | sets `Akv` and `Akt` to the prescribed fields (`ANA_VMIX` with `VMIX_PRESCRIBED_Huaiyu`) | 2-D runs |

Each `.F` override is the complete `src` file of that commit with the
change added; `ana_vmix.h` replaces the `src` file of the same name. For
another commit of ROMS, apply the same changes to that commit's
files: copying these files over a different `src` also replaces whatever
else differs in them.

## partit

`partit.F` is `Tools-Roms/partit.F` with one fix. The original passes
compression settings to `nf90_def_var` for every variable, netCDF rejects
them for scalar variables, and the scalars (`spherical`, `xl`, `el`, `hc`,
`theta_s`, ...) are then missing from the partitioned grid and initial
files, so the model stops in `get_grid`. Install the fixed file once per
clone:

    rm -f $ROMS_ROOT/Tools-Roms/partit.f $ROMS_ROOT/Tools-Roms/partit.o
    cp partit.F $ROMS_ROOT/Tools-Roms/partit.F
    cd $ROMS_ROOT/Tools-Roms && make partit

Keep the `.F` extension: the Tools Makefile preprocesses only `.F` sources.

## Cluster setup

Clone the `development` branch of ucla-roms on the cluster, build the NHMG
library (`cd NHMG && make`) and the tools (`cd Tools-Roms && make`, for
`partit` and `ncjoin`), and install `partit.F` as above. The generated
scripts call `partit` and `ncjoin` by name and submit with `--export=ALL`,
so `$ROMS_ROOT/Tools-Roms` must be on the `PATH` of the shell that submits
the jobs. The generated run script stops with a message if `ROMS_ROOT` has
no built NHMG library.

## Settings

`newexp_defaults.m` holds the settings shared by the generators:

- `preferred_output_root`: where experiment directories are written
- `fallback_output_root`: used when the preferred root cannot be created or
  written; `generated_experiments/` in this folder
- `roms_root`: path of the ucla-roms clone on the cluster, exported as
  `ROMS_ROOT` by the run script
- `slurm_account`, `slurm_email`: SLURM allocation and mail address; empty
  values leave those directives out of the generated scripts. A bare user
  name as `slurm_email` delivers only to the cluster-local mailbox.
- `roms_tools_dir`, `postprocess_dir`: optional folders added to the MATLAB
  path

The values in the file are the author's; change them before use.

## MATLAB requirements

- NetCDF support (`netcdf.*`, `ncwrite`, `ncread`)
- `smoothn.m` (in this folder), used when the ridge crest is smoothed
- the `topo` helper class (`topo.m`, in this folder). The generator uses
  `topo.zlevs3` and `topo.uvp_masks`, and with the along-ridge wind option
  also `topo.u2rho` and `topo.v2rho`.
- not in this folder: `cmocean`, and for the along-ridge wind option the
  `haxby` colormap. They must be on the MATLAB path, for example through
  `roms_tools_dir` or `postprocess_dir`.

## Typical workflow

1. Set the parameters near the top of `NewExp_Neptune_ridge.m` and run it
   in MATLAB.
2. Inspect the generated experiment directory (in the output root, or in
   `generated_experiments/` when the preferred root was unavailable).
3. Copy the directory to the cluster and submit `run_roms_auto` with
   `sbatch`.

The generator changes MATLAB's working directory to the experiment
directory.

## Generated batch scripts

`run_roms_auto` (written by `generate_run_script.m`):

- compiles the model and partitions the input files when needed
- restarts from the newest `neptune_rst.<timestamp>.<rank>.nc` set and
  writes its name into `neptune.in`
- resubmits itself until the target calendar year is reached
- stops the chain when the model log reports a blow-up (ROMS exits with
  status 0 even then, so the log is scanned) and writes a `BLOWUP_DETECTED`
  file that also blocks queued continuations; remove the file to run again
- moves older restart sets into `backup/` and keeps the newest set in the
  run directory
- submits `do_ncjoin` when the target year is reached

`do_ncjoin` (written by `generate_ncjoin_script.m`) joins the partitioned
output files. While the run is still active it leaves the newest set of
each output type alone.

## Backup_old_tree

`Backup_old_tree/` holds the generators for the 2024-08 snapshot of
UCLA-ROMS; see `README_backup.md` there.
