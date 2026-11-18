# Idealized ridge experiments for UCLA-ROMS

MATLAB scripts that generate complete experiment directories for idealized
mid-ocean-ridge simulations with UCLA-ROMS
([CESR-lab/ucla-roms](https://github.com/CESR-lab/ucla-roms), `development`
branch). An experiment directory holds:

- the ROMS grid, initial-condition and forcing NetCDF files
- the ROMS namelist (`neptune.in`) and option files (`*.opt`)
- the SLURM run and join scripts and the partition script, for the Expanse cluster
- a Makefile and the source overrides the experiment is compiled with

## Entry points

- `NewExp_Neptune_ridge.m`: generates one 3-D ridge experiment. 
- `run_ridge_ensemble.m`: generates a set of one-at-a-time perturbation
  experiments through the generator. 
- `NewExp_Neptune_ridge_2D.m`: generates a 2-D (along-ridge-uniform) run
  from an existing 3-D experiment, with zero wind stress and the vertical
  mixing prescribed from the 3-D run's time- and along-ridge mean. 

The experiment-specific parameters (forcing, topography, thermal expansion
coefficient, resolution) are set near the top of `NewExp_Neptune_ridge.m`;
the vertical stretching and the initial temperature profile are set further
down, in the sections `generate and plot vertical grid` and `Initial
Stratification Setup`. The run-control parameters are near its end: `time_step` and `ndtfast`
(barotropic step = `time_step/ndtfast`), `runtime_request` (SLURM wall time
per job segment) and `restart_stop_year` (the year at which the automatic
resubmission stops).

## Typical workflow

1. Set the parameters near the top of `NewExp_Neptune_ridge.m` and run it in MATLAB.
2. Inspect the generated experiment directory (in the output root, or in
   `generated_experiments/` when the preferred root was unavailable).
3. Copy the directory to the cluster and submit `run_roms_auto` with `sbatch`.

## Cluster Configuration

The cluster-related files are tailored for the [Expanse cluster](https://www.sdsc.edu/systems/expanse/user_guide.html) at UC San Diego (UCSD).

Obtaining computational hours on Expanse is relatively straightforward by applying through the [ACCESS program](https://allocations.access-ci.org/).

##  Resources

- Additional UCLA-ROMS tools are available at:  
  👉 [https://github.com/nmolem/ucla-tools](https://github.com/nmolem/ucla-tools)
