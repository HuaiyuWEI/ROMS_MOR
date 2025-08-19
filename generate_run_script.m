function generate_run_script(num_nodes, account_name,user_name, filename)
    % Default filename if not provided
    if nargin < 4
        filename = 'run_roms';
    end

    fid = fopen(filename, 'w');
    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#SBATCH --job-name="roms"\n\n');
    fprintf(fid, '# %%j=job_number and %%N gives nodelist output="wec_real.%%j.%%N.out"\n');
    fprintf(fid, '#SBATCH --output="log.%%j.%%N.log"\n\n');
    fprintf(fid, '# Can only use a max of 32 nodes on ''compute'' partition:\n');
    fprintf(fid, '#SBATCH --partition=compute\n');
    fprintf(fid, '#SBATCH --nodes=%d\n\n', num_nodes);
    fprintf(fid, '# Request number of cores (Expanse has 128 cores per node)\n');
    fprintf(fid, '#SBATCH --ntasks-per-node=128\n');
    fprintf(fid, '#SBATCH --account=%s\n', account_name);
    fprintf(fid, '#SBATCH --export=ALL\n\n');
    % fprintf(fid, '# Memory settings (default: 1GB per core)\n');
    % fprintf(fid, '#SBATCH --mem-per-cpu=2G\n\n');
    fprintf(fid, '# Max run time is 48 hours:\n');
    fprintf(fid, '#SBATCH -t 48:00:00\n\n');
    fprintf(fid, '# Email notifications:\n');
    fprintf(fid, '#SBATCH --mail-user=%s\n', user_name);
    fprintf(fid, '#SBATCH --mail-type=ALL\n\n');
    fprintf(fid, '#-----------------------------------------------------------------\n\n');
    fprintf(fid, '# MVAPICH2 flags:\n');
    fprintf(fid, 'export MV2_USE_RDMA_CM=0\n');
    fprintf(fid, 'export MV2_IBA_HCA=mlx5_2\n');
    fprintf(fid, 'export MV2_DEFAULT_PORT=1\n\n');
    fprintf(fid, 'module purge\n');
    fprintf(fid, 'module load slurm\n');
    fprintf(fid, 'module load cpu/0.15.4  intel/19.1.1.217  mvapich2/2.3.4\n');
    fprintf(fid, 'module load netcdf-c/4.7.4\n');
    fprintf(fid, 'module load netcdf-fortran/4.5.3\n\n');
    fprintf(fid, 'cd Neptune_input/\n');
    fprintf(fid, 'sh do_partit.sh\n');
    fprintf(fid, 'cd ..\n');
    fprintf(fid, 'make\n');
    fprintf(fid, 'srun --mpi=pmi2 -n %d roms neptune.in\n', num_nodes * 128);

    fclose(fid);

    fprintf('SLURM batch script written to "%s"\n', filename);
end
