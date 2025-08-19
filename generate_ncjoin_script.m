function generate_ncjoin_script(account_name,user_name,filename)
% generate_ncjoin_script([filename])
% Writes the requested ncjoin SLURM script in the current directory.
% Default: 'ncjoin_mpi50.sh'

    if nargin < 3 || isempty(filename)
        filename = 'do_ncjoin';
    end

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open "%s" for writing.', filename);
    end

    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#SBATCH --job-name="ncj_mpi50"\n\n');

    fprintf(fid, '#       %%j=job_number and %%N gives nodelist output="wec_real.%%j.%%N.out"\n');
    fprintf(fid, '#SBATCH --output="ncjoin_mpi50.%%j.%%N.log"\n\n');

    fprintf(fid, '#       If using less than 128 cores then the partition is ''shared''\n');
    fprintf(fid, '#       or ''large-shared'' if extra memory is needed (only use if ''shared'' failed due to memory error).\n');
    fprintf(fid, '#       Note, the job charge costs 4x more for ''large-shared'' than ''shared'' or ''compute'', use allocation wisely!\n');
    fprintf(fid, '#SBATCH --partition=shared\n');
    fprintf(fid, '#SBATCH --nodes=1\n\n');

    fprintf(fid, '#       Request number of cores: (Expanse has 128 cores per node)\n');
    fprintf(fid, '#       If using ''compute'', you will be charged 128 x number of nodes requested! Regardless of ntasks=per-node,\n');
    fprintf(fid, '#       which is very important to understand to avoid being overcharged and wasting resources.\n');
    fprintf(fid, '#       For ncjoin_mpi, not recommend to go higher than 50, as diminishing returns.\n');
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
    fprintf(fid, '#SBATCH --account=%s\n', account_name);
    fprintf(fid, '#SBATCH --export=ALL\n');
    fprintf(fid, '#       Time duration requested for run:\n');
    fprintf(fid, '#SBATCH -t 24:00:00\n');
    fprintf(fid, '#SBATCH --mail-user=%s\n', user_name);
    fprintf(fid, '#SBATCH --mail-type=ALL\n\n');

    fprintf(fid, '#-----------------------------------------------------------------\n\n');

    fprintf(fid, '# Flags needed for mvapich2:\n');
    fprintf(fid, 'export MV2_USE_RDMA_CM=0\n');
    fprintf(fid, 'export MV2_IBA_HCA=mlx5_2\n');
    fprintf(fid, 'export MV2 DEFAULT PORT=1\n\n');

    fprintf(fid, 'module purge\n');
    fprintf(fid, 'module load slurm\n');
    fprintf(fid, 'module load cpu/0.15.4  intel/19.1.1.217  mvapich2/2.3.4\n');
    fprintf(fid, 'module load netcdf-c/4.7.4\n');
    fprintf(fid, 'module load netcdf-fortran/4.5.3\n');
    fprintf(fid, 'ncjoin -d neptune_avg.*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_dia.*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2000*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2001*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2002*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2003*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2004*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2005*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2006*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2007*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2008*.*.nc &\n');
    fprintf(fid, 'ncjoin -d neptune_his.2009*.*.nc &\n');

    fclose(fid);

    % Make executable on Unix/macOS
    if isunix || ismac
        fileattrib(filename, '+x', 'a');
    end

    fprintf('Wrote "%s" in the current directory.\n', filename);
end
