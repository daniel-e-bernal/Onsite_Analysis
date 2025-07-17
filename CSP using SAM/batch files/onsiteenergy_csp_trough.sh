#!/bin/bash 
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --account=iedo00onsite
#SBATCH --ntasks=50 # CPUs requested for job 
#SBATCH --mem-per-cpu=2000 # Request 2GB per core.
#SBATCH --job-name=testing_run
#SBATCH --partition=shared
#SBATCH --array=1-10

module load julia
julia run_ssc_trough.jl