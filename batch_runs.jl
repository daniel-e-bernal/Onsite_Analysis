#!/bin/bash                                                                                                             
#SBATCH --job-name=iedo_tech_potential_chunk01
#SBATCH --output=R-%x.%j.out
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --time=1:00:00
#SBATCH --account=iedo00onsite
#SBATCH --mail-user dbernal@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL
module load julia
activate /scratch/dbernal/"something.jl"
julia --project=. "something.jl"
export TMPDIR=/scratch/dbernal/sc_tmp/
srun -ntasks=1 /scratch/dbernal/path/to/julia /scratch/dbernal/script_julia.jl 1000 0