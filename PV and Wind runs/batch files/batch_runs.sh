#!/bin/bash                                                                                                             
#SBATCH --job-name=iedo_tech_potential_chunk01
#SBATCH --output=R-%x.%j.out
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --time=1:00:00
#SBATCH --account=iedo00onsite
#SBATCH --mail-user daniel.bernal@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL
module purge
export TMPDIR=/scratch/dbernal/sc_tmp/
module load julia/1.8.5-generic_linux
julia --project=/scratch/dbernal/Wind_PV_Tech_Potential/ 
include("tests_v3.jl")
tech_potential_reopt(num_sites=1000, index_start=1, batch_number=1)
