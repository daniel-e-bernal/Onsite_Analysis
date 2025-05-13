#!/bin/bash                                                                                                             
#SBATCH --job-name=iedo_tech_potential_test_100
#SBATCH --output=R-%x.%j.out
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --mem=50GB
#SBATCH --time=1:00:00
#SBATCH --account=iedo00onsite
#SBATCH --mail-user daniel.bernal@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL
module purge
export TMPDIR=/scratch/dbernal/sc_tmp/
module load julia/1.8.5-generic-linux
julia --project=/scratch/dbernal/Wind_PV_Tech_Potential/ -e '
    include("./Wind_PV_Tech_Potential/tests_v3.jl")
    tech_potential_reopt(num_sites=100,index_start=1,batch_number=1)'
