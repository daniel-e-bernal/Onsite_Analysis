
1. login to kestrel and open a terminal. Navigate to your scratch folder and open any files you may want to edit during your run (you wont be able to open the once the interactive node has started)
    ```
    cd ../../scratch/dbernal
    code Wind_PV_Tech_Potential/tests_v3.jl
    ```
2. request an interactive node session for 1 hour with 10 GB of memory
    ```
    salloc --time=60 --account=iedo00onsite --ntasks=1 --partition=shared --mem=10GB
    ```
3. once the session has started, load the necessary modules and active your environment
    ```
    module purge
    export TMPDIR=/scratch/dbernal/sc_tmp/
    module load julia/1.8.5-generic-linux
    julia --project=/scratch/dbernal/Wind_PV_Tech_Potential/
    include("tests_v3.jl")
    ```
4. try to run your code for 10 sites
    ```
    tech_potential_reopt(num_sites=10,index_start=1,batch_number=1)
    ```
5. if no bugs arise, exit the interactive session
    ```
    exit
    ```