
1. login to kestrel and open a terminal. Navigate to your scratch folder and make a folder for the batch scripts.
    ```
    cd ../../scratch/dbernal
    mkdir batch_files
    ```
2. copy the batch script for 100 sites iedo_batch_test_100sites.sh and to the folder: /scratch/dbernal/batch_files
3. submit the batch job for 100 sites
    ```
    sbatch /scratch/dbernal/batch_files/iedo_batch_test_100sites.sh
    ```
4. make note of the job id (write it down somewhere). you can now log out of kestrel.
5. you'll get an email when the job has started, if it failed, or when it finished. 
    - if it finishes successfully (exit code 0) - take note of the simulation time that is stated in the email that says the job finished.
    - if it fails - to see the slurm log that will have error information, 
        1. login to kestrel and open a terminal
            ```
            cd ../../scratch/dbernal
            ```
        2. open the error log which will have the filename R-<job_name>.<job_id>.out. <job_name> is the name of the job specified in the .sh file you submitted. <job_id> is the number you wrote down in step 4. for example:
            ```
            code R-iedo_tech_potential_test_100.<job_id>.out
            ```