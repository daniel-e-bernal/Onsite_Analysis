function copy_unique_files(source_dir::String, target_dir::String)
    # Get a list of files in the source and target directories
    source_files = readdir(source_dir)  # List of files in source directory
    target_files = readdir(target_dir)  # List of files in target directory

    for file in source_files
        source_path = joinpath(source_dir, file)
        target_path = joinpath(target_dir, file)
        
        # Only copy if it's a file (not a folder) and if it does not exist in the target directory
        if isfile(source_path) && !(file in target_files)
            try
                cp(source_path, target_path)
                println("Copied: $file to $target_dir")
            catch e
                println("Failed to copy $file: ", e)
            end
        end
    end
end

# Example usage
source_dir = "C:/Users/dbernal/OneDrive - NREL/General - IEDO Onsite Energy/Data/Batch Results/PV_results/"
target_dir = "C:/GitRepos/Onsite_Energy_temporary/results/PV/results/"
copy_unique_files(source_dir, target_dir)
