using DataFrames
using CSV 

# Define path to xlsx file
file_storage_location = "./results/validation_run_result.csv"

# Write DataFrame to a new or existing CSV file
if isfile(file_storage_location)
    # Append new data to the CSV file
    open(file_storage_location, "a") do io
        CSV.write(io, df, header=false)  # Append without headers
    end
else
    # Write DataFrame to a new CSV file
    CSV.write(file_storage_location, df)
end