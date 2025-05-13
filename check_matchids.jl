using CSV, DataFrames, FilePaths, DelimitedFiles


folder_path_e = "C:/Users/dbernal/Downloads/Facility Load Profiles/" #the electric loads folder path
traits_file = "Load Facility Traits Set "
set_file = "Load Profiles Set "

#the files are numbered, therefore, we'll set a counter 
counter = 0
max_counter = 33
a = collect(counter:max_counter)

loads_w_matchids = []
#loads_w_matchids = Set()

for i in a
    #we get the traits path first 
    file_path_e = joinpath(folder_path_e, "$traits_file$(string(i)).csv")
    file_e = CSV.read(file_path_e, DataFrame)
    filter!(row -> !ismissing(row.MatchID), file_e)
    #println(file_e)
    
    matchids_list = file_e.MatchID
    unique_ids = unique(matchids_list)
    uniqueids_len = length(unique_ids)
    #println("Length of unique IDs for iteration $i is: ", uniqueids_len)
    runs = collect(1:uniqueids_len)

    for i in runs
        b = unique_ids[i]
        println(b)
        b = String(b)
        println(b)
        if !(b in loads_w_matchids)
            push!(loads_w_matchids, b)
        end
    end
    # Extract and add unique MatchIDs to the set
    #union!(loads_w_matchids, unique(file_e.MatchID))
end

loads_w_matchids = collect(loads_w_matchids)
println(length(loads_w_matchids))
df = DataFrame(MatchID=loads_w_matchids)

# Define path to csv file
file_storage_location = "C:/Users/dbernal/Downloads/filter_w_this.csv"

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