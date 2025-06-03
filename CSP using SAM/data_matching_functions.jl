"""
This file below contains the functions necessary for the data matching process.
The functions are:
    A) electric load profile match 
    B) ng consumption match
    C) resource data match **may be updated**
"""

using DataFrames, FilePaths, CSV, XLSX, DelimitedFiles, Pickle

# Get load data 
function get_load_data(;
    match_id::Any,
    folder_path_e::String, #the electric loads folder path
    folder_path_ng::String, #the ng loads folder path
    traits_file::String = "Load Facility Traits Set ",
    set_file::String = "Load Profiles Set ")
    
    #prepare to intake any type of String
    if match_id isa String31
        match_id = String(match_id)  # Convert String31 to String
    elseif match_id isa AbstractString
        match_id = String(match_id)  # Handle other strings
    elseif match_id isa String15
        match_id = String(match_id)
    else
        match_id
    end
    #the files are numbered, therefore, we'll set a counter 
    counter = 0
    max_counter = 35

    while true
        #check counter 
        if counter == max_counter
            return 0 #exit loop
        end

        #we get the traits path first 
        file_path_e = joinpath(folder_path_e, "$traits_file$(string(counter)).csv")
        file_e = CSV.read(file_path_e, DataFrame)
        file_path_ng = joinpath(folder_path_ng, "Manufacturing Parcel Data - Natural Gas Estimates.csv")
        file_ng = CSV.read(file_path_ng, DataFrame)
        #electric_estimated_or_random = []
        #find the MatchID in the electric loads file 
        find_match_id_row = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_e)
        if isempty(find_match_id_row)
            println("Got stuck trying to find a match file", counter)
            counter += 1
            continue
        else #so if a match was found 
            simu_id = find_match_id_row.simulationID[1]
            simu_id = lpad(simu_id, 6, '0')
            #estimated_or_real = find_match_id_row.energy_estimated_or_random[1]
            #append!(electric_estimated_or_random, "$estimated_or_real")
            println("Got the simu_id = ", simu_id)
        end
        ng_load = [] 
        #find the MatchID in the natural gas file
        find_match_id_row_ng = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_ng)
        if isempty(find_match_id_row_ng)
            println("Did not find a MatchID for the natural gas consumption.")
            append!(ng_load, 0)
            #append!(ng_estimated_or_random, "NaN")
        else #so if a match was found 
            annual_ng_mmbtu = find_match_id_row_ng.Estimated_Annual_Natural_Gas_MMBtu[1]
            #estimated_or_real_ng = find_match_id_row_ng.Natural_Gas_Estimated_E_or_Random_R[1]
            append!(ng_load, annual_ng_mmbtu)
            #append!(ng_estimated_or_random, "$estimated_or_real_ng")
            println("The annual MMBtu consumption is  = ", annual_ng_mmbtu)
        end

        #get the file path for the set file that contains the hourly loads 
        file_path_s = joinpath(folder_path_e, "$set_file$(string(counter)).csv")
        file_s = CSV.read(file_path_s, DataFrame)
        electric_hourly_load = []
        append!(electric_hourly_load, file_s[!, simu_id])
        return electric_hourly_load,  ng_load
    end
end

# Get load data 
function get_load_data2(;
    match_id::Any,
    folder_path_e::String, #the electric loads folder path
    folder_path_ng::String, #the ng loads folder path
    traits_file::String = "Load Facility Traits Set ",
    set_file::String = "Load Profiles Set ")
    
    #prepare to intake any type of String
    if match_id isa String31
        match_id = String(match_id)  # Convert String31 to String
    elseif match_id isa AbstractString
        match_id = String(match_id)  # Handle other strings
    elseif match_id isa String15
        match_id = String(match_id)
    else
        match_id
    end
    #the files are numbered, therefore, we'll set a counter 
    counter = 0
    max_counter = 35

    while true
        #check counter 
        if counter == max_counter
            return 0 #exit loop
        end

        #we get the traits path first 
        file_path_e = joinpath(folder_path_e, "$traits_file$(string(counter)).csv")
        file_e = CSV.read(file_path_e, DataFrame)
        file_path_ng = joinpath(folder_path_ng, "Manufacturing Parcel Data - Natural Gas Estimates.csv")
        file_ng = CSV.read(file_path_ng, DataFrame)
        #electric_estimated_or_random = []
        #find the MatchID in the electric loads file 
        find_match_id_row = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_e)
        if isempty(find_match_id_row)
            println("Got stuck trying to find a match file", counter)
            counter += 1
            continue
        else #so if a match was found 
            simu_id = find_match_id_row.simulationID[1]
            simu_id = lpad(simu_id, 6, '0')
            #estimated_or_real = find_match_id_row.energy_estimated_or_random[1]
            #append!(electric_estimated_or_random, "$estimated_or_real")
            println("Got the simu_id = ", simu_id)
        end
        ng_load = [] 
        #find the MatchID in the natural gas file
        find_match_id_row_ng = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_ng)
        if isempty(find_match_id_row_ng)
            println("Did not find a MatchID for the natural gas consumption.")
            append!(ng_load, 0)
            #append!(ng_estimated_or_random, "NaN")
        else #so if a match was found 
            annual_ng_mmbtu = find_match_id_row_ng.Estimated_Annual_Natural_Gas_MMBtu[1]
            #estimated_or_real_ng = find_match_id_row_ng.Natural_Gas_Estimated_E_or_Random_R[1]
            append!(ng_load, annual_ng_mmbtu)
            #append!(ng_estimated_or_random, "$estimated_or_real_ng")
            println("The annual MMBtu consumption is  = ", annual_ng_mmbtu)
        end

        #get the file path for the set file that contains the hourly loads 
        file_path_s = joinpath(folder_path_e, "$set_file$(string(counter)).csv")
        file_s = CSV.read(file_path_s, DataFrame)
        electric_hourly_load = []
        append!(electric_hourly_load, file_s[!, simu_id])
        return electric_hourly_load,  ng_load
    end
end

#select production factor series based on match_id
function select_prod_factor(; 
    match_id::Any,
    roof_prod_factor_folder::String = pv_roof_prod_factors,
    ground_prod_fixed_factor_folder::String = pv_ground_fixed_prod_factors,
    ground_prod_axis_factor_folder::String = pv_ground_axis_prod_factors)

    #prepare to intake any type of String
    if match_id isa String31
        match_id = String(match_id)  # Convert String31 to String
    elseif match_id isa AbstractString
        match_id = String(match_id)  # Handle other strings
    elseif match_id isa String15
        match_id = String(match_id)
    else
        match_id
    end

    get_roof_file = joinpath(roof_prod_factor_folder, "$match_id.csv")
    get_ground_fixed_file = joinpath(ground_prod_fixed_factor_folder, "$match_id.csv")
    get_ground_axis_file = joinpath(ground_prod_axis_factor_folder, "$match_id.csv")

    # Check if either file exists
    roof_exists = isfile(get_roof_file)
    ground_fixed_exists = isfile(get_ground_fixed_file)
    ground_axis_exists = isfile(get_ground_axis_file)

    if !roof_exists
        error("No production factor file found for match_id: $match_id in roof_prod_factor_folder")
    end
    if !ground_fixed_exists && !ground_axis_exists
        error("No production factor file found for match_id: $match_id in ground_prod_fixed_factor_folder OR ground_prod_axis_factor_folder")
    end

    #read csv file and convert it to array, each file only has 1 column no header
    try
        data_ground_fixed = readdlm(get_ground_fixed_file, ',', Float64)  # Read only the first column
        data_ground_fixed = vec(data_ground_fixed)
        data_ground_axis = 0
        if ground_axis_exists == true
            data_ground_axis = readdlm(get_ground_axis_file, ',', Float64)  # Read only the first column
            data_ground_axis = vec(data_ground_axis)
        end
        #data_ground = collect(data_ground)  # Convert DataFrame column to an array
        data_roof = readdlm(get_roof_file, ',', Float64)
        data_roof = vec(data_roof)
        #data_roof = CSV.read(get_roof_file, DataFrame)[:, 1]  # Read only the first column
        #data_roof = collect(data_roof)  # Convert DataFrame column to an array
        return data_roof, data_ground_fixed, data_ground_axis
    catch e
        error("Error reading file $get_ground_fixed_file, $get_ground_axis_file and $get_roof_file: $e")
    end
end

"""
The functions below are primarily used for wind runs.
"""

function get_wind_parameters2(;
    match_id::Any,
    folder_path_pickl::String, #the wind resource parameters pickle file
    turbine_step::AbstractString #accepts a string that suggests the hub height ["24m", "37m", "60m", "80m", "115m"]
    )
    # Locate the file that matches the match_id
    file_path = nothing
    expected_file_name = match_id * "-2012-" * turbine_step * ".pkl"
    files = readdir(folder_path_pickl)
    if expected_file_name in files
        file_path = joinpath(folder_path_pickl, expected_file_name)
    end

    # Handle case where file was not found
    if file_path === nothing
        println("Site $match_id does not have wind resource data")
        return nothing
    end

    # Load the pickle file
    dict_pckl = Pickle.npyload(file_path)
    resources_count = length(dict_pckl["fields"])
    time_series_resource_length = length(dict_pckl["data"])  # should be 8760

    # Check if the file has 4 or 8 fields
    if resources_count == 4
        # Extract data for 4 fields
        wind_m_per_sec = []
        wind_direction_degrees = []
        temperature_celsius = []
        atmos_pressure = []
        println("Single hub height")
        for i in 1:time_series_resource_length
            append!(wind_m_per_sec, Float64(dict_pckl["data"][i][3]))
            append!(wind_direction_degrees, Float64(dict_pckl["data"][i][4]))
            append!(temperature_celsius, Float64(dict_pckl["data"][i][1]))
            append!(atmos_pressure, Float64(dict_pckl["data"][i][2]))  # Fixed to use 4th element
        end

        return Dict(
            "wind_m_per_sec_1" => wind_m_per_sec,
            "wind_direction_degrees_1" => wind_direction_degrees,
            "temperature_celsius_1" => temperature_celsius,
            "atmos_pressure_1" => atmos_pressure,
            "heights" => 1
        )

    elseif resources_count == 8
        # Extract data for 8 fields (two heights)
        height_1 = dict_pckl["heights"][1]
        height_2 = dict_pckl["heights"][5]  # Assuming the 5th index is relevant for height 2
        println("Getting resources for two hub heights")
        wind_m_per_sec_height_1 = []
        wind_direction_degrees_height_1 = []
        temperature_celsius_height_1 = []
        atmos_pressure_height_1 = []

        wind_m_per_sec_height_2 = []
        wind_direction_degrees_height_2 = []
        temperature_celsius_height_2 = []
        atmos_pressure_height_2 = []

        for i in 1:time_series_resource_length
            append!(wind_m_per_sec_height_1, Float64(dict_pckl["data"][i][3]))
            append!(wind_direction_degrees_height_1, Float64(dict_pckl["data"][i][4]))
            append!(temperature_celsius_height_1, Float64(dict_pckl["data"][i][1]))
            append!(atmos_pressure_height_1, Float64(dict_pckl["data"][i][2]))  # Fixed to use 4th element

            append!(wind_m_per_sec_height_2, Float64(dict_pckl["data"][i][7]))
            append!(wind_direction_degrees_height_2, Float64(dict_pckl["data"][i][8]))
            append!(temperature_celsius_height_2, Float64(dict_pckl["data"][i][5]))
            append!(atmos_pressure_height_2, Float64(dict_pckl["data"][i][6]))
        end

        return Dict(
            "wind_m_per_sec_1" => wind_m_per_sec_height_1,
            "wind_direction_degrees_1" => wind_direction_degrees_height_1,
            "temperature_celsius_1" => temperature_celsius_height_1,
            "atmos_pressure_1" => atmos_pressure_height_1,
            "wind_m_per_sec_2" => wind_m_per_sec_height_2,
            "wind_direction_degrees_2" => wind_direction_degrees_height_2,
            "temperature_celsius_2" => temperature_celsius_height_2,
            "atmos_pressure_2" => atmos_pressure_height_2,
            "heights" => 2
        )

    else
        println("Incomplete information for site $match_id")
        return nothing
    end
end


#=
I need to create a function that reads the match id in the LC_facility_parcels_NREL_05_18_25_TEST_CSP.csv file.
    Then, using the MatchID, go into the path folder C:\Users\dbernal\Downloads\ITO Load Profiles\ to extract the matching load profile.
    Once it finds the matching load profile, it should create an individual csv file with that load profile naming the file with the MatchID.
=#
function create_load_profile_csv( 
    match_id::Any,
    folder_path_e::String = "C:/Users/dbernal/Downloads/ITO Load Profiles/", #the electric loads folder path
    traits_file::String = "Load Facility Traits Set ",
    set_file::String = "Load Profiles Set ",
    output_folder::String = "./tests/Load Profiles/",
    override_existing_file::Bool = true) #output folder for the csv file

    if !override_existing_file
        println("Not overriding existing file")
        return
    end
    println("Overriding existing file")

    #prepare to intake any type of String
    if match_id isa String31
        match_id = String(match_id)  # Convert String31 to String
    elseif match_id isa AbstractString
        match_id = String(match_id)  # Handle other strings
    elseif match_id isa String15
        match_id = String(match_id)
    else
        match_id
    end
    
    #the files are numbered, therefore, we'll set a counter     
    counter = 0
    max_counter = 34

    #create empty DataFrame to store the electric load profile
    df = DataFrame()

    while true
        #check counter 
        if counter == max_counter
            return 0 #exit loop
        end

        #we get the traits path first 
        file_path_e = joinpath(folder_path_e, "$traits_file$(string(counter)).csv")
        file_e = CSV.read(file_path_e, DataFrame)

        #empty simu_id 
        simu_id = nothing

        #find the MatchID in the electric loads file 
        find_match_id_row = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_e)
        if isempty(find_match_id_row)
            println("Got stuck trying to find a match file", counter)
            counter += 1
            continue
        else #so if a match was found 
            simu_id = find_match_id_row.simulationID[1]
            simu_id = lpad(simu_id, 6, '0')
            #estimated_or_real = find_match_id_row.energy_estimated_or_random[1]
            #append!(electric_estimated_or_random, "$estimated_or_real")
            println("Got the simu_id = ", simu_id)
        end
        println("Onto the next step")
        #get the file path for the set file that contains the hourly loads 
        file_path_s = joinpath(folder_path_e, "$set_file$(string(counter)).csv")
        file_s = CSV.read(file_path_s, DataFrame)
        #electric_hourly_load = []
        #append!(electric_hourly_load, file_s[!, simu_id])
        df[!, :Load] = file_s[!, simu_id]
        save_to_csv(df, output_folder, match_id)
        println("The load profile for $match_id has been created in the folder $output_folder")
        break #exit loop
    end
end
    
#save to CSV file function
function save_to_csv(data::DataFrame, file_path::String, name::String)
    #create file path
    file = joinpath(file_path, "$(name)_electric.csv")

    # Write DataFrame to a new or existing CSV file
    if isfile(file)
        # Append new data to the CSV file
        open(file, "a") do io
            CSV.write(io, df, header=true)  # Append without headers
        end
    else
        # Write DataFrame to a new CSV file
        CSV.write(file, data)
    end
end

