using JSON, REopt, JuMP, DelimitedFiles, DataFrames, CSV, Pickle

#necessary functions 
function convert_string_to_float(entry::Any)
    #println("Entry: $entry")
    #println("Type: $(typeof(entry))")
    
    if entry isa String31
        entry = String(entry)  # Convert String31 to String
        num = round(parse(Float64, entry), digits=4)
        return num
    elseif entry isa AbstractString
        num = round(parse(Float64, entry), digits=4)  # Handle other strings
        return num
    elseif entry isa Float64
        return entry  # Already a Float64
    else
        #println("Unhandled type: $(typeof(entry))")
        return NaN  # Return NaN for unsupported types
    end
end

#get the max size
function hub_height(size_class::AbstractString)
    if size_class == "Bergey Excel 15"
        return "24m"
    elseif size_class == "Northern Power Systems 100"
        return "37m"
    elseif size_class == "Vestas V-47"
        return "60m"
    elseif size_class == "GE 1.5 MW"
        return "80m"
    elseif size_class == "Bespoke 6 MW 170"
        return "115m"
    else
        println("Not correct size class insertion")
    end
end

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
#columns to select from csv file for parcel data
cols = [:MatchID, :naicsCode, :place_name, :latitude, :longitude, :state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, :critical_habs_2_int, :wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :wind_ground_area, :under_1_acre, :turbine_0, :"turbine_0: Maximum Wind Capacity [kW]", :"turbine_0: # turbines", :turbine_1, :"turbine_1: Maximum Wind Capacity [kW]", :"turbine_1: # turbines", :turbine_2, :"turbine_2: Maximum Wind Capacity [kW]", :"turbine_2: # turbines", :turbine_3, :"turbine_3: Maximum Wind Capacity [kW]", :"turbine_3: # turbines", :turbine_4, :"turbine_4: Maximum Wind Capacity [kW]", :"turbine_4: # turbines", :turbine_5, :"turbine_5: Maximum Wind Capacity [kW]", :"turbine_5: # turbines"]
# Define renaming rules
rename_rules = Dict(
    :"turbine_0: Maximum Wind Capacity [kW]" => :turbine_0_max_cap_kW,
    :"turbine_0: # turbines" => :turbine_0_num_turbines,
    :"turbine_1: Maximum Wind Capacity [kW]" => :turbine_1_max_cap_kW,
    :"turbine_1: # turbines" => :turbine_1_num_turbines,
    :"turbine_2: Maximum Wind Capacity [kW]" => :turbine_2_max_cap_kW,
    :"turbine_2: # turbines" => :turbine_2_num_turbines,
    :"turbine_3: Maximum Wind Capacity [kW]" => :turbine_3_max_cap_kW,
    :"turbine_3: # turbines" => :turbine_3_num_turbines,
    :"turbine_4: Maximum Wind Capacity [kW]" => :turbine_4_max_cap_kW,
    :"turbine_4: # turbines" => :turbine_4_num_turbines,
    :"turbine_5: Maximum Wind Capacity [kW]" => :turbine_5_max_cap_kW,
    :"turbine_5: # turbines" => :turbine_5_num_turbines  # Handles string version of the name
)
new_cols = [:MatchID, :naicsCode, :place_name, :latitude, 
:longitude, :state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, 
:critical_habs_2_int, :wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :wind_ground_area, :under_1_acre, 
:turbine_0, :turbine_0_max_cap_kW, :turbine_0_num_turbines, 
:turbine_1, :turbine_1_max_cap_kW, :turbine_1_num_turbines, 
:turbine_2, :turbine_2_max_cap_kW, :turbine_2_num_turbines, 
:turbine_3, :turbine_3_max_cap_kW, :turbine_3_num_turbines, 
:turbine_4, :turbine_4_max_cap_kW, :turbine_4_num_turbines,
:turbine_5, :turbine_5_max_cap_kW, :turbine_5_num_turbines,]

function read_csv_parcel_file(file_path::String)
    # Read the CSV into a DataFrame
    initial_df = CSV.read(file_path, DataFrame)
    
    # Rename the columns using the rename_rules
    rename!(initial_df, rename_rules)
    #println(initial_df)

    #get the selected cols from the df
    df = initial_df[:, new_cols]
    df = df[.!ismissing.(df.MatchID), :]
    df = df[.!ismissing.(df.wind_ground_area), :]
    return df 
end

function select_wind_turbines(row::DataFrame, max_load::Real)
    # Initialize variables to store results
    selected_class = nothing
    selected_turbines = 0
    total_capacity = 0.0
    capacity_factor = 0.35
        
    # Loop through turbine classes in the row
    for i in 0:5  # Iterates over turbine_0 to turbine_5
        # Get the turbine class, total capacity, and number of turbines
        turbine_class = row[1, Symbol("turbine_$i")]
        println("The turbine class is ", turbine_class)
        
        total_capacity_class = row[1, Symbol("turbine_$(i)_max_cap_kW")]
        num_turbines = row[1, Symbol("turbine_$(i)_num_turbines")]

        if !ismissing(total_capacity_class) && !ismissing(num_turbines)
            total_capacity_class = total_capacity_class[1]
            num_turbines = num_turbines[1]
            
            println("Total capacity class: ", total_capacity_class)
            println("Total turbine count: ", num_turbines)
            
            # Calculate the capacity per turbine
            if num_turbines > 0
                capacity_per_turbine = total_capacity_class / num_turbines
                println("Capacity per turbine is: ", capacity_per_turbine)
                
                # Calculate max capacity based on load
                max_size_based_on_load = max_load / (capacity_factor * 8760)
                println("The max size (kW) based on load is ", max_size_based_on_load)
                
                counter = 0 # for reducing turbine count

                # Reduce number of turbines one by one
                while counter < num_turbines
                    turbines_left = max(0, num_turbines - counter)
                    total_capacity = turbines_left * capacity_per_turbine
                    println("Turbines left: ", turbines_left, " | Total capacity: ", total_capacity)
                    
                    if total_capacity <= max_size_based_on_load && total_capacity > 0
                        # Return when condition is satisfied
                        selected_class = turbine_class
                        selected_turbines = turbines_left
                        println("Selected Turbine Class: ", selected_class)
                        println("Number of Turbines: ", selected_turbines)
                        println("Total Capacity (kW): ", total_capacity)
                        return Dict(
                            "Selected Turbine Class" => selected_class,
                            "Number of Turbines" => selected_turbines,
                            "Total Capacity (kW)" => total_capacity
                        )
                    else
                        counter += 1
                    end
                end
            else
                println("Skipping turbine class ", turbine_class, " due to zero turbines.")
            end
        else
            println("Data is missing for turbine ", turbine_class)
        end
    end
    
    # If no valid configuration found
    return "No valid configuration found for the given max load."
end

function select_wind_turbines2(row::DataFrame, max_load::Real, cap_factor::Real, inputs::Dict)
    # Initialize variables to store results
    selected_class = nothing
    selected_turbines = 0
    total_capacity = 0.0
    capacity_factor = 0.35
    if cap_factor !== nothing
        capacity_factor = cap_factor
    end
        
    # Loop through turbine classes in the row
    for i in 0:5  # Iterates over turbine_0 to turbine_5
        # Get the turbine class, total capacity, and number of turbines
        turbine_class = row[1, Symbol("turbine_$i")]
        #println("The turbine class is ", turbine_class)
        
        total_capacity_class = row[1, Symbol("turbine_$(i)_max_cap_kW")]
        num_turbines = row[1, Symbol("turbine_$(i)_num_turbines")]

        if !ismissing(total_capacity_class) && !ismissing(num_turbines) && i == 0
            total_capacity_class = total_capacity_class[1]
            num_turbines = num_turbines[1]
            
            #println("Total capacity class: ", total_capacity_class)
            #println("Total turbine count: ", num_turbines)
            
            # Calculate the capacity per turbine
            if num_turbines > 0
                capacity_per_turbine = total_capacity_class / num_turbines
                #println("Capacity per turbine is: ", capacity_per_turbine)
                
                # Calculate max capacity based on load
                max_size_based_on_load = max_load / (capacity_factor * 8760)
                #println("The max size (kW) based on load is ", max_size_based_on_load)
                
                counter = 0 # for reducing turbine count

                # Reduce number of turbines one by one
                while counter < num_turbines
                    turbines_left = max(0, num_turbines - counter)
                    total_capacity = turbines_left * capacity_per_turbine
                    #println("Turbines left: ", turbines_left, " | Total capacity: ", total_capacity)
                    
                    if total_capacity <= max_size_based_on_load && total_capacity > 0
                        # Return when condition is satisfied
                        selected_class = turbine_class
                        selected_turbines = turbines_left
                        #println("Selected Turbine Class: ", selected_class)
                        #println("Number of Turbines: ", selected_turbines)
                        #println("Total Capacity (kW): ", total_capacity)
                        return Dict(
                            "Selected Turbine Class" => selected_class,
                            "Number of Turbines" => selected_turbines,
                            "Total Capacity (kW)" => total_capacity
                        )
                    else
                        counter += 1
                    end
                end
            else
                println("Skipping turbine class ", turbine_class, " due to zero turbines.")
            end
        elseif !ismissing(total_capacity_class) && !ismissing(num_turbines) && i > 0
            total_capacity_class = total_capacity_class[1]
            num_turbines = num_turbines[1]
            
            #println("Total capacity class: ", total_capacity_class)
            #println("Total turbine count: ", num_turbines)
            
            # Calculate the capacity per turbine
            if num_turbines > 0
                capacity_per_turbine = total_capacity_class / num_turbines
                #println("Capacity per turbine is: ", capacity_per_turbine)

                input_data_s = copy(inputs)

                input_dat_s["Wind"]["size_class"] = turbine_class
                
                # Calculate max capacity based on load
                max_size_based_on_load = max_load / (capacity_factor * 8760)
                #println("The max size (kW) based on load is ", max_size_based_on_load)
                
                counter = 0 # for reducing turbine count

                # Reduce number of turbines one by one
                while counter < num_turbines
                    turbines_left = max(0, num_turbines - counter)
                    total_capacity = turbines_left * capacity_per_turbine
                    #println("Turbines left: ", turbines_left, " | Total capacity: ", total_capacity)
                    
                    if total_capacity <= max_size_based_on_load && total_capacity > 0
                        # Return when condition is satisfied
                        selected_class = turbine_class
                        selected_turbines = turbines_left
                        #println("Selected Turbine Class: ", selected_class)
                        #println("Number of Turbines: ", selected_turbines)
                        #println("Total Capacity (kW): ", total_capacity)
                        return Dict(
                            "Selected Turbine Class" => selected_class,
                            "Number of Turbines" => selected_turbines,
                            "Total Capacity (kW)" => total_capacity
                        )
                    else
                        counter += 1
                    end
                end
            else
                println("Skipping turbine class ", turbine_class, " due to zero turbines.")
            end
        else
            println("Data is missing for turbine ", turbine_class)
        end
    end
    
    # If no valid configuration found
    return "No valid configuration found for the given max load."
end
function select_num_turbines(row::DataFrame, max_load::Real, cap_factor::Real, step::Integer)
    # Initialize variables to store results
    selected_class = nothing
    selected_turbines = 0
    total_capacity = 0.0
    capacity_factor = 0.35
    if cap_factor !== nothing
        capacity_factor = cap_factor
    end
        
    # Get the turbine class, total capacity, and number of turbines
    turbine_class = row[1, Symbol("turbine_$step")]
    #println("The turbine class is ", turbine_class)
    
    total_capacity_class = row[1, Symbol("turbine_$(step)_max_cap_kW")]
    num_turbines = row[1, Symbol("turbine_$(step)_num_turbines")]

    total_capacity_class = total_capacity_class[1]
    num_turbines = num_turbines[1]

    if num_turbines > 0

        capacity_per_turbine = total_capacity_class / num_turbines

        # Calculate max capacity based on load
        max_size_based_on_load = max_load / (capacity_factor * 8760)
        println("The max size (kW) based on load is ", max_size_based_on_load)
        
        counter = 0 # for reducing turbine count

        # Reduce number of turbines one by one
        while counter < num_turbines
            turbines_left = max(0, num_turbines - counter)
            total_capacity = turbines_left * capacity_per_turbine
            println("Turbines left: ", turbines_left, " | Total capacity: ", total_capacity)
            
            if total_capacity <= max_size_based_on_load && total_capacity > 0
                # Return when condition is satisfied
                selected_class = turbine_class
                selected_turbines = turbines_left
                println("Selected Turbine Class: ", selected_class)
                println("Number of Turbines: ", selected_turbines)
                println("Total Capacity (kW): ", total_capacity)
                return Dict(
                    "Selected Turbine Class" => selected_class,
                    "Number of Turbines" => selected_turbines,
                    "Total Capacity (kW)" => total_capacity
                )
            else
                counter += 1
            end
        end
    else
        println("No valid configuration found for the given max load.")
        return nothing
    end
end

function get_wind_parameters(;
    match_id::Any,
    folder_path_pickl::String #the wind resource parameters pickle file
    )
    # Locate the file that matches the match_id
    matching_files = String[]
    files = readdir(folder_path_pickl)
    for f in files
        if occursin(string(match_id), f)
            file_path = joinpath(folder_path_pickl, f)
            push!(matching_files, file_path)
            println("The file path for site $match_id is ", file_path)
        end
    end

    # Handle case where file was not found
    if isempty(matching_files)
        println("Site $match_id does not have wind resource data")
        return nothing
    end

    #store results
    wind_resources = Dict()

    #process each individual file
    for file in matching_files
        println("Processing $file")
        # Load the pickle file
        dict_pckl = Pickle.npyload(file)
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

            parts = split(file, "-") # Split by '-'
            height_text = split(parts[end], ".")[1] # Take the last part, i.e. the turbine height + m

            wind_resources[height_text] = Dict()
            wind_resources[height_text]["wind_m_per_sec_1"] = wind_m_per_sec
            wind_resources[height_text]["wind_direction_degrees_1"] = wind_direction_degrees
            wind_resources[height_text]["temperature_celsius_1"] = temperature_celsius
            wind_resources[height_text]["atmos_pressure_1"] = atmos_pressure
            wind_resources[height_text]["heights"] = 1

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

            parts = split(file, "-") # Split by '-'
            height_text = split(parts[end], ".")[1] # Take the last part, i.e. the turbine height + m

            wind_resources[height_text] = Dict()
            wind_resources[height_text]["wind_m_per_sec_1"] = wind_m_per_sec_height_1
            wind_resources[height_text]["wind_direction_degrees_1"] = wind_direction_degrees_height_1
            wind_resources[height_text]["temperature_celsius_1"] = temperature_celsius_height_1
            wind_resources[height_text]["atmos_pressure_1"] = atmos_pressure_height_1
            wind_resources[height_text]["wind_m_per_sec_2"] = wind_m_per_sec_height_2
            wind_resources[height_text]["wind_direction_degrees_2"] = wind_direction_degrees_height_2
            wind_resources[height_text]["temperature_celsius_2"] = temperature_celsius_height_2
            wind_resources[height_text]["atmos_pressure_2"] = atmos_pressure_height_2
            wind_resources[height_text]["heights"] = 2

        else
            println("Incomplete information for site $match_id")
            return nothing
        end
    end
    return wind_resources
end
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

# Folder paths
#get load data from IEDO Teams
electric_load_folder_path = "C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/Load Profiles/"
load_traits_text = "Load Facility Traits Set"
load_data_text = "Load Facility Set"
ng_load_folder_path = "C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/NG Consumption/"

#wind resource pickle files path
wind_resource_path = "C:/Users/dbernal/Documents/GitHub/wind_resource_data/"

#Set-up inputs file for PV runs 
data_file = "wind_runs.json"
input_data = JSON.parsefile("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/Input Resources/$data_file")

#parcel file path in IEDO Teams 
parcel_file = "C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/Input Resources/wind_sitelist_3x7_no_sites_with_exclusions.csv"

#get data from CSV file for parcel data 
data = read_csv_parcel_file(parcel_file)