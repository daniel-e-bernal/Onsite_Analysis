using CSV, DataFrames, FilePaths, Pickle, DelimitedFiles

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
    max_counter = 34

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
            #println("Got stuck trying to find a match file", counter)
            counter += 1
            continue
        else #so if a match was found 
            simu_id = find_match_id_row.simulationID[1]
            simu_id = lpad(simu_id, 6, '0')
            #estimated_or_real = find_match_id_row.energy_estimated_or_random[1]
            #append!(electric_estimated_or_random, "$estimated_or_real")
            #println("Got the simu_id = ", simu_id)
        end
        ng_load = [] 
        #find the MatchID in the natural gas file
        find_match_id_row_ng = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_ng)
        if isempty(find_match_id_row_ng)
            #println("Did not find a MatchID for the natural gas consumption.")
            append!(ng_load, 0)
            #append!(ng_estimated_or_random, "NaN")
        else #so if a match was found 
            annual_ng_mmbtu = find_match_id_row_ng.Estimated_Annual_Natural_Gas_MMBtu[1]
            #estimated_or_real_ng = find_match_id_row_ng.Natural_Gas_Estimated_E_or_Random_R[1]
            append!(ng_load, annual_ng_mmbtu)
            #append!(ng_estimated_or_random, "$estimated_or_real_ng")
            #println("The annual MMBtu consumption is  = ", annual_ng_mmbtu)
        end

        #get the file path for the set file that contains the hourly loads 
        file_path_s = joinpath(folder_path_e, "$set_file$(string(counter)).csv")
        file_s = CSV.read(file_path_s, DataFrame)
        electric_hourly_load = []
        append!(electric_hourly_load, file_s[!, simu_id])
        return electric_hourly_load,  ng_load
    end
end

function get_load_data(;
    match_id::Any,
    folder_path_e::String, #the electric loads folder path
    folder_path_ng::String, #the ng loads folder path
    traits_file::String = "Load Facility Traits Set ",
    set_file::String = "Load Profiles Set ",
    base_load_path::String)
    
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
    file_path_ng = joinpath(folder_path_ng, "Manufacturing Parcel Data - Natural Gas Estimates_updated_Mar2025.csv")
    file_ng = CSV.read(file_path_ng, DataFrame)
    file_path_b = joinpath(base_load_path, "baseload_sizes_all.csv")
    file_b = CSV.read(file_path_b, DataFrame)

    while true
        #check counter 
        if counter == max_counter
            return 0 #exit loop
        end

        #we get the traits path first 
        file_path_e = joinpath(folder_path_e, "$traits_file$(string(counter)).csv")
        file_e = CSV.read(file_path_e, DataFrame)
        #make global variable for simu_id
        simu_id = 0
        #find the MatchID in the electric loads file 
        find_match_id_row = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_e)
        if isempty(find_match_id_row)
            #println("Got stuck trying to find a match file", counter)
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
        #find the simulation ID in the baseload file under the column name :Profile and the baseload kW is under the column "LDC 0.2 Val"
        baseload_kw = []
        simu_id2 = parse(Int, simu_id) #convert the simu_id to an Int
        find_match_id_row_b = filter(row -> !ismissing(row.Profile) && row.Profile == simu_id2, file_b)
        if isempty(find_match_id_row_b)
            println("Did not find a MatchID for the baseload consumption.")
            append!(baseload_kw, 0)
        else #so if a match was found 
            baseload_size = find_match_id_row_b[!, "LDC 0.2 Val"][1]
            append!(baseload_kw, baseload_size)
            println("The baseload kw capacity is  = ", baseload_size)
        end
        #get the file path for the set file that contains the hourly loads 
        file_path_s = joinpath(folder_path_e, "$set_file$(string(counter)).csv")
        file_s = CSV.read(file_path_s, DataFrame)
        electric_hourly_load = []
        append!(electric_hourly_load, file_s[!, simu_id])
        return electric_hourly_load,  ng_load, baseload_kw
    end
end

function select_num_turbines(
                            row::DataFrame, 
                            max_load::Real, 
                            cap_factor::Real, 
                            step::Integer, 
                            option::String, 
                            baseload_kw::Real,
                            method::String = "normal"
                            )
    # Initialize variables to store results
    selected_class = nothing
    selected_turbines = 0
    total_capacity = 0.0
    capacity_factor = 0.35
    if cap_factor !== nothing
        capacity_factor = cap_factor
    else
        @warn "Using arbitrary 35% capacity factor instead of ", cap_factor, "% - calculated by SAM"
    end
        
    # Get the turbine class, total capacity, and number of turbines
    turbine_class = row[1, Symbol("turbine_$step")]
    #println("The turbine class is ", turbine_class)
    
    total_capacity_class = row[1, Symbol("turbine_$(step)_max_cap_kW")]
    num_turbines = row[1, Symbol("turbine_$(step)_num_turbines")]

    total_capacity_class = total_capacity_class[1]
    num_turbines = num_turbines[1]

    if option == "A" && method == "normal" #size up to baseload (kW)
        if num_turbines > 0 
            # Calculate how many turbines are needed to suffice the baseload size (kW) for option A
            capacity_per_turbine = total_capacity_class / num_turbines
            counter = 0 # for reducing turbine count
            #reduce number of turbines one by one 
            while counter < num_turbines 
                turbines_left = max(0, num_turbines - counter)
                total_capacity = turbines_left * capacity_per_turbine
                if total_capacity <= baseload_kw && total_capacity > 0
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
            println("No turbines available for option A.")
            return nothing
        end
    elseif option == "A" && method == "+1" #size up to baseload (kW) + 1 turbine
        if num_turbines == 1
            #if the total number of turbines you can fit is 1, then just output that single turbine for that size class
            return Dict(
                        "Selected Turbine Class" => turbine_class,
                        "Number of Turbines" => num_turbines,
                        "Total Capacity (kW)" => total_capacity_class
                    )
        elseif num_turbines > 1
            # Calculate how many turbines are needed to suffice the baseload size (kW) for option A
            capacity_per_turbine = total_capacity_class / num_turbines
            counter = 0 # for reducing turbine count
            #reduce number of turbines one by one 
            while counter <= num_turbines 
                turbines_left = max(0, num_turbines - counter)
                total_capacity = turbines_left * capacity_per_turbine
                if total_capacity <= baseload_kw && total_capacity >= 0
                    # Return when condition is satisfied
                    selected_class = turbine_class
                    selected_turbines = turbines_left
                    selected_total_capacity = total_capacity 
                    if counter > 0
                        #add one more turbine
                        selected_turbines = turbines_left + 1
                        #adjust total capacity for adding 1 more turbine of that size class
                        selected_total_capacity = total_capacity + capacity_per_turbine
                    end
                    println("Selected Turbine Class: ", selected_class)
                    println("Number of Turbines: ", selected_turbines)
                    println("Total Capacity (kW): ", selected_total_capacity)
                    return Dict(
                        "Selected Turbine Class" => selected_class,
                        "Number of Turbines" => selected_turbines,
                        "Total Capacity (kW)" => selected_total_capacity
                    )
                else
                    counter += 1
                end
            end
        else 
            println("No turbines available for option A.")
            return nothing
        end
    elseif option == "B" && method == "normal" #max turbine capacity based on total load (kWh)
        if num_turbines > 0

            capacity_per_turbine = total_capacity_class / num_turbines

            # Calculate max capacity based on load
            max_size_based_on_load = max_load / (capacity_factor * 8760)
            #println("The max size (kW) based on load is ", max_size_based_on_load)
            
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
    elseif option == "B" && method == "+1" #max turbine capacity based on total load (kWh) + 1 turbine
        if num_turbines == 1
            #if the total number of turbines you can fit is 1, then just output that single turbine for that size class
            return Dict(
                "Selected Turbine Class" => turbine_class,
                "Number of Turbines" => num_turbines,
                "Total Capacity (kW)" => total_capacity_class
            )
        elseif num_turbines > 1
             #if the total number of turbines is more than one, then go along with the drop down loop
            capacity_per_turbine = total_capacity_class / num_turbines
    
            # Calculate max capacity based on load
            max_size_based_on_load = max_load / (capacity_factor * 8760)
            println("The max size (kW) based on load is ", max_size_based_on_load)
            
            counter = 0 # for reducing turbine count
    
            # Reduce number of turbines one by one
            while counter <= num_turbines
                turbines_left = max(0, num_turbines - counter)
                total_capacity = turbines_left * capacity_per_turbine
                println("Turbines left: ", turbines_left, " | Total capacity: ", total_capacity)
                
                if total_capacity <= max_size_based_on_load && total_capacity >= 0
                    # Return when condition is satisfied
                    selected_class = turbine_class
                    selected_turbines = turbines_left 
                    selected_total_capacity = total_capacity 
                    if counter > 0
                        #add one more turbine
                        selected_turbines = turbines_left + 1
                        #adjust total capacity for adding 1 more turbine of that size class
                        selected_total_capacity = total_capacity + capacity_per_turbine
                    end
                    println("Selected Turbine Class: ", selected_class)
                    println("Number of Turbines: ", selected_turbines)
                    println("Total Capacity (kW): ", selected_total_capacity)
                    return Dict(
                        "Selected Turbine Class" => selected_class,
                        "Number of Turbines" => selected_turbines,
                        "Total Capacity (kW)" => selected_total_capacity
                    )
                else
                    counter += 1
                end
            end
        else
            println("No valid configuration found for the given max load.")
            return nothing
        end
    elseif option == "C" #max turbine capacity based on space
        if num_turbines > 0
            #start calculating based on total land 
            return Dict(
                "Selected Turbine Class" => turbine_class,
                "Number of Turbines" => num_turbines,
                "Total Capacity (kW)" => total_capacity
            )
        else 
            println("No turbines available for option C.")
            return nothing
        end
    end
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
