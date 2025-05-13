# From https://researchcomputing.princeton.edu/support/knowledge-base/julia#distributed
# https://github.com/Arpeggeo/julia-distributed-computing

using Distributed, SlurmClusterManager

# Env instantiation
@everywhere begin
    using Pkg
    Pkg.activate(@__DIR__)
end

# Function declarations, module imports
@everywhere begin
    
    using JSON, REopt, JuMP, DelimitedFiles, DataFrames, CSV, Pickle

    ENV["NREL_DEVELOPER_API_KEY"]="gAXbkyLjfTFEFfiO3YhkxxJ6rkufRaSktk40ho4x"

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
    #columns to select from csv file for parcel data
    new_cols = [:MatchID, :naicsCode, :place_name, :latitude, 
    :longitude, :parcel_latitude, :parcel_longitude, :state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, 
    :critical_habs_2_int, :wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :wind_ground_area, :under_1_acre, 
    :turbine_0, :turbine_0_max_cap_kW, :turbine_0_num_turbines, 
    :turbine_1, :turbine_1_max_cap_kW, :turbine_1_num_turbines, 
    :turbine_2, :turbine_2_max_cap_kW, :turbine_2_num_turbines, 
    :turbine_3, :turbine_3_max_cap_kW, :turbine_3_num_turbines, 
    :turbine_4, :turbine_4_max_cap_kW, :turbine_4_num_turbines,
    :turbine_5, :turbine_5_max_cap_kW, :turbine_5_num_turbines,]

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

    function select_wind_turbines(row::DataFrame, max_load::Real, cap_factor::Real)
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

    function get_wind_parameters(;
        match_id::Any,
        folder_path_pickl::String #the wind resource parameters pickle file
        )
        # Locate the file that matches the match_id
        file_path = nothing
        files = readdir(folder_path_pickl)
        for f in files
            if occursin(string(match_id), f)
                file_path = joinpath(folder_path_pickl, f)
                break
            end
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
    
    # Add common site information
    function run_site(i::Int)
            
        #store results
        analysis_runs = DataFrame()
        emissions = DataFrame() 
        #for i in sites_iter
        input_data_site = copy(input_data)

        #get the standard inputs 
        #conversion for latitude
        latitude = data[!, :parcel_latitude][i]
        #latitude = String(latitude)
        latitude = convert_string_to_float(latitude)
        #converstion for longitude
        longitude = data[!, :parcel_longitude][i]
        #longitude = String(longitude)
        longitude = convert_string_to_float(longitude)
        #conversion for land_acres
        land_acres = data[!, :wind_ground_area][i]
        if ismissing(land_acres) || isnan(land_acres)
            land_acres = 0
        else 
            land_acres = round(land_acres / 4046.86, digits=4) #conversion from m2 to acres
        end
        
        #get MatchID in data DataFrame to start getting other data from loads 
        match_id = data[!, :MatchID][i]
        println("The site's MatchID is: ", match_id)

        #create the file name for this specific site + Wind (the other version is PV)
        file_name = string(match_id, "_Wind")
        site_load_info = get_load_data(match_id=match_id, folder_path_e=electric_load_folder_path, folder_path_ng=ng_load_folder_path)
        load_vector = site_load_info[1]
        
        #println(typeof(load_vector))
        ng_annual_mmbtu = site_load_info[2]
        
        input_data_site["ElectricLoad"]["loads_kw"] = load_vector
        input_data_site["ElectricLoad"]["annual_kwh"] = sum(load_vector)
        annual_load = sum(load_vector)
        
        #srmer_co2e_c <- for emissions reduction
        hourly_fuel = ng_annual_mmbtu[1] / 8760
        hourly_fuel = fill(hourly_fuel, 8760)
        input_data_site["DomesticHotWaterLoad"]["fuel_loads_mmbtu_per_hour"] = hourly_fuel

        #get wind parameters to create production factor 
        wind_parameters_dict = get_wind_parameters(match_id=match_id, folder_path_pickl=wind_resource_path)
        wind_speed = wind_parameters_dict["wind_m_per_sec_1"]
        #println("The length of wind speed is ", length(wind_speed), " and the type of object is ", typeof(wind_speed))
        #println("The first element for wind speed is ", wind_speed[1], " and the type of object for the element is ",typeof(wind_speed[1]))
        site_temp = wind_parameters_dict["temperature_celsius_1"]
        atmos_press = wind_parameters_dict["atmos_pressure_1"]
        wind_direction = wind_parameters_dict["wind_direction_degrees_1"]
        height_quantity = wind_parameters_dict["heights"]
        wind_speed_2 = []
        site_temp_2 = []
        atmos_press_2 = []
        wind_direction_2 = []
        if height_quantity == 1
            println("just one height")
        else 
            println("multiple heights")
            wind_speed_2 = wind_parameters_dict["wind_m_per_sec_2"]
            site_temp_2 = wind_parameters_dict["temperature_celsius_2"]
            atmos_press_2 = wind_parameters_dict["atmos_pressure_2"]
            wind_direction_2 = wind_parameters_dict["wind_direction_degrees_2"]
        end

        #input into the site specific data
        input_data_site["Site"]["latitude"] = latitude
        input_data_site["Site"]["longitude"] = longitude
        input_data_site["Site"]["land_acres"] = land_acres
        
        wind_speed = convert(Vector{Float64}, wind_speed)
        wind_speed_2 = convert(Vector{Float64}, wind_speed_2)

        wind_direction = convert(Vector{Float64}, wind_direction)
        wind_direction_2 = convert(Vector{Float64}, wind_direction_2)

        atmos_press = convert(Vector{Float64}, atmos_press)
        atmos_press_2 = convert(Vector{Float64}, atmos_press_2)

        site_temp = convert(Vector{Float64}, site_temp)
        site_temp_2 = convert(Vector{Float64}, site_temp_2)

        #input arrays to get production factor for the largest possible turbine
        input_data_site["Wind"]["wind_meters_per_sec"] = height_quantity == 1 ? wind_speed : [wind_speed, wind_speed_2]
        input_data_site["Wind"]["wind_direction_degrees"] = height_quantity == 1 ? wind_direction : [wind_direction, wind_direction_2]
        input_data_site["Wind"]["pressure_atmospheres"] = height_quantity == 1 ? atmos_press : [atmos_press, atmos_press_2]
        input_data_site["Wind"]["temperature_celsius"] = height_quantity == 1 ? site_temp : [site_temp, site_temp_2]
        #input_data_site["Wind"]["wind_meters_per_sec"] = wind_speed 
        #input_data_site["Wind"]["wind_direction_degrees"] = wind_direction 
        #input_data_site["Wind"]["pressure_atmospheres"] = atmos_press 
        #input_data_site["Wind"]["temperature_celsius"] = site_temp 
    
        #wind.production_factor_series
        wind_class = data[!, :turbine_0][i]
        println(wind_class)
        if wind_class == "Nothern Power Systems 100"
            wind_class = "Northern Power Systems 100"
            println("Corrected wind class name is ", wind_class)
        end
        sleep(10)
        max_cap_num_of_turbines = data[!, :turbine_0_num_turbines][i] #this value should be number of turbines for the largest size class possible 
        max_kw_based_on_space = data[!, :turbine_0_max_cap_kW][i] #this value should be size class * num of turbines of that size class that can fit in the plot 
        input_data_site["Wind"]["size_class"] = String(wind_class)
        input_data_site["Wind"]["max_kw"] = max_kw_based_on_space
        #input_data_site["Wind"]["min_kw"] = max_kw_based_on_space * 0.99
        println("The max size for the wind system will be ", input_data_site["Wind"]["max_kw"])
        """ Below is attaining the REopt inputs related to aer_gen_co2e_c emissions to calculate BAU emissions."""
        input_data_site["ElectricUtility"]["cambium_metric_col"] = "aer_gen_co2e_c"
        s1 = Scenario(input_data_site)
        println("Scenario passed for $match_id")
        sleep(3)
        inputs1 = REoptInputs(s1)
        println("REoptInputs passed ofr $match_id")

        #getting max size based on annual load (using capacity factor)
        wind_prod_data = inputs1.production_factor["Wind", :].data
        println("maximum ", maximum(wind_prod_data))
        println("minimum ", minimum(wind_prod_data))
        #wind_prod_data = [x[1] for x in wind_prod_data]
        #println("The type of data in the production factor series is ", typeof(wind_prod_data), " and the length is ", length(wind_prod_data))
        sleep(3)
        wind_capacity_factor = sum(wind_prod_data) / 8760 #actual capacity factor
        println("The site $match_id capacity factor for largest possible turbine is: ", wind_capacity_factor)
        wind_max_size_based_on_load = input_data_site["ElectricLoad"]["annual_kwh"] / (wind_capacity_factor * 8760)
        println("The site's $match_id max size (kW) based on load is: ", wind_max_size_based_on_load)

        #data
        site_row = data[i, :]
        site_row = DataFrame(site_row)

        #initialize 
        new_wind_class = "Kept original"
        new_max_wind_capacity = 0
        wind_size = max_kw_based_on_space
        num_turbines = max_cap_num_of_turbines
        new_wind_data = select_wind_turbines(site_row, annual_load, wind_capacity_factor)

        if wind_max_size_based_on_load >= max_kw_based_on_space 
            println("===========================================================================================")
            println("Moving forward with pre-determined size class and quantity of turbines of that class")
            println("===========================================================================================")
            sleep(5)
        else
            println("===========================================================================================")
            println("selecting new turbine class")
            println("===========================================================================================")
            return nothing
            """
            new_wind_class = new_wind_data["Selected Turbine Class"]
            println("The new wind turbine class for $match_id is ", new_wind_class)
            if new_wind_class == "Nothern Power Systems 100"
                new_wind_class = "Northern Power Systems 100"                 
            end
            input_data_site["Wind"]["wind_meters_per_sec"] = []
            input_data_site["Wind"]["wind_direction_degrees"] = [] 
            input_data_site["Wind"]["pressure_atmospheres"] = [] 
            input_data_site["Wind"]["temperature_celsius"] = []
            new_max_wind_capacity = new_wind_data["Total Capacity (kW)"]
            wind_size = new_max_wind_capacity
            input_data_site["Wind"]["size_class"] = new_wind_class
            input_data_site["Wind"]["max_kw"] = new_max_wind_capacity
            #input_data_site["Wind"]["min_kw"] = new_max_wind_capacity * 0.99
            num_turbines = new_wind_data["Number of Turbines"]"""
        end
        
        """
        #getting max size based on annual load (using capacity factor)
        wind_prod_data = inputs1.production_factor
        wind_capacity_factor = sum(wind_prod_data) / 8760 #actual capacity factor
        wind_max_size_based_on_load = input_data_site["ElectricLoad"]["annual_kwh"] / (wind_capacity_factor * 8760)
        #println("Max size for PV on ground (kW) based on annual load for # i is: ", PV_ground_max_size_based_on_load)
        #now getting max size based on space
        wind_max_size_based_on_space = land_acres * (1/wind_capacity_factor) #inputs1.max_sizes["ground_mount"] #max size based on space 
        #println("Max size for PV on roof (kW) based on roof space for # i is: ", PV_roof_max_size_based_n_space)
        inputs1.max_sizes["ground_mount"] = PV_ground_max_size_based_on_space
        #create global variables for PV sizes on roof and ground 
        wind_size = 0 #kW
        """
       
        """ Below is attaining the REopt inputs related to srmer_co2e_c emissions to calculate BAU emissions."""

        input_data_site["ElectricUtility"]["cambium_metric_col"] = "srmer_co2e_c"
        
        s2 = Scenario(input_data_site)
        inputs2 = REoptInputs(s2)

        #will get Wind related and necessary variables 
        wind_prod_factor_series = inputs2.production_factor["Wind", :].data #this gets the production factor series for 
        
        wind_prod_kwh_series = wind_prod_factor_series * wind_size
        wind_load_minus_prod_kwh_series = load_vector - wind_prod_kwh_series
        wind_export_kwh_series = [] #8760 intervals
        grid_supplied_kwh = [] #8760 intervals
        wind_serving_load_total_series = [] #8760 intervals 
        #place values in export series array 
        for i in eachindex(wind_load_minus_prod_kwh_series)
            if wind_load_minus_prod_kwh_series[i] <= 0
                push!(grid_supplied_kwh, 0)
                push!(wind_serving_load_total_series, load_vector[i])
                excess = wind_prod_kwh_series[i] - load_vector[i]
                push!(wind_export_kwh_series, excess)
            else
                grid_s = load_vector[i] - wind_prod_kwh_series[i]
                push!(grid_supplied_kwh, grid_s)
                push!(wind_serving_load_total_series, wind_prod_kwh_series[i])
                push!(wind_export_kwh_series, 0)
            end
        end
        
        analysis_runs = DataFrame(
            #identifier information
            MatchID = match_id,
            place_name = data[!, :place_name][i],
            naicsCode = data[!, :naicsCode][i],
            state = data[!, :state][i],
            DAC = data[!, :cejst_dac_int][i],
            latitude = latitude,
            longitude = longitude,
            land_space = land_acres,
            annual_electric = annual_load,
            max_wind_class_based_on_space = wind_class,
            max_num_turbines = max_cap_num_of_turbines,
            max_wind_kw_based_on_space = max_kw_based_on_space,
            new_wind_size_class = new_wind_class,
            new_wind_cap = wind_size,
            new_wind_class_num_turbines = new_wind_class == "Kept original" ? 0 : num_turbines,
            wind_production_total = sum(wind_prod_kwh_series),
            wind_production_series = [wind_prod_kwh_series],
            wind_export_total_kwh = sum(wind_export_kwh_series),
            wind_export_series_kwh = [wind_export_kwh_series],
            annual_grid_supplied_kwh = sum(grid_supplied_kwh),
            grid_supplied_kwh_series = [grid_supplied_kwh],
            wind_serving_load_kwh_series = [wind_serving_load_total_series],
            wind_serving_load_kwh_total = sum(wind_serving_load_total_series),
            ng_annual_consumption = ng_annual_mmbtu
        )

        bau_inputs1 = REopt.BAUInputs(inputs1)

        BAU_emissions_aer_total = bau_inputs1.s.site.bau_emissions_lb_CO2_per_year #this is an aggregate total electric + ng
        #println("The total BAU emissions for # $i is: ", BAU_emissions_aer_total)
        BAU_grid_emissions_aer_total = bau_inputs1.s.site.bau_grid_emissions_lb_CO2_per_year #this is grid specific total emissions 
        BAU_grid_emissions_aer_series = bau_inputs1.s.electric_utility.emissions_factor_series_lb_CO2_per_kwh #this is an 8760 series for grid emissions

        bau_inputs2 = REopt.BAUInputs(inputs2)

        BAU_emissions_srmer_total = bau_inputs2.s.site.bau_emissions_lb_CO2_per_year #this is an aggregate total electric + ng
        BAU_grid_emissions_srmer_total = bau_inputs2.s.site.bau_grid_emissions_lb_CO2_per_year #this is grid specific total emissions 
        BAU_grid_emissions_srmer_series = bau_inputs2.s.electric_utility.emissions_factor_series_lb_CO2_per_kwh #this is an 8760 series for grid-electric related emissions 

        #add in lrmer emissions 
        input_data_site["ElectricUtility"]["cambium_metric_col"] = "lrmer_co2e_c"

        s3 = Scenario(input_data_site)
        inputs3 = REoptInputs(s3)
        bau_inputs3 = REopt.BAUInputs(inputs3)

        BAU_emissions_lrmer_total = bau_inputs3.s.site.bau_emissions_lb_CO2_per_year #this is an aggregate total electric + ng
        BAU_grid_emissions_lrmer_total = bau_inputs3.s.site.bau_grid_emissions_lb_CO2_per_year #this is grid specific total emissions 
        BAU_grid_emissions_lrmer_series = bau_inputs3.s.electric_utility.emissions_factor_series_lb_CO2_per_kwh #this is an 8760 series for grid-electric related emissions 

        #now get the total aer emissions based on PV ground and roof sizes (kW) and prod factor series and emissions series
        srmer_emissions_delta = sum((BAU_grid_emissions_srmer_series .* wind_prod_kwh_series))
        aer_minus_srmer_w_tech_emissions_site = BAU_emissions_aer_total - srmer_emissions_delta #append 
        aer_minus_srmer_w_tech_emissions_percent_change_site = 1 - (aer_minus_srmer_w_tech_emissions_site ./ BAU_emissions_aer_total) #append
        aer_minus_srmer_w_tech_emissions_grid = BAU_grid_emissions_aer_total - srmer_emissions_delta #append 
        aer_minus_srmer_w_tech_emissions_percent_change_grid = 1 - (aer_minus_srmer_w_tech_emissions_grid ./ BAU_grid_emissions_aer_total) #append

        #get net load 
        net_load = load_vector - wind_prod_kwh_series
        if sum(net_load) < 0
            println("The total annual load is ", sum(load_vector))
            println("The total wind production is ", sum(wind_prod_kwh_series))
            println("The net load is ", sum(net_load))
        end

        #now get srmer_co2e_c metrics
        srmer_emissions_w_tech = sum(net_load .* BAU_grid_emissions_srmer_series) #this is appended to resultant df 
        srmer_emissions_w_tech_percent_change_site = 1 - (srmer_emissions_w_tech / BAU_emissions_srmer_total)  #appended 
        srmer_emissions_w_tech_percent_change_grid = 1 - (srmer_emissions_w_tech / BAU_grid_emissions_srmer_total)  #appended 

        #now get aer_gen_co2e_c metrics
        aer_emissions_w_tech = sum(net_load .* BAU_grid_emissions_aer_series) #append 
        aer_emissions_w_tech_percent_change_site = 1 - (aer_emissions_w_tech / BAU_emissions_aer_total) #append 
        aer_emissions_w_tech_percent_change_grid = 1 - (aer_emissions_w_tech / BAU_grid_emissions_aer_total) #append 
        
        #now get lrmer_co2e_c metrics
        lrmer_emissions_w_tech = sum(net_load .* BAU_grid_emissions_lrmer_series) #append 
        #println("The LRMER emissions w tech is ", lrmer_emissions_w_tech)
        #println("The BAU emissions total including NG is ", BAU_emissions_lrmer_total)
        #println("The BAU GRID only emissions are ", BAU_grid_emissions_lrmer_total)
        lrmer_emissions_w_tech_percent_change_site = 1 - (lrmer_emissions_w_tech / BAU_emissions_lrmer_total) #append 
        lrmer_emissions_w_tech_percent_change_grid = 1 - (lrmer_emissions_w_tech / BAU_grid_emissions_lrmer_total) #append 

        """
        Emissions 1 is BAU aer emissions - srmer tech (multiplying tech production w srmer rates, emissions delta variable)
        Emissions 2 is BAU aer emissions - aer net load w tech 
        Emissions 3 is BAU srmer emissions - srmer net load w tech
        """

        emissions = DataFrame(
            #identifier information
            MatchID = data[!, :MatchID][i],
            place_name = data[!, :place_name][i],
            #emissions information
            #aer minus srmer
            BAU_emissions_aer_gen_co2e_c_wo_tech_no_ng = BAU_grid_emissions_aer_total, #store BAU emissions before any tech is deployed
            emissions_srmer_from_tech = srmer_emissions_delta, #store emission difference when tech is deployed using srmer, net load
            RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech_no_ng = aer_minus_srmer_w_tech_emissions_grid, #store BAU emissions calculated with aer_gen_co2e_c - srmer emissions from net load of tech deployment
            PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech_no_ng = aer_minus_srmer_w_tech_emissions_percent_change_grid, #percent 
            RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech_w_ng = aer_minus_srmer_w_tech_emissions_site, #store BAU emissions calculated with aer_gen_co2e_c - srmer emissions from net load of tech deployment
            PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech_w_ng = aer_minus_srmer_w_tech_emissions_percent_change_site, #percent 
            #aer
            RESULT_BAU_emissions_aer_minus_aer_emissions_w_tech_no_ng = aer_emissions_w_tech, #store BAU emissions calculated with aer_gen_co2e_c - aer_gen_co2e_c emissions from net load of tech deployment 
            PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech_no_ng = aer_emissions_w_tech_percent_change_grid, #percent 
            BAU_emissions_aer_gen_co2e_c_wo_tech_w_ng = BAU_emissions_aer_total,
            PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech_w_ng = aer_emissions_w_tech_percent_change_site,
            #srmer
            BAU_emissions_srmer_co2e_wo_tech_no_ng = BAU_grid_emissions_srmer_total, #store BAU emissions if not using aer but using srmer for the base case 
            RESULT_bau_emissions_srmer_minus_srmer_emissions_w_tech_no_ng = srmer_emissions_w_tech, #emissions using srmer for base case and tech deployment case 
            PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech_no_ng = srmer_emissions_w_tech_percent_change_grid, #percent
            BAU_emissions_srmer_co2e_wo_tech_w_ng = BAU_emissions_srmer_total,
            PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech_w_ng = srmer_emissions_w_tech_percent_change_site, #percent  
            #lrmer 
            BAU_emissions_lrmer_co2e_wo_tech_no_ng = BAU_grid_emissions_lrmer_total, #store BAU emissions 
            RESULT_bau_emissions_lrmer_minus_lrmer_emissions_w_tech_no_ng = lrmer_emissions_w_tech, #emissions using lrmer for base case and tech deployment case 
            PERCENT_CHANGE_from_bau_emissions_lrmer_lrmer_emissions_w_tech_no_ng = lrmer_emissions_w_tech_percent_change_grid, #percent 
            BAU_emissions_lrmer_co2e_wo_tech_w_ng = BAU_emissions_lrmer_total,
            PERCENT_CHANGE_from_bau_emissions_lrmer_lrmer_emissions_w_tech_w_ng = lrmer_emissions_w_tech_percent_change_site #percent 
        )

        println("Completed runs number $i")
        
        write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/analysis_runs/", "$(file_name)_analysis_runs.json"), JSON.json(analysis_runs))
        write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/emissions/", "$(file_name)_emissions.json"), JSON.json(emissions))
        #write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/inputs_REopt/", "$(file_name)_inputs_REopt.json"), JSON.json(inputs2))
        write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/inputs_site/", "$(file_name)_inputs_data_site.json"), JSON.json(input_data_site))
        sleep(1.0)

        # Set up extractable json file with all inputs to put onto DataFrame
        #inputs_all = JSON.parsefile(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/inputs_REopt/", "$(file_name)_inputs_REopt.json"))
        #input_data_dic = JSON.parsefile(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/inputs_site/", "$(file_name)_inputs_data_site.json"))

        df = DataFrame(
            MatchID = analysis_runs[1, :MatchID],
            NAICS = analysis_runs[1, :naicsCode],
            name = analysis_runs[1, :place_name],
            state = analysis_runs[1, :state],
            DAC = analysis_runs[1, :DAC],
            input_Latitude = (analysis_runs[1, :latitude]), 
            input_Longitude = (analysis_runs[1, :longitude]), 
            input_landacres = (analysis_runs[1, :land_space]), 
            input_annual_electric_load_kWh = (round.(analysis_runs[1, :annual_electric], digits=0)), 
            input_annual_ng_load_mmbtu = (round.(analysis_runs[1, :ng_annual_consumption], digits=0)),
            #Windinputs
            wind_max_possible_cap = (analysis_runs[1, :max_wind_kw_based_on_space]),
            wind_max_size_class = (analysis_runs[1, :max_wind_class_based_on_space]),
            wind_max_size_class_num_turbines = (analysis_runs[1, :max_num_turbines]),
            input_wind_size_class = (analysis_runs[1, :new_wind_size_class]),
            input_wind_num_of_turbines = (analysis_runs[1, :new_wind_class_num_turbines]), 
            input_wind_total_total_cap = (analysis_runs[1, :new_wind_cap]),
            #total metrics 
            wind_energy_exported_total = (round.(analysis_runs[1, :wind_export_total_kwh], digits=0)),
            wind_annual_kwh_energy_production_avg_total = (round.(analysis_runs[1, :wind_production_total], digits=0)),
            wind_serving_total_load_kwh = (round.(analysis_runs[1, :wind_serving_load_kwh_total], digits=0)),
            #all other metrics
            Grid_Electricity_Supplied_kWh_annual = (round.(analysis_runs[1, :annual_grid_supplied_kwh], digits=0)),
            #AER minus SRMER
            BAU_Total_Annual_Emissions_lbs_CO2e_aer_gen_site = (round.(emissions[1, :BAU_emissions_aer_gen_co2e_c_wo_tech_w_ng], digits=4)),
            emissions_srmer_delta_lbs_CO2e_w_tech = (round.(emissions[1, :emissions_srmer_from_tech], digits=4)),
            Total_Annual_Emissions_lbs_CO2_aer_gen_minus_srmer_w_tech_site = (round.(emissions[1, :RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech_w_ng], digits=4)),
            Emission_Reduction_Fraction_aer_srmer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech_w_ng], digits=4)),
            BAU_Total_Annual_Emissions_lbs_CO2e_aer_gen_elec_only = (round.(emissions[1, :BAU_emissions_aer_gen_co2e_c_wo_tech_no_ng], digits=4)),
            Total_Annual_Emissions_lbs_CO2_aer_gen_minus_srmer_w_tech_elec_only = (round.(emissions[1, :RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech_no_ng], digits=4)),
            Emission_Reduction_Fraction_aer_srmer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech_no_ng], digits=4)),
            #AER 
            Emission_Reduction_Fraction_aer_aer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech_w_ng], digits=4)),
            Total_Annual_Emissions_lbs_CO2_aer_gen_w_tech_elec_only = (round.(emissions[1, :RESULT_BAU_emissions_aer_minus_aer_emissions_w_tech_no_ng], digits=4)),
            Emission_Reduction_Fraction_aer_aer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech_no_ng], digits=4)),
            #SRMER
            BAU_Total_Annual_Emissions_lbs_CO2e_srmer_site = (round.(emissions[1, :BAU_emissions_srmer_co2e_wo_tech_w_ng], digits=4)),
            Emission_Reduction_Fraction_srmer_srmer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech_w_ng], digits=4)),
            BAU_Total_annual_emissions_lbs_CO2e_srmer_elec_only = (round.(emissions[1, :BAU_emissions_srmer_co2e_wo_tech_no_ng], digits=4)),
            Total_Annual_Emissions_lbs_CO2_srmer_w_tech_site = (round.(emissions[1, :RESULT_bau_emissions_srmer_minus_srmer_emissions_w_tech_no_ng], digits=4)),
            Emission_reduction_fraction_srmer_srmer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech_no_ng], digits=4)),
            #LRMER 
            BAU_Total_Annual_Emissions_lbs_CO2e_lrmer_site = (round.(emissions[1, :BAU_emissions_lrmer_co2e_wo_tech_w_ng], digits=4)),
            Emission_Reduction_Fraction_lrmer_lrmer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_lrmer_lrmer_emissions_w_tech_w_ng], digits=4)),
            BAU_Total_annual_emissions_lbs_CO2e_lrmer_elec_only = (round.(emissions[1, :BAU_emissions_lrmer_co2e_wo_tech_no_ng], digits=4)),
            Total_Annual_Emissions_lbs_CO2_lrmer_w_tech_elec_only = (round.(emissions[1, :RESULT_bau_emissions_lrmer_minus_lrmer_emissions_w_tech_no_ng], digits=4)),
            Emission_reduction_fraction_lrmer_lrmer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_lrmer_lrmer_emissions_w_tech_no_ng], digits=4))
        )
        #println(df)
        # Define path to csv file
        file_storage_location = joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/results/", "$(file_name)_run_result.csv")

        # Write DataFrame to a new or existing CSV file
        if isfile(file_storage_location)
            # Append new data to the CSV file
            open(file_storage_location, "w") do io
                CSV.write(io, df, header=true)  # Append without headers
            end
        else
            # Write DataFrame to a new CSV file
            CSV.write(file_storage_location, df)
        end
    end
end

println("Number of cores: ", nprocs())
println("Number of workers: ", nworkers())

# each worker gets its id, process id and hostname
for i in workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    println(id, " " , pid, " ", host)
end

## Read the files
scenarios = read_csv_parcel_file(parcel_file)
match_id = scenarios[!, :MatchID]
evaluated = readdir("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/results/")

@info size(scenarios)
list_m = [5, 9, 18, 32, 38, 42, 60, 81, 92, 94, 130, 131, 132, 152, 165, 197, 250, 275, 287, 300, 362, 398, 402, 407, 431, 446, 471, 487, 499, 515, 520, 527, 532, 587, 617, 630, 638, 660, 662, 675, 695, 700, 711, 717, 727, 729, 731, 734, 737, 752, 760, 763, 774, 778, 810, 829, 837, 840, 846, 859, 874, 883, 897, 909, 911, 920, 927, 936, 967, 994, 1008, 1025, 1076, 1152, 1229, 1267, 1316, 1329, 1349, 1363, 1392, 1404, 1407, 1421, 1439, 1449, 1453, 1456, 1479, 1481, 1499, 1562, 1572, 1588, 1608, 7113, 7118, 7211, 7332, 7338, 8321, 8326, 8353, 8364, 8894, 8898, 8902, 8917, 8957, 8971]
@time pmap([728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740]) do i
    fname = string(match_id[i], "_Wind_run_result.csv")
    try
        if !(fname in evaluated)
            # Pass a vector of values for each site.
            #println(i, "and type of object is: ", typeof(i))
            @time run_site(i)
        end
    catch e
    	@info e
        @warn "Error for " scenarios[!, 1][i]
        sleep(3)
    end
end

# remove the workers
for i in workers()
    rmprocs(i)
end