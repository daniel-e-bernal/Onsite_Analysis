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
    
    using JSON, REopt, JuMP, DelimitedFiles, DataFrames, CSV, XLSX

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
    #defining the data to calculate GCR
    energy_loss_data = [
        ("5%", "bifacial", -0.560, 0.133, 40.2, 0.70, -0.00268, 0.361),
        ("5%", "monofacial", -0.550, 0.138, 43.4, 0.71, -0.00282, 0.388),
        ("10%", "bifacial",  -0.485, 0.171, 46.2, 0.72, -0.00437, 0.575),
        ("10%", "monofacial", -0.441, 0.198, 48.7, 0.74, -0.00476, 0.621),
        ("15%", "bifacial",  -0.414, 0.207, 49.9, 0.74, -0.00576, 0.762),
        ("15%", "monofacial", -0.371, 0.208, 51.5, 0.75, -0.00633, 0.825)
    ]
    #creating a dataframe from the data in energy_loss_data
    energy_loss_df = DataFrame(
        energy_loss = [row[1] for row in energy_loss_data],
        solar_type = [row[2] for row in energy_loss_data],
        P = [row[3] for row in energy_loss_data],
        k = [row[4] for row in energy_loss_data],
        α_0 = [row[5] for row in energy_loss_data],
        GCR_0 = [row[6] for row in energy_loss_data],
        m = [row[7] for row in energy_loss_data],
        b = [row[8] for row in energy_loss_data] 
    )
    #println(energy_loss_df)
    
    #(0: Ground Mount Fixed (Open Rack); 1: Rooftop, Fixed; 2: Ground Mount 1-Axis Tracking)
    function GCR_eq(;
        latitude::Real, 
        array_type::Integer, # this comes from the array_type in the REopt.pv core function 
        energy_loss::Integer = 5, # have the option of 5%, 10%, or 15%
        solar_type::String = "monofacial" #bifacial panels will produce greater amounts of energy 
        )
        #convert the necessary inputs to filter the data by
        energy_loss_str = string(energy_loss) * "%"
        solar_type = lowercase(solar_type) == "monofacial" ? "monofacial" : "bifacial"
    
        #filter the dataframe based on the energy_loss and solar_type
        row = energy_loss_df[(energy_loss_df.energy_loss .== energy_loss_str) .& 
                             (energy_loss_df.solar_type .== solar_type), :]
    
        α = latitude
        α_0 = row.α_0[1]
        k = row.k[1]
        GCR_0 = row.GCR_0[1]
        P = row.P[1]
        m = row.m[1]
        b = row.b[1]
    
        if array_type == 1 #this means roof-top fixed
            GCR = 0.4
        elseif array_type == 0 # this is ground mounted fixed
            GCR = P / (1 + exp(-k * (α - α_0))) + GCR_0  #k, α_0, GCR_0, and P are fitted parameters... α = latitude
        
        elseif array_type == 2 # this is one-axis tracking
            GCR = m*latitude + b
        else
            error("Invalid array_type. Must be (0: Ground Mount Fixed (Open Rack); 1: Rooftop, Fixed; 2: Ground Mount 1-Axis Tracking)")
        end
        return GCR  
    end
    
    """
    Below is a work in progress that is completed.
    It calculates the tilt based on optimized tilt angle from work done by Tonita et al 2023. 
    We use Tonita et al's figures and data which was sent over to Daniel Bernal.
    """
    
    #below we have the tilt adjustment factor data. The rows correspond to the specific latitude range.
    #the columns correspond to the GCR, which is the following: [1.0, 0.667, 0.50, 0.40, 0.285714, 0.20, 0.1, 0.002]
    #defining the data to calculate tilt optimum
    tilt_adjustment_factor_data = [
        (74.6973,	-60,	-60,	-49.2,	-39.5,	-25.7,	-19.2,	-13.7,	-12.3),
        (69.1169,	-55,	-55,	-42,	-33.3,	-23.3,	-17.1,	-12.2,	-10.5),
        (67.5696,	-54.6,	-53.9,	-43.7,	-34,	-24.5,	-18.6,	-13.1,	-10.7),
        (63.7467,	-50,	-50,	-42,	-33.7,	-22.5,	-16.5,	-10.7,	-7.7),
        (62.8084,	-50,	-49.7,	-38.8,	-28.8,	-20.8,	-16.9,	-9.6,	-6.6),
        (62.454,	-45,	-45,	-37.4,	-27.9,	-22,	-18.1,	-13.2,	-9.8),
        (50.4452,	-35,	-35,	-31.5,	-22.2,	-13.1,	-7.5,	-4.9,	-4.3),
        (49.8954,	-35,	-35,	-30.5,	-21.7,	-10.8,	-7.9,	-4, -3.3),
        (49.2827,	-35,	-35,	-30.4,	-22.2,	-15.6,	-13.6,	-11.7,	-11.1),
        (49.1852,	-35,	-35,	-30.4,	-21.1,	-13.7,	-9.3,	-6.1,	-5.3),
        (45.4201,	-30,	-30,	-21.8,	-14.8,	-8.5,	-4.9,	-2.7,	-2),
        (44.6509,	-30,	-30,	-20.8,	-14.7,	-9.8,	-7.9,	-5.3,	-4.6),
        (43.615,	-30,	-30,	-18.1,	-12.8,	-9.6,	-7.7,	-6.5,	-6),
        (40.4406,	-25,	-25,	-17.4,	-12.2,	-9.4,	-7.9, -6.2,	-5.6),
        (40.015,	-25,	-24.8,	-12.3,	-9.7,	-3.5,	-1.5,	0.3,	1),
        (39.7392,	-25,	-23.9,	-11.8,	-9.1,	-4.1,	-2.4,	-1,	-0.5),
        (36.8516,	-20,	-19.3,	-13.2,	-7.9,	-6,	-4.7,	-3.7,	-3.2),
        (35.1495,	-20,	-18.2,	-11.1,	-7.7,	-5.8,	-4.9,	-3.4,	-2.9),
        (35.0844,	-20,	-17.7,	-9.4,	-4.4,	-3,	-1.9,	-0.6,	-0.2),
        (33.4484,	-20,	-16.3,	-7.4,	-5.1,	-2.1,	-1,	-0.2,	0.2),
        (29.9509,	-15,	-14.7,	-7.3,	-5.7,	-3,	-2.2,	-1.4,	-0.8),
        (29.7601,	-15,	-13.3,	-8.2,	-5.7,	-4.1,	-3,	-2.1,	-1.5),
        (28.6434,	-15,	-11.7,	-4.1,	-3.6,	-0.5,	0.1,	0.8,	1.2),
        (25.7617,	-10,	-10,	-4.3,	-3.1,	-1.9,	-0.7,	0.1,	0.7),
        (25.6866,	-10,	-10,	-4.2,	-3.3,	-2.2,	-1,	-0.3,	0.3),
        (24.8091,	-10,	-7,	    -4.1,	-2.3,	-0.4,	0.1,	0.9,	1.6),
        (20.9674,	-5,	    -5,	    -3,	    -2.1,	-1.1,	-0.5,	0.1,	0.6),
        (20.6752,	-5,	    -5,	    -1.4,	-0.2,	0.7,	1.4,	2,	2.6),
        (16.7516,	-4.6,	-3.4,	-2.1,	-1.1,	0.7,	1.9,	2.7,	3.5)
    ]
    #creating a dataframe from the data in tilt_adjustment_factor_data
    tilt_adjustment_factor_df = DataFrame(
        latitude_input = [row[1] for row in tilt_adjustment_factor_data],
        GCR_1 = [row[2] for row in tilt_adjustment_factor_data], #GCR_1 = 1
        GCR_2 = [row[3] for row in tilt_adjustment_factor_data], #GCR_2 = 2/3 ~ 0.6667
        GCR_3 = [row[4] for row in tilt_adjustment_factor_data], #GCR_3 = 0.5
        GCR_4 = [row[5] for row in tilt_adjustment_factor_data], #GCR_4 = 0.4
        GCR_5 = [row[6] for row in tilt_adjustment_factor_data], #GCR_5 = 0.285714
        GCR_6 = [row[7] for row in tilt_adjustment_factor_data], #GCR_6 = 0.2
        GCR_7 = [row[8] for row in tilt_adjustment_factor_data], #GCR_7 = 0.1
        GCR_8 = [row[9] for row in tilt_adjustment_factor_data] #GCR_8 = 0.002
    )
    #println(tilt_adjustment_factor_df)
    
    gcr_data_vector = [
        (GCR=1.000, row_spacing=2, name="GCR_1"),
        (GCR=0.667, row_spacing=3, name="GCR_2"),
        (GCR=0.500, row_spacing=4, name="GCR_3"),
        (GCR=0.400, row_spacing=5, name="GCR_4"),
        (GCR=0.285714, row_spacing=7, name="GCR_5"),
        (GCR=0.200, row_spacing=10, name="GCR_6"),
        (GCR=0.100, row_spacing=20, name="GCR_7"),
        (GCR=0.002, row_spacing=1000, name="GCR_8")
    ] 
    
    #get closest GCR
    function closest_GCR(;
        GCR::Float64,
        gcr_data = gcr_data_vector
        )
        # Find the tuple in `gcr_data` that has the GCR closest to GCR
        closest_tuple = gcr_data[1]
        min_diff = abs(GCR - closest_tuple.GCR)
    
        for gcr_tuple in gcr_data
            current_diff = abs(GCR - gcr_tuple.GCR)
            if current_diff < min_diff
                min_diff = current_diff
                closest_tuple = gcr_tuple
            end
        end
    
        return closest_tuple
    end
    
    #get the optimal tilt 
    function tilt_pv(;
        latitude::Real,
        GCR::Real,
        array_type::Integer, #(0: Ground Mount Fixed (Open Rack); 1: Rooftop, Fixed; 2: Ground Mount 1-Axis Tracking)
        tilt_adjustment_factor_df = tilt_adjustment_factor_df
        )
        gcr_val = GCR
        #println(gcr_val)
        if GCR < 0 || GCR >= 1
            println("Invalid GCR value. Available options are:", GCR_cols)
            return nothing
        end
        
        if array_type == 2 #one axis tracking
            tilt = 0
            return tilt 
        elseif array_type == 1 #fixed tilt roof mounted
            tilt = 15.00
            return tilt
        elseif array_type == 0 #ground mounted fixed tilt
            # Get the closest GCR and the corresponding tilt adjustment factor
            closest_gcr_tuple = closest_GCR(GCR=gcr_val)  #this returns the tuple (GCR, row_spacing, name)
            closest_gcr_value = closest_gcr_tuple.name #provides string form of the column name that would be in the df with the adjustment factor
            #println(closest_gcr_value)
            tilt_column = Symbol(closest_gcr_value)  #generate the column name dynamically
            #println(tilt_column)
            #find the index of the latitude that is closest to the given latitude
            lat_diff = abs.(tilt_adjustment_factor_df.latitude_input .- latitude)
            #println(lat_diff)
            closest_lat_idx = argmin(lat_diff)  #get index of closest latitude
            #println(closest_lat_idx)
            #extract the adjustment factor from the dataframe for the closest latitude and GCR
            tilt_adjustment = tilt_adjustment_factor_df[closest_lat_idx, tilt_column]
            #println(tilt_adjustment)
            # Output the tilt based on the GCR adjustment factor
            tilt = latitude + tilt_adjustment
            return tilt
        end
        return tilt
    end
    
    function power_density(;
        array_type::Integer #(0: Ground Mount Fixed (Open Rack); 1: Rooftop, Fixed; 2: Ground Mount 1-Axis Tracking)
        )
        if array_type == 0 #Ground Mount Fixed, acres/kW
            power_density = 0.0031
            return power_density
        elseif  array_type == 1 #rooftop fixed, kW/ft^2
            power_density = 0.017
            return power_density
        else                #ground mount 1-axis tracking, acres/kW
            power_density = 0.0042
            return power_density
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
    cols = [:MatchID, :naicsCode, :place_name, :latitude, :longitude, :state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, :critical_habs_2_int, :wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :rooftop_area_m2, :solarPV_ground_area, :wind_ground_area, :under_1_acre]
    function read_csv_parcel_file(file_path::String)
        # Read the CSV into a DataFrame
        initial_df = CSV.read(file_path, DataFrame)
        
        #get the selected cols from the df
        df = initial_df[:, cols]
        df = df[df.state .!= "PR", :]
        df = df[df.state .!= "VI", :]
        df = df[.!ismissing.(df.latitude), :]
        df = df[.!isnan.(df.latitude), :]
        df = df[.!ismissing.(df.longitude), :]
        df = df[.!isnan.(df.longitude), :]
        df = df[.!ismissing.(df.MatchID), :]
        df = df[.!ismissing.(df.rooftop_area_m2) .& .!ismissing.(df.solarPV_ground_area), :]
        return df 
    end
    
    #select production factor series based on match_id
    function select_prod_factor(; 
        match_id::Any,
        roof_prod_factor_folder::String = pv_roof_prod_factors,
        ground_prod_factor_folder::String = pv_ground_prod_factors)
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
        get_ground_file = joinpath(ground_prod_factor_folder, "$match_id.csv")
    
        # Check if either file exists
        roof_exists = isfile(get_roof_file)
        ground_exists = isfile(get_ground_file)
    
        if !roof_exists
            error("No production factor file found for match_id: $match_id in roof_prod_factor_folder")
        end
        if !ground_exists
            error("No production factor file found for match_id: $match_id in ground_prod_factor_folder")
        end
    
        #read csv file and convert it to array, each file only has 1 column no header
        try
            data_ground = readdlm(get_ground_file, ',', Float64)  # Read only the first column
            data_ground = vec(data_ground)
            #data_ground = collect(data_ground)  # Convert DataFrame column to an array
            data_roof = readdlm(get_roof_file, ',', Float64)
            data_roof = vec(data_roof)
            #data_roof = CSV.read(get_roof_file, DataFrame)[:, 1]  # Read only the first column
            #data_roof = collect(data_roof)  # Convert DataFrame column to an array
            return data_ground, data_roof
        catch e
            error("Error reading file $get_ground_file and $get_roof_file: $e")
        end
    end

    #read directories for prod factors 
    #pv_roof_prod_dir = readdir("C:/Users/dbernal/Documents/GitHub/Public_REopt_analysis/pvwatts_roof_csvs/")
    #pv_ground_prod_dir = readdir("C:/Users/dbernal/Documents/Github/Public_REopt_analysis/pvwatts_ground_csvs/")

    # Folder paths
    #get load data from IEDO Teams
    electric_load_folder_path = "C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/Load Profiles/"
    load_traits_text = "Load Facility Traits Set"
    load_data_text = "Load Facility Set"
    ng_load_folder_path = "C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/NG Consumption/"

    # File imports
    pv_roof_prod_factors = "C:/Users/dbernal/Documents/GitHub/Public_REopt_analysis/pvwatts_roof_csvs/"
    pv_ground_prod_factors = "C:/Users/dbernal/Documents/GitHub/Public_REopt_analysis/pvwatts_ground_csvs/"

    #Set-up inputs file for PV runs 
    data_file = "solar_runs_v2.json"
    input_data = JSON.parsefile("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/Input Resources/$data_file")

    #parcel file path in IEDO Teams 
    #parcel_file = "C:/Users/dbernal/Documents/GitHub/Public_REopt_analysis/Input Resources/LC_facility_parcels_NREL_11_27.csv"
    parcel_file = "C:/Users/dbernal/OneDrive - NREL/Non-shared files/IEDO/Onsite Energy Program/Analysis Team/results/rerun_solar3_negatives_ground.csv"
    
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
        latitude = data[!, :latitude][i]
        #latitude = String(latitude)
        latitude = convert_string_to_float(latitude)
        #converstion for longitude
        longitude = data[!, :longitude][i]
        #longitude = String(longitude)
        longitude = convert_string_to_float(longitude)
        #conversion for land_acres
        land_acres = data[!, :solarPV_ground_area][i]
        if ismissing(land_acres) || isnan(land_acres)
            land_acres = 0
        else 
            land_acres = round(land_acres / 4046.86, digits=4) #conversion from m2 to acres
            println("The land acres constraint is (acres): ", land_acres)
        end
        #roofspace conversion
        roof_space = data[!, :rooftop_area_m2][i]
        if ismissing(roof_space) || isnan(roof_space)
            roof_space = 0
        elseif isa(roof_space, String) 
            roof_space = String(roof_space)
            roof_space = convert_string_to_float(roof_space)
            roof_space = round(roof_space * 10.7639, digits=4) #conversion from m2 to ft2
        else 
            roof_space = round(roof_space * 10.7639, digits=4) #conversion from m2 to ft2
        end
        
        #get MatchID in data DataFrame to start getting other data from loads 
        match_id = data[!, :MatchID][i]
        println("The site's Match ID is: ", match_id)

        #create the file name for this specific site + PV (the other version is Wind)
        file_name = string(match_id, "_PV")
        site_load_info = get_load_data(match_id=match_id, folder_path_e=electric_load_folder_path, folder_path_ng=ng_load_folder_path)
        load_vector = site_load_info[1]
        
        #println(typeof(load_vector))
        ng_annual_mmbtu = site_load_info[2]
        #println(typeof(ng_annual_mmbtu))

        input_data_site["ElectricLoad"]["loads_kw"] = load_vector
        input_data_site["ElectricLoad"]["annual_kwh"] = sum(load_vector)
        #println(input_data_site["ElectricLoad"]["annual_kwh"])
        #srmer_co2e_c <- for emissions reduction
        hourly_fuel = ng_annual_mmbtu[1] / 8760
        hourly_fuel = fill(hourly_fuel, 8760)
        input_data_site["DomesticHotWaterLoad"]["fuel_loads_mmbtu_per_hour"] = hourly_fuel

        #set variables outside for loop
        PV_ground_power_density = 0
        PV_roof_power_density = 0
        #input into the site specific data
        input_data_site["Site"]["latitude"] = latitude
        input_data_site["Site"]["longitude"] = longitude
        input_data_site["Site"]["land_acres"] = land_acres
        input_data_site["Site"]["roof_squarefeet"] = roof_space

        for name in [1, 2] #["ground_mount", "roof_fixed"]
            if input_data_site["PV"][name]["name"] == "ground_mount"
                #assign over 1 MW threshold of changing from ground mount FIXED to ground mount one-axis tracking 
                input_data_site["PV"][name]["array_type"] = land_acres < 4.2 ? 0 : 2 #if over 6 acres one-axis tracking
                array_type_i = input_data_site["PV"][name]["array_type"]
                #println(array_type_i)
                input_data_site["PV"][name]["acres_per_kw"] = power_density(array_type=array_type_i)
                #input_data_site["PV"][name]["kw_per_square_foot"] = 0 #only ground
                PV_ground_power_density =  input_data_site["PV"][name]["acres_per_kw"]
                #assign gcr 
                gcr_pv = GCR_eq(latitude=latitude, array_type=array_type_i)
                input_data_site["PV"][name]["gcr"] = gcr_pv
                input_data_site["PV"][name]["tilt"] = tilt_pv(latitude=latitude, GCR=gcr_pv, array_type=array_type_i)
                #input_data_site["PV"][name]["production_factor_series"] = select_prod_factor(match_id=match_id)[1]
                #prod_factor = select_prod_factor(match_id=match_id)[1]
                #println("The type object for the prod_factor is ", typeof(prod_factor))
            elseif input_data_site["PV"][name]["name"] == "roof_fixed"
                input_data_site["PV"][name]["array_type"] = 1 #roof-fixed rack
                array_type_i = input_data_site["PV"][name]["array_type"]
                #println(array_type_i)
                #input_data_site["PV"][name]["acres_per_kw"] = 0 #only roof
                input_data_site["PV"][name]["kw_per_square_foot"] = power_density(array_type=array_type_i)
                PV_roof_power_density = input_data_site["PV"][name]["kw_per_square_foot"]
                #assign gcr 
                gcr_pv = GCR_eq(latitude=latitude, array_type=array_type_i)
                input_data_site["PV"][name]["gcr"] = gcr_pv
                input_data_site["PV"][name]["tilt"] = tilt_pv(latitude=latitude, GCR=gcr_pv, array_type=array_type_i)
                #input_data_site["PV"][name]["production_factor_series"] = select_prod_factor(match_id=match_id)[2]
            end
        end
        #println("The power density for ground mounted PV is: ", input_data_site["PV"][1]["acres_per_kw"])
        #sleep(10)
        """ Below is attaining the REopt inputs related to aer_gen_co2e_c emissions to calculate BAU emissions."""
        input_data_site["ElectricUtility"]["cambium_metric_col"] = "aer_gen_co2e_c"
        s1 = Scenario(input_data_site)
        inputs1 = REoptInputs(s1)

        #getting max size based on annual load (using capacity factor)
        PV_prod_data = inputs1.production_factor
        PV_ground_capacity_factor = sum(PV_prod_data["ground_mount", :]) / 8760 #actual capacity factor
        PV_roof_capacity_factor = sum(PV_prod_data["roof_fixed", :]) / 8760 #actual capacity factor
        PV_ground_max_size_based_on_load = input_data_site["ElectricLoad"]["annual_kwh"] / (PV_ground_capacity_factor * 8760)
        PV_roof_max_size_based_on_load = input_data_site["ElectricLoad"]["annual_kwh"] / (PV_roof_capacity_factor * 8760)
        println("Max size for PV on ground (kW) based on annual load for # is: ", PV_ground_max_size_based_on_load)
        println("Max size for PV on roof (kW) based on annual load for # is: ", PV_roof_max_size_based_on_load)
        #now getting max size based on space
        PV_ground_max_size_based_on_space = land_acres * (1/PV_ground_power_density) #inputs1.max_sizes["ground_mount"] #max size based on space 
        PV_roof_max_size_based_n_space = roof_space * (PV_roof_power_density) #inputs1.max_sizes["roof_fixed"] #max size based on space 
        #println("Max size for PV on ground (kW) based on ground space for # $i is: ", PV_ground_max_size_based_on_space)
        #println("Max size for PV on roof (kW) based on roof space for # $i is: ", PV_roof_max_size_based_n_space)
        inputs1.max_sizes["ground_mount"] = PV_ground_max_size_based_on_space
        println("The PV ground max size based on space is (acres): ", PV_ground_max_size_based_on_space)
        inputs1.max_sizes["roof_fixed"] = PV_roof_max_size_based_n_space
        println("The PV roof max size based on space is: (acres)", PV_roof_max_size_based_n_space)
        #create global variables for PV sizes on roof and ground 
        roof_PV_size = 0
        ground_PV_size = 0
        #now re-set the sizes for PV on the roof and ground, with roof as the priority
        if PV_roof_max_size_based_on_load <= PV_roof_max_size_based_n_space
            println("PV roof max size based on load for # $i is less than or equal to PV max size based on space.")
            input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_on_load * 0.99
            roof_PV_size = PV_roof_max_size_based_on_load
            input_data_site["PV"][2]["max_kw"] = roof_PV_size
            #since the roof size for the total load was calculated above and it is smaller than the total based on roof size, then we do not have any more PV to deploy and the ground PV is 0
            ground_PV_size_remainder = 0 
            input_data_site["PV"][1]["max_kw"] = ground_PV_size_remainder
        elseif PV_roof_max_size_based_on_load > PV_roof_max_size_based_n_space
            println("PV roof max size based on load for # is greater than PV max size based on space.")
            #fix the minimum kW for the roof 
            input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_n_space * 0.99
            #identify the max kW for ground remainder after knowing the total load 
            roof_PV_size = PV_roof_max_size_based_n_space
            input_data_site["PV"][2]["max_kw"] = roof_PV_size
            roof_remainder = PV_roof_max_size_based_on_load - PV_roof_max_size_based_n_space
            ground_remainder = 0
            if PV_ground_max_size_based_on_space < PV_ground_max_size_based_on_load
                ground_remainder = roof_remainder > PV_ground_max_size_based_on_space ? PV_ground_max_size_based_on_space : roof_remainder
            else
                ground_remainder = roof_remainder > PV_ground_max_size_based_on_load ? PV_ground_max_size_based_on_load : roof_remainder
            end
            println("The max size for ground: ", ground_remainder)
            input_data_site["PV"][1]["max_kw"] = ground_remainder
            ground_PV_size = input_data_site["PV"][1]["max_kw"]
            println("Assinged max ground PV size (kW): ", ground_PV_size)
            input_data_site["PV"][1]["min_kw"] = ground_PV_size * 0.99

            input_data_site["PV"][1]["array_type"] = ground_PV_size <= 1356 ? 0 : 2 #if over 1 MW -> 2, one-axis tracking
            array_type_i = input_data_site["PV"][1]["array_type"]
            #println("The array type is: ", array_type_i)

            input_data_site["PV"][1]["acres_per_kw"] = power_density(array_type=array_type_i)
            PV_ground_power_density =  input_data_site["PV"][1]["acres_per_kw"]
            #assign gcr 
            gcr_pv = GCR_eq(latitude=latitude, array_type=array_type_i)
            input_data_site["PV"][1]["gcr"] = gcr_pv
            input_data_site["PV"][1]["tilt"] = tilt_pv(latitude=latitude, GCR=gcr_pv, array_type=array_type_i)
        else 
            println("Issue with assigning values to roof and ground PV for # $i.")
        end
        """ Below is attaining the REopt inputs related to srmer_co2e_c emissions to calculate BAU emissions."""

        input_data_site["ElectricUtility"]["cambium_metric_col"] = "srmer_co2e_c"
        
        s2 = Scenario(input_data_site)
        inputs2 = REoptInputs(s2)

        #will get PV related and necessary variables 
        PV_ground_prod_factor_series = inputs2.production_factor["ground_mount",:].data #this gets the production factor series for the ground PV
        #println("The object type of PV_ground_prod_factor_series for # $i is: ", typeof(PV_ground_prod_factor_series))
        PV_roof_prod_factor_series = inputs2.production_factor["roof_fixed", :].data #this gets the production factor series for the roof PV 
        #println("The set PV ground size (kW) for # $i is: ", ground_PV_size)
        #println("The set PV roof size (kW) for # $i is: ", roof_PV_size)
        PV_ground_prod_kwh_series = PV_ground_prod_factor_series * ground_PV_size
        PV_roof_prod_kwh_series = PV_roof_prod_factor_series * roof_PV_size
        PV_production_total_kwh_series = PV_ground_prod_kwh_series + PV_roof_prod_kwh_series
        #println("The object type of PV_ground_prod_kwh_series for # $i is a ", typeof(PV_ground_prod_kwh_series))
        PV_ground_load_minus_prod_kwh_series = load_vector - PV_ground_prod_kwh_series
        PV_roof_load_minus_prod_kwh_series = load_vector - PV_roof_prod_kwh_series
        PV_load_minus_total_prod_kwh_series = load_vector - (PV_roof_prod_kwh_series + PV_ground_prod_kwh_series)
        PV_ground_export_kwh_series = [] #8760 intervals
        PV_roof_export_kwh_series = [] #8760 intervals
        grid_supplied_kwh = [] #8760 intervals
        PV_serving_load_total_series = [] #8760 intervals 
        #place values in export series array 
        for i in eachindex(PV_ground_load_minus_prod_kwh_series)
            if PV_ground_load_minus_prod_kwh_series[i] <= 0
                append!(PV_ground_export_kwh_series, abs(PV_ground_load_minus_prod_kwh_series[i]))
            else
                append!(PV_ground_export_kwh_series, 0)
            end
        end
        for i in eachindex(PV_roof_load_minus_prod_kwh_series)
            if PV_roof_load_minus_prod_kwh_series[i] <= 0
                append!(PV_roof_export_kwh_series, abs(PV_roof_load_minus_prod_kwh_series[i]))
            else
                append!(PV_roof_export_kwh_series, 0)
            end
        end
        for i in eachindex(PV_load_minus_total_prod_kwh_series)
            if PV_load_minus_total_prod_kwh_series[i] <= 0
                append!(grid_supplied_kwh, 0)
                append!(PV_serving_load_total_series, load_vector[i])
            else
                append!(grid_supplied_kwh, abs(PV_load_minus_total_prod_kwh_series[i]))
                append!(PV_serving_load_total_series, PV_production_total_kwh_series[i])
            end
        end 
        #add export series array 
        PV_total_export_kwh_series = PV_roof_export_kwh_series + PV_ground_export_kwh_series
        analysis_runs = DataFrame(
            #identifier information
            MatchID = data[!, :MatchID][i],
            place_name = data[!, :place_name][i],
            naicsCode = data[!, :naicsCode][i],
            state = data[!, :state][i],
            DAC = data[!, :cejst_dac_int][i],
            PV_production_total = sum(PV_ground_prod_kwh_series + PV_roof_prod_kwh_series),
            PV_ground_production_total = sum(PV_ground_prod_kwh_series),
            PV_roof_production_total = sum(PV_roof_prod_kwh_series),
            PV_production_series = [(PV_ground_prod_kwh_series + PV_roof_prod_kwh_series)],
            PV_ground_production_series = [PV_ground_prod_kwh_series],
            PV_roof_production_series = [PV_roof_prod_kwh_series],
            PV_export_total_kwh = sum(PV_total_export_kwh_series),
            PV_export_series_kwh = [PV_total_export_kwh_series],
            annual_grid_supplied_kwh = sum(grid_supplied_kwh),
            grid_supplied_kwh_series = [grid_supplied_kwh],
            PV_serving_load_kwh_series = [PV_serving_load_total_series],
            PV_serving_load_kwh_total = sum(PV_serving_load_total_series),
            ng_annual_consumption = ng_annual_mmbtu,
            PV_max_possible_size_kw_ground = PV_ground_max_size_based_on_space,
            PV_max_possible_size_kw_roof = PV_roof_max_size_based_n_space
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
        srmer_emissions_delta = sum((BAU_grid_emissions_srmer_series .* PV_ground_prod_kwh_series) + (BAU_grid_emissions_srmer_series .* PV_roof_prod_kwh_series))
        aer_minus_srmer_w_tech_emissions_site = BAU_emissions_aer_total - srmer_emissions_delta #append 
        aer_minus_srmer_w_tech_emissions_percent_change_site = 1 - (aer_minus_srmer_w_tech_emissions_site ./ BAU_emissions_aer_total) #append
        aer_minus_srmer_w_tech_emissions_grid = BAU_grid_emissions_aer_total - srmer_emissions_delta #append 
        aer_minus_srmer_w_tech_emissions_percent_change_grid = 1 - (aer_minus_srmer_w_tech_emissions_grid ./ BAU_grid_emissions_aer_total) #append

        #get net load 
        net_load = load_vector - (PV_ground_prod_kwh_series + PV_roof_prod_kwh_series)

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
        
        write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/analysis_runs/", "$(file_name)_analysis_runs.json"), JSON.json(analysis_runs))
        write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/emissions/", "$(file_name)_emissions.json"), JSON.json(emissions))
        write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/inputs_REopt/", "$(file_name)_inputs_REopt.json"), JSON.json(inputs2))
        write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/inputs_site/", "$(file_name)_inputs_data_site.json"), JSON.json(input_data_site))
        sleep(1.0)

        # Set up extractable json file with all inputs to put onto DataFrame
        inputs_all = JSON.parsefile(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/inputs_REopt/", "$(file_name)_inputs_REopt.json"))
        input_data_dic = JSON.parsefile(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/inputs_site/", "$(file_name)_inputs_data_site.json"))

        df = DataFrame(
            MatchID = analysis_runs[1, :MatchID],
            NAICS = analysis_runs[1, :naicsCode],
            name = analysis_runs[1, :place_name],
            state = analysis_runs[1, :state],
            DAC = analysis_runs[1, :DAC],
            input_Latitude = [round(inputs_all["s"]["site"]["latitude"], digits=4)],
            input_Longitude = [round(inputs_all["s"]["site"]["longitude"], digits=4)],
            input_roofsqft = [round(inputs_all["s"]["site"]["roof_squarefeet"], digits=4)],
            input_landacres = [round(inputs_all["s"]["site"]["land_acres"], digits=4)],
            input_annual_electric_load_kWh = [sum(inputs_all["s"]["electric_load"]["loads_kw"])],
            input_annual_ng_load_mmbtu = (analysis_runs[1, :ng_annual_consumption]),
            #ground PV inputs
            input_PV_ground_location = [input_data_dic["PV"][1]["location"]],
            input_PV_array_type_ground = [input_data_dic["PV"][1]["array_type"]],
            input_PV_ground_gcr = [input_data_dic["PV"][1]["gcr"]],
            input_PV_ground_tilt = [input_data_dic["PV"][1]["tilt"]],
            input_PV_ground_power_density = [input_data_dic["PV"][1]["acres_per_kw"]],
            PV_size_kw_ground = [round(input_data_dic["PV"][1]["max_kw"], digits=4)],
            PV_annual_kwh_energy_production_avg_ground = (round.(analysis_runs[1, :PV_ground_production_total], digits=0)),
            PV_max_possible_size_kw_ground = (round.(analysis_runs[1, :PV_max_possible_size_kw_ground], digits=4)),
            #roof PV inputs
            input_PV_roof_location = [input_data_dic["PV"][2]["location"]],
            input_PV_array_type_roof = [input_data_dic["PV"][2]["array_type"]],
            input_PV_roof_gcr = [input_data_dic["PV"][2]["gcr"]],
            input_PV_roof_tilt = [input_data_dic["PV"][2]["tilt"]],
            input_PV_roof_power_density = [input_data_dic["PV"][2]["kw_per_square_foot"]],
            PV_size_kw_roof = [round(input_data_dic["PV"][2]["max_kw"], digits=4)],
            PV_annual_kwh_energy_production_avg_roof = (round.(analysis_runs[1, :PV_roof_production_total], digits=0)),
            PV_max_possible_size_kw_roof = (round.(analysis_runs[1, :PV_max_possible_size_kw_roof], digits=4)),
            #total metrics 
            PV_energy_exported_total = (round.(analysis_runs[1, :PV_export_total_kwh], digits=0)),
            PV_annual_kwh_energy_production_avg_total = (round.(analysis_runs[1, :PV_production_total], digits=0)),
            PV_serving_total_load_kwh = (round.(analysis_runs[1, :PV_serving_load_kwh_total], digits=0)),
            #all other metrics
            Grid_Electricity_Supplied_kWh_annual = (round.(analysis_runs[1, :annual_grid_supplied_kwh], digits=0)),
            #AER minus SRMER
            BAU_Total_Annual_Emissions_lbs_CO2e_aer_gen_site = (round.(emissions[1, :BAU_emissions_aer_gen_co2e_c_wo_tech_w_ng], digits=4)),
            emissions_srmer_delta_lbs_CO2e_w_tech = (round.(emissions[1, :emissions_srmer_from_tech], digits=4)),
            Total_Annual_Emissions_lbs_CO2_aer_gen_minus_srmer_w_tech_site = (round.(emissions[1, :RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech_w_ng], digits=4)),
            Emission_Reduction_Fraction_aer_srmer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech_w_ng], digits=2)),
            BAU_Total_Annual_Emissions_lbs_CO2e_aer_gen_elec_only = (round.(emissions[1, :BAU_emissions_aer_gen_co2e_c_wo_tech_no_ng], digits=4)),
            Total_Annual_Emissions_lbs_CO2_aer_gen_minus_srmer_w_tech_elec_only = (round.(emissions[1, :RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech_no_ng], digits=4)),
            Emission_Reduction_Fraction_aer_srmer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech_no_ng], digits=2)),
            #AER 
            Emission_Reduction_Fraction_aer_aer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech_w_ng], digits=2)),
            Total_Annual_Emissions_lbs_CO2_aer_gen_w_tech_elec_only = (round.(emissions[1, :RESULT_BAU_emissions_aer_minus_aer_emissions_w_tech_no_ng], digits=4)),
            Emission_Reduction_Fraction_aer_aer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech_no_ng], digits=2)),
            #SRMER
            BAU_Total_Annual_Emissions_lbs_CO2e_srmer_site = (round.(emissions[1, :BAU_emissions_srmer_co2e_wo_tech_w_ng], digits=4)),
            Emission_Reduction_Fraction_srmer_srmer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech_w_ng], digits=2)),
            BAU_Total_annual_emissions_lbs_CO2e_srmer_elec_only = (round.(emissions[1, :BAU_emissions_srmer_co2e_wo_tech_no_ng], digits=4)),
            Total_Annual_Emissions_lbs_CO2_srmer_w_tech_site = (round.(emissions[1, :RESULT_bau_emissions_srmer_minus_srmer_emissions_w_tech_no_ng], digits=4)),
            Emission_reduction_fraction_srmer_srmer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech_no_ng], digits=4)),
            #LRMER 
            BAU_Total_Annual_Emissions_lbs_CO2e_lrmer_site = (round.(emissions[1, :BAU_emissions_lrmer_co2e_wo_tech_w_ng], digits=4)),
            Emission_Reduction_Fraction_lrmer_lrmer_site = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_lrmer_lrmer_emissions_w_tech_w_ng], digits=2)),
            BAU_Total_annual_emissions_lbs_CO2e_lrmer_elec_only = (round.(emissions[1, :BAU_emissions_lrmer_co2e_wo_tech_no_ng], digits=4)),
            Total_Annual_Emissions_lbs_CO2_lrmer_w_tech_elec_only = (round.(emissions[1, :RESULT_bau_emissions_lrmer_minus_lrmer_emissions_w_tech_no_ng], digits=4)),
            Emission_reduction_fraction_lrmer_lrmer_elec_only = (round.(emissions[1, :PERCENT_CHANGE_from_bau_emissions_lrmer_lrmer_emissions_w_tech_no_ng], digits=4))
        )
        #println(df)
        # Define path to csv file
        file_storage_location = joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/results/", "$(file_name)_run_result.csv")

        # Write DataFrame to a new or existing CSV file
        if isfile(file_storage_location)
            # Append new data to the CSV file
            open(file_storage_location, "a") do io
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
evaluated = readdir("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/results/")

@info size(scenarios)
        
@time pmap(1:50) do i
    fname = string(match_id[i], "_PV_run_result.csv")
    try
        if !(fname in evaluated) #|| (fname in evaluated)
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