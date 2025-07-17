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
    include("./solar_functions.jl")

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
        max_counter = 35
    
        while true
            #check counter 
            if counter == max_counter
                return 0 #exit loop
            end
    
            #we get the traits path first 
            file_path_e = joinpath(folder_path_e, "$traits_file$(string(counter)).csv")
            file_e = CSV.read(file_path_e, DataFrame)
            file_path_ng = joinpath(folder_path_ng, "Manufacturing Parcel Data - Natural Gas Estimates_updated_Mar2025.csv")
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
    cols = [:MatchID, :naicsCode, :place_name, :parcel_latitude, :parcel_longitude, :latitude, :longitude, :state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, :critical_habs_2_int, :wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :rooftop_area_m2, :solarPV_ground_area, :wind_ground_area, :under_1_acre]
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

    #read directories for prod factors 
    #pv_roof_prod_dir = readdir("C:/Users/dbernal/Documents/GitHub/Public_REopt_analysis/pvwatts_roof_csvs/")
    #pv_ground_prod_dir = readdir("C:/Users/dbernal/Documents/Github/Public_REopt_analysis/pvwatts_ground_csvs/")

    # Folder paths
    #get load data from IEDO Teams
    electric_load_folder_path = "C:/Users/dbernal/Downloads/ITO Load Profiles/"
    load_traits_text = "Load Facility Traits Set"
    load_data_text = "Load Facility Set"
    ng_load_folder_path = "C:/Users/dbernal/OneDrive - NREL/Non-shared files/IEDO/Onsite Energy Program/Analysis Team/"

    # File imports
    pv_roof_prod_factors = "C:/Users/dbernal/Documents/GitHub (esque)/Onsite Data/PVWatts_/PVWatts_/pvwatts_roof_csvs/"
    pv_ground_fixed_prod_factors = "C:/Users/dbernal/Documents/GitHub (esque)/Onsite Data/PVWatts_/PVWatts_/pvwatts_ground_fixed_csvs/"
    pv_ground_axis_prod_factors = "C:/Users/dbernal/Documents/GitHub (esque)/Onsite Data/PVWatts_/PVWatts_/pvwatts_ground_axis_csvs/"

    #Set-up inputs file for PV runs 
    data_file = "solar_runs_v2.json"
    input_data = JSON.parsefile("./Input Resources/$data_file")

    #parcel file path in IEDO Teams 
    parcel_file = "C:/Users/dbernal/Documents/Github (esque)/Onsite Data/updated_1_28_2025/LC_facility_parcels_NREL_01_28_25.csv"
    #parcel_file = "C:/Users/dbernal/OneDrive - NREL/Non-shared files/IEDO/Onsite Energy Program/Analysis Team/results/rerun_solar3_negatives_ground.csv"
    
    #get data from CSV file for parcel data 
    data = read_csv_parcel_file(parcel_file)

    #ENV["NREL_DEVELOPER_API_KEY"]="gAXbkyLjfTFEFfiO3YhkxxJ6rkufRaSktk40ho4x"
    
    # Add common site information
    function run_site(i::Int, option::String)
        
        #check if the sizing option that was entered is valid 
        if !(option in ["A", "B", "C"])
            error("Invalid sizing option. Must be 'A', 'B', or 'C'.")
        end

        #store results
        analysis_runs = DataFrame()
        emissions = DataFrame() 
        #for i in sites_iter
        input_data_site = copy(input_data)
 
        #conversion for latitude
        latitude = 0
        latitude = data[!, :parcel_latitude][i]
        #latitude = String(latitude)
        latitude = convert_string_to_float(latitude)
        #converstion for longitude
        longitude = 0
        longitude = data[!, :parcel_longitude][i]
        #longitude = String(longitude)
        longitude = convert_string_to_float(longitude)
        #conversion for land_acres
        land_acres = 0
        land_acres = data[!, :solarPV_ground_area][i]
        if ismissing(land_acres) || isnan(land_acres)
            land_acres = 0
        else 
            land_acres = round(land_acres / 4046.86, digits=4) #conversion from m2 to acres
            println("The land acres constraint is (acres): ", land_acres)
        end
        #roofspace conversion
        roof_space = 0
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
        
        ng_annual_mmbtu = site_load_info[2]

        input_data_site["ElectricLoad"]["loads_kw"] = load_vector
        input_data_site["ElectricLoad"]["annual_kwh"] = sum(load_vector)
        #println(input_data_site["ElectricLoad"]["annual_kwh"])
        #srmer_co2e_c <- for emissions reduction
        hourly_fuel = ng_annual_mmbtu[1] / 8760
        hourly_fuel = fill(hourly_fuel, 8760)
        input_data_site["DomesticHotWaterLoad"]["fuel_loads_mmbtu_per_hour"] = hourly_fuel

        #attain production factors for the site
        prod_factors = select_prod_factor(match_id=match_id)
        roof_fixed_prod_factor_series = prod_factors[1]
        ground_fixed_prod_factor_series = prod_factors[2]
        ground_axis_prod_factor_series = prod_factors[3]

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
                input_data_site["PV"][name]["array_type"] = land_acres < 4.2 ? 0 : 2 #if over 4.2 acres one-axis tracking
                array_type_i = input_data_site["PV"][name]["array_type"]
                #println(array_type_i)
                input_data_site["PV"][name]["acres_per_kw"] = power_density(array_type=array_type_i)
                #input_data_site["PV"][name]["kw_per_square_foot"] = 0 #only ground
                PV_ground_power_density =  input_data_site["PV"][name]["acres_per_kw"]
                #assign gcr 
                gcr_pv = GCR_eq(latitude=latitude, array_type=array_type_i)
                input_data_site["PV"][name]["gcr"] = gcr_pv
                input_data_site["PV"][name]["tilt"] = tilt_pv(latitude=latitude, GCR=gcr_pv, array_type=array_type_i)
                input_data_site["PV"][name]["production_factor_series"] = ground_axis_prod_factor_series == 0 ? ground_fixed_prod_factor_series : ground_axis_prod_factor_series
                #println("The type object for the prod_factor is ", typeof(prod_factor))
            elseif input_data_site["PV"][name]["name"] == "roof_fixed"
                input_data_site["PV"][name]["array_type"] = 1 #roof-fixed rack
                #input_data_site["PV"][name]["acres_per_kw"] = 0 #only roof
                input_data_site["PV"][name]["kw_per_square_foot"] = power_density(array_type=1)
                PV_roof_power_density = input_data_site["PV"][name]["kw_per_square_foot"]
                #assign gcr 
                gcr_pv = GCR_eq(latitude=latitude, array_type=1)
                input_data_site["PV"][name]["gcr"] = gcr_pv
                input_data_site["PV"][name]["tilt"] = tilt_pv(latitude=latitude, GCR=gcr_pv, array_type=1)
                input_data_site["PV"][name]["production_factor_series"] = roof_fixed_prod_factor_series
            end
        end
        println("Sized PV")
        """ Below is attaining the REopt inputs related to aer_gen_co2e_c emissions to calculate BAU emissions."""
        input_data_site["ElectricUtility"]["cambium_metric_col"] = "aer_gen_co2e_c"
        s1 = Scenario(input_data_site)
        inputs1 = REoptInputs(s1)

        #getting the pv sizing options now 
        pv_sizing_results = pv_size_option(sizing_option=option, 
            roof_prod_series=roof_fixed_prod_factor_series,
            ground_fixed_prod_series=ground_fixed_prod_factor_series,
            ground_axis_prod_series=ground_axis_prod_factor_series,
            baseload=baseload_kw,
            load_kw=load_vector,
            land_avail=land_acres,
            roof_avail=roof_space,
            data_dict=input_data_site
        )
        roof_PV_size = pv_sizing_results[1]
        ground_PV_size = pv_sizing_results[2]
        one_axis_deployment = pv_sizing_results[3]
        input_data_site = pv_sizing_results[4]
        """ Below is attaining the REopt inputs related to srmer_co2e_c emissions to calculate BAU emissions."""
        println("emissions start")
        input_data_site["ElectricUtility"]["cambium_metric_col"] = "srmer_co2e_c"
        
        s2 = Scenario(input_data_site)
        inputs2 = REoptInputs(s2)
        println("Successfully ran inputs2 for $i")
        PV_ground_prod_factor_series = one_axis_deployment == false ? ground_fixed_prod_factor_series : ground_axis_prod_factor_series

        #calculate actual production series (kWh)
        PV_ground_prod_kwh_series = PV_ground_prod_factor_series * ground_PV_size #8760 series
        PV_roof_prod_kwh_series = roof_fixed_prod_factor_series * roof_PV_size #8760 series
        PV_production_total_kwh_series = PV_ground_prod_kwh_series + PV_roof_prod_kwh_series #8760 series
        println("Successfully calculated kwh series for $i")
        #difference between load to production 
        PV_load_minus_total_prod_kwh_series = load_vector - PV_production_total_kwh_series #8760 series
        println("Successfully got the net load - PV load minus total production for $i")
        #initialize empty arrays
        PV_export_kwh_series = [] #8760 intervals
        grid_supplied_kwh = [] #8760 intervals
        PV_serving_load_total_series = [] #8760 intervals 

        #calculate values in export series array 
        for i in eachindex(PV_load_minus_total_prod_kwh_series)
            if PV_load_minus_total_prod_kwh_series[i] <= 0
                push!(grid_supplied_kwh, 0)
                push!(PV_serving_load_total_series, load_vector[i])
                excess = PV_production_total_kwh_series[i] - load_vector[i]
                push!(PV_export_kwh_series, excess)
            else
                grid_s = load_vector[i] - PV_production_total_kwh_series[i]
                push!(grid_supplied_kwh, grid_s)
                push!(PV_serving_load_total_series, PV_production_total_kwh_series[i])
                push!(PV_export_kwh_series, 0)
            end
        end 
        println("Successfully calculated export and grid supplied kWh series for $i")
        analysis_runs = DataFrame(
            #identifier information
            MatchID = data[!, :MatchID][i],
            place_name = data[!, :place_name][i],
            naicsCode = data[!, :naicsCode][i],
            state = data[!, :state][i],
            DAC = data[!, :cejst_dac_int][i],
            latitude = latitude,
            longitude = longitude,
            land_space = land_acres,
            roof_sqft = roof_space,
            annual_elec_load = sum(load_vector),
            PV_production_total = sum(PV_production_total_kwh_series),
            PV_ground_production_total = sum(PV_ground_prod_kwh_series),
            PV_roof_production_total = sum(PV_roof_prod_kwh_series),
            PV_production_series = [PV_production_total_kwh_series],
            PV_ground_production_series = [PV_ground_prod_kwh_series],
            PV_roof_production_series = [PV_roof_prod_kwh_series],
            PV_export_total_kwh = sum(PV_export_kwh_series),
            PV_export_series_kwh = [PV_export_kwh_series],
            annual_grid_supplied_kwh = sum(grid_supplied_kwh),
            grid_supplied_kwh_series = [grid_supplied_kwh],
            PV_serving_load_kwh_series = [PV_serving_load_total_series],
            PV_serving_load_kwh_total = sum(PV_serving_load_total_series),
            ng_annual_consumption = ng_annual_mmbtu,
            PV_max_possible_size_kw_ground = PV_ground_max_size_based_on_space,
            PV_max_possible_size_kw_roof = PV_roof_max_size_based_n_space
        )
        println("Successfully created df for analysis_runs for $i")
        bau_inputs1 = REopt.BAUInputs(inputs1)

        BAU_emissions_aer_total = bau_inputs1.s.site.bau_emissions_lb_CO2_per_year #this is an aggregate total electric + ng
        println("The total BAU emissions for # $i is: ", BAU_emissions_aer_total)
        BAU_grid_emissions_aer_total = bau_inputs1.s.site.bau_grid_emissions_lb_CO2_per_year #this is grid specific total emissions 
        BAU_grid_emissions_aer_series = bau_inputs1.s.electric_utility.emissions_factor_series_lb_CO2_per_kwh #this is an 8760 series for grid emissions

        bau_inputs2 = REopt.BAUInputs(inputs2)
        println("Successfully ran inputs2 for $i")
        BAU_emissions_srmer_total = bau_inputs2.s.site.bau_emissions_lb_CO2_per_year #this is an aggregate total electric + ng
        BAU_grid_emissions_srmer_total = bau_inputs2.s.site.bau_grid_emissions_lb_CO2_per_year #this is grid specific total emissions 
        BAU_grid_emissions_srmer_series = bau_inputs2.s.electric_utility.emissions_factor_series_lb_CO2_per_kwh #this is an 8760 series for grid-electric related emissions 

        #add in lrmer emissions 
        input_data_site["ElectricUtility"]["cambium_metric_col"] = "lrmer_co2e_c"

        s3 = Scenario(input_data_site)
        inputs3 = REoptInputs(s3)
        bau_inputs3 = REopt.BAUInputs(inputs3)
        println("Successfuly got bau_inputs 3 for $match_id")
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
        net_load = PV_load_minus_total_prod_kwh_series #basically load_vector - (PV_ground_prod_kwh_series + PV_roof_prod_kwh_series)

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
        
        write(joinpath("./results/PV/analysis_runs/", "$(file_name)_analysis_runs.json"), JSON.json(analysis_runs))
        write(joinpath("./results/PV/emissions/", "$(file_name)_emissions.json"), JSON.json(emissions))
        #write(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/inputs_REopt/", "$(file_name)_inputs_REopt.json"), JSON.json(inputs2))
        write(joinpath("./results/PV/inputs_site/", "$(file_name)_inputs_data_site.json"), JSON.json(input_data_site))
        #sleep(1.0)

        # Set up extractable json file with all inputs to put onto DataFrame
        #inputs_all = JSON.parsefile(joinpath("C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/PV/inputs_REopt/", "$(file_name)_inputs_REopt.json"))
        input_data_dic = JSON.parsefile(joinpath("./results/PV/inputs_site/", "$(file_name)_inputs_data_site.json"))

        df = DataFrame(
            MatchID = analysis_runs[1, :MatchID],
            NAICS = analysis_runs[1, :naicsCode],
            name = analysis_runs[1, :place_name],
            state = analysis_runs[1, :state],
            DAC = analysis_runs[1, :DAC],
            input_Latitude = ((analysis_runs[1, :latitude])),
            input_Longitude = ((analysis_runs[1, :longitude])),
            input_roofsqft = (round.(analysis_runs[1, :roof_sqft], digits=2)),
            input_landacres = (round.(analysis_runs[1, :land_space], digits=4)),
            input_annual_electric_load_kWh = (round.(analysis_runs[1, :annual_elec_load], digits=0)),
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
            alignment_factor_for_tech_serving_load = ((round.(analysis_runs[1, :PV_serving_load_kwh_total], digits=0)) / (round.(analysis_runs[1, :annual_elec_load], digits=0))),
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
        file_storage_location = joinpath("./results/PV/results/", "$(file_name)_run_result.csv")

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
evaluated = readdir("./results/PV/results/")

@info size(scenarios)
        
@time pmap(1:30) do i
    fname = string(match_id[i], "_PV_run_result.csv")
    try
        if !(fname in evaluated) #|| (fname in evaluated)
            # Pass a vector of values for each site.
            #println(i, "and type of object is: ", typeof(i))
            @time run_site(i, "B")
            #sleep(1)
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