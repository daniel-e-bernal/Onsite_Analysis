"""
The information below is the running variables for the parabolic troughts (PTC) CSP onsite scenarios.
    csp_type = "trough"
    facility_id = 100
    option = "C"
    peak_power = 150.0
    annual_demand = 100.0*8760
    println(string("Annual Energy Demand [Mwh]: ", string(round(Int,annual_demand))))
    available_area = 300.0
    rated_power = run_ssc_options(csp_type,facility_id,option,peak_power,annual_demand,available_area)
    print(rated_power)
"""

using CSV
using DataFrames
using DelimitedFiles
using HTTP
using JSON
using Dates
using XLSX
using FilePaths
include("./data_matching_function.jl")
include("./run_ssc.jl")
include("./emissions_functions.jl")

#columns to select from csv file for parcel data
cols = [:MatchID, :naicsCode, :place_name, :parcel_latitude, :parcel_longitude, :latitude, :longitude, 
:state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, :critical_habs_2_int, 
:wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :solarPV_ground_area, :wind_ground_area,
:earthquake_hazard, :desert_tortoise, :acec, :coastline_5mi, :critical_bird_areas]
function read_csv_parcel_file(file_path::String)
    # Read the CSV into a DataFrame
    initial_df = CSV.read(file_path, DataFrame)
    
    #get the selected cols from the df
    df = initial_df[:, cols]
    df = df[df.state .!= "PR", :]
    df = df[df.state .!= "VI", :]
    df = df[.!ismissing.(df.parcel_longitude), :]
    df = df[.!isnan.(df.parcel_longitude), :]
    df = df[.!ismissing.(df.MatchID), :]
    df = df[.!ismissing.(df.wind_ground_area) .& .!ismissing.(df.solarPV_ground_area), :]
    return df 
end

# Folder paths
#get load data from Load Profiles folder
electric_load_folder_path = "C:/Users/dbernal/Downloads/ITO Load Profiles/"
#electric_load_folder_path = "./tests/Load Profiles/"
ng_load_folder_path = "C:/Users/dbernal/OneDrive - NREL/Non-shared files/IEDO/Onsite Energy Program/Analysis Team/"

#sitelist csv
sitelist_csv_path = "./tests/Sitelist/LC_facility_parcels_NREL_05_18_25_TEST_CSP.csv"

#get data 
data = read_csv_parcel_file(sitelist_csv_path)

#create a function to run the CSP model specifying the sizing option and the site index
function run_csp(i::Int, option::String)

    #check that the option is a String and that it is one of the options A, B, or C 
    if !(option in ["A", "B", "C"])
        error("Option must be one of A, B, or C")
    end

    #get MatchID in data DataFrame to start getting other data from loads 
    match_id = data[!, :MatchID][i]
    match_id = String(match_id)
    #println("The site's Match ID is: ", match_id)

    #get load data, use_test=true is set so that we are using the load profiles within the tests folder not from the entire load profile dataset
    site_load_info = get_load_data(match_id=match_id, folder_path_e=electric_load_folder_path, folder_path_ng=ng_load_folder_path)
    load_vector = site_load_info[1]
    ng_annual_mmbtu = site_load_info[2]

    #use data 
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

    #because this is a lfr, the minimum land we need is 4.537 acres, if the land is less than that then we'll exit function
    if land_acres < 2.970 || data[!, :earthquake_hazard][i] == true || data[!, :desert_tortoise][i] == true || data[!, :acec][i] == true || data[!, :coastline_5mi][i] == true || data[!, :wind_exclusion][i] == true || data[!, :DOD_airspace_int][i] == true 
        @warn "Land area is less than 2.970 acres, exiting function."
        return nothing
    end

    #find out if there is a weather file for this site, if there is, use it, if not, download the weather file and store it
    #check if the weather file exists
    weather_file_csv = joinpath(@__DIR__, "weatherfiles", "weatherfile_$(match_id)_wdir.csv")
    if !isfile(weather_file_csv)
        #if the file does not exist, download it
        @info "Weather file for MatchID $(match_id) does not exist. Downloading..."
        get_weatherdata(lat=latitude, lon=longitude, debug=true, api_dv="", facility_id=match_id)
    else
        @info "Weather file for MatchID $(match_id) already exists."
    end


    #create the RESULT file name for this specific site + CSP type "_trough" (the other versios are "_mst" and "_fresnel")
    f_name = string("result_", match_id, "_lfr.csv")
    result_file_name = joinpath("./results/fresnel/option $option/$f_name")  

    csp_type = "fresnel"
    facility_id = match_id #this is how the weather file is pulled, aka MatchID
    # option is already defined in the function
    peak_power = maximum(load_vector / 1000) # from the load profile
    annual_demand = sum(load_vector / 1000) # from the load profile 
    println(string("Annual Energy Demand [Mwh]: ", string(round(Int,annual_demand))))
    available_area = land_acres # acres gotten above 
    #= run_ssc_options includes:
        csp_type - established as lfr or "fresnel", 
        facility_id = match_id, 
        option - established at the top of the function,
        peak_power - got from load_vector, 
        annual_demand - got from load_vector, 
        available_area = land_acres
    =#
    results_dict = run_ssc_options(csp_type, facility_id, option, peak_power, annual_demand, available_area)
    println("The length of the results_dict is: ", length(results_dict))
    #= net electricity produced series
    positive values are net production 
    negative values are net consumption
    =#
    net_electricity_produced_series = results_dict["Electricity Produced [MW]"] * 1000 # convert MWh to kW
    println("The sum of the net electricity produced series is: ", sum(net_electricity_produced_series))
    #initialize empty arrays
    export_series_kW = [] #8760 intervals
    grid_supplied_kW = [] #8760 intervals
    csp_serving_load_series = [] #8760 intervals

    #= calculate the values for each timestep such that we know A) how much energy was produced at that timestep, B) how much energy was supplied by grid, and C) how much was exported from CST to grid.
        When the net electricity produced is positive then we know the CST plant is producing more energy than consuming.
            If Load - Production > 0, then no export 
            If Load - Production < 0, then export = Production - Load or abs(Load - Production)
        When the net electricity produced is negative then we know the CST plant is consuming more energy than producing.
            Then, grid supplied = load + abs(Production) and no export
    =#
    for i in eachindex(net_electricity_produced_series)
        if net_electricity_produced_series[i] >= 0
            load_minus_production = load_vector[i] - net_electricity_produced_series[i]
            if load_minus_production >= 0
                push!(grid_supplied_kW, load_minus_production)
                push!(csp_serving_load_series, net_electricity_produced_series[i])
                push!(export_series_kW, 0)
            else # the load minus production < 0, which means the CSP plant is producing more than the load 
                push!(grid_supplied_kW, 0) #csp is producing more than load, therefore the grid does not supply any 
                push!(csp_serving_load_series, load_vector[i]) #csp is serving the entire facility's load + its own load 
                push!(export_series_kW, abs(load_minus_production)) #excess
            end
        else # net_electricity_produced_series[i] < 0
            grid_s = load_vector[i] + abs(net_electricity_produced_series[i]) #the load that the facility needs to be met + load CSP plant needs 
            push!(grid_supplied_kW, grid_s)
            push!(csp_serving_load_series, 0) #csp is not serving the load at all 
            push!(export_series_kW, 0) #since csp is not serving the load, then there is obviously no export 
        end
    end
    "Sample below of what is occurring above"
    # net_electricity_produced_series = [0, -10, -5, 10, 15, 20, 20, 15, 20, -2, -5, 0]
    # load_vector =                     [5,   5, 10, 10, 15, 15, 10, 10,  5,  5,  0, 0]
    # export =                          [0,   0,  0,  0,  0,  5, 10,  5, 15,  0,  0, 0]
    # grid supplied =                   [5,  15, 15,  0,  0,  0,  0,  0,  0,  7,  5, 0]
    # csp serving the load =            [0,   0,  0, 10, 15, 15, 10, 10,  5,  0,  0, 0]

    "Emissions calculations below"
    #first get emissions series for the location 
    emissions_dict = emissions_calc(latitude, longitude)
    if "emissions_factor_series_lb_CO2_per_kwh" in keys(emissions_dict)
        emissions_series = emissions_dict["emissions_factor_series_lb_CO2_per_kwh"]
    else
        emissions_series = emissions_dict["data_series"]
    end
    #println("The sum of the emissions series is: ", sum(emissions_series))
    #BAU emissions 
    bau_emissions_series = load_vector .* emissions_series #lbs CO2e per hour
    bau_total_annual_emissions = sum(bau_emissions_series) #lbs CO2e per year

    #now get future emissions using the net_electricity_produced_series
    future_emissions_series = grid_supplied_kW .* emissions_series #lbs CO2e per hour
    future_total_annual_emissions = sum(future_emissions_series) #lbs CO2e per year
    
    #emissions for ng = 117.03 lbs CO2e per MMBtu
    #ng_file_path = "./tests/NG Emissions/Manufacturing Parcel Data - Natural Gas Estimates.csv"
    ng_emissions = ng_annual_mmbtu[1] * 117.03 #lbs CO2e per year
    
    # Create a DataFrame to store results for this scenario
    df = DataFrame()

    # Add new columns with values (each as a single-element vector)
    df[!, :MatchID] = [match_id]
    df[!, :NAICS] = [data[!, :naicsCode][i]]
    df[!, :state] = [data[!, :state][i]]
    df[!, :DAC] = [data[!, :cejst_dac_int][i]]
    df[!, :input_Latitude] = [latitude]
    df[!, :input_Longitude] = [longitude]
    df[!, :input_roofsqft] = ["NA"]
    df[!, :input_landacres] = [land_acres]
    df[!, :input_annual_electric_load_kWh] = [sum(load_vector)]
    df[!, :input_annual_ng_load_mmbtu] = [ng_annual_mmbtu]
    df[!, :csp_type] = ["lfr"]
    df[!, :option] = [option]
    df[!, :cst_size_kWe] = [results_dict["Rated Power [MW]"] * 1000] # convert MW to kW
    df[!, :cst_solar_kWth] = [results_dict["Rated Power [MW]"] * 1000 * results_dict["Solar Multiple [-]"]]
    df[!, :cst_tes_hrs] = [12]
    df[!, :annual_kwh_energy_production] = [results_dict["Annual Electricity, Net [MWh]"] * 1000] # convert MWh to kWh
    df[!, :cst_size_acres] = [results_dict["Total Area [acre]"]]
    df[!, :solar_multiple] = [results_dict["Solar Multiple [-]"]]
    df[!, :capacity_factor] = [results_dict["Capacity Factor [-]"]]
    df[!, :csp_serving_load] = [sum(csp_serving_load_series)] 
    df[!, :Grid_Electricity_Supplied_kWh_annual] = [sum(grid_supplied_kW)] 
    df[!, :BAU_Total_Annual_Emissions_lbs_CO2e_lrmer_site] = [bau_total_annual_emissions + ng_emissions]
    df[!, :BAU_Total_Annual_Emissions_lbs_CO2e_lrmer_elec_only] = [bau_total_annual_emissions]
    df[!, :Total_Annual_Emissions_lbs_CO2e_lrmer_site] = [future_total_annual_emissions + ng_emissions]
    df[!, :Total_Annual_Emissions_lbs_CO2e_lrmer_elec_only] = [future_total_annual_emissions]

    # Save the DataFrame to CSV
    CSV.write(result_file_name, df)
    println("The scenario number $i was completed.")

end

## Read the data
scenarios = read_csv_parcel_file(sitelist_csv_path)
match_id = scenarios[!, :MatchID]
#end_run = length(match_id)
evaluatedB = readdir("./results/fresnel/option B/")
evaluatedC = readdir("./results/fresnel/option C/")

#get ids 
#task_id = ENV["SLURM_ARRAY_TASK_ID"]
#task_id_int = parse(Int, task_id)
#task_id_int = 1

#get corresponding MatchID
#match = match_id[task_id_int]
for i in 1:10
    fname = string("result_", match_id[i], "_lfr", ".csv") #change i to task_id_int
    task_id_int = i
    #fnameB = string("result_", match_id[i], "_mst", ".csv") #change i to task_id_int
    #fname3 = string("result_", match_id[i], "_fresnel", ".csv") #change i to task_id_int
    try
        if !(fname in evaluatedB) 
            @time run_csp(task_id_int, "B") # change i to task_id_int
        end
    catch e
        @info e
        @warn "Error encountered for MatchID $(match_id[task_id_int]) at scenario index $task_id_int: $(e)"
        # Adding a delay to prevent rapid retries or excessive logging during error handling.
        sleep(3)
    end


    fname = string("result_", match_id[i], "_lfr", ".csv") #change i to task_id_int
    #fnameB = string("result_", match_id[i], "_mst", ".csv") #change i to task_id_int
    #fname3 = string("result_", match_id[i], "_fresnel", ".csv") #change i to task_id_int
    try
        if !(fname in evaluatedC)
            @time run_csp(task_id_int, "C") #change i to task_id_int
        end
    catch e
        @info e
        @warn "Error encountered for MatchID $(match_id[task_id_int]) at scenario index $task_id_int: $(e)"
        # Adding a delay to prevent rapid retries or excessive logging during error handling.
        sleep(3)
    end
end