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

using CSV, DataFrames, DelimitedFiles

#columns to select from csv file for parcel data
cols = [:MatchID, :naicsCode, :place_name, :parcel_latitude, :parcel_longitude, :latitude, :longitude, 
:state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, :critical_habs_2_int, 
:wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :solarPV_ground_area, :wind_ground_area]
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
electric_load_folder_path = "./tests/Load Profiles/"
ng_load_folder_path = "C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/NG Consumption/"

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

    #df to store results
    analysis_runs = DataFrame()
    emissions = DataFrame()

    #get load data 

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

    #SAM takes in area in acres 
    #


    #get MatchID in data DataFrame to start getting other data from loads 
    match_id = data[!, :MatchID][i]
    #println("The site's Match ID is: ", match_id)

    #create the file name for this specific site + CSP type "_ptc" (the other versios are "_pt" and "_lfr")
    file_name = string(match_id, "_ptc")
    site_load_info = get_load_data(match_id=match_id, folder_path_e=electric_load_folder_path, folder_path_ng=ng_load_folder_path)
    load_vector = site_load_info[1]
    
    #leave for later
    ng_annual_mmbtu = site_load_info[2]

    csp_type = "trough"
    facility_id = match_id #this is how the weather file is pulled, aka MatchID
    # option is already defined in the function
    peak_power = maximum(load_vector / 1000) # from the load profile
    annual_demand = sum(load_vector / 1000) # from the load profile 100.0*8760
    println(string("Annual Energy Demand [Mwh]: ", string(round(Int,annual_demand))))
    available_area = land_acres # acres gotten above 
    #= run_ssc_options includes:
        csp_type - established as ptc or "trough", 
        facility_id = match_id, 
        option - established at the top of the function,
        peak_power - got from load_vector, 
        annual_demand - got from load_vector, 
        available_area = land_acres
    =#
    rated_power = run_ssc_options(csp_type, facility_id, option, peak_power, annual_demand, land_acres)
    print(rated_power)

end

## Read the data
scenarios = read_csv_parcel_file(parcel_file)
match_id = scenarios[!, :MatchID]
evaluated = readdir("./results/trough/")

#get ids 
task_id = ENV["SLURM_ARRAY_TASK_ID"]
task_id_int = parse(Int, task_id)

# Loop through the scenarios
fname = string(match_id[task_id_int], ".csv")
try
    if !(fname in evaluated)
        @time run_csp(task_id_int, "B")
        @time run_csp(task_id_int, "C")
    end
catch e
    @info e
    @warn "Error for " scenarios[!, 1][i]
    sleep(3)
end