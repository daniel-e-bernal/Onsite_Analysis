
using DataFrames
using Statistics
using DelimitedFiles
using CSV
using XLSX
using REopt
using FilePaths
using Distributions
using StatsBase
using JSON
using JuMP
using Xpress
using Plots

# Function to safely extract values from JSON with default value if key is missing
function safe_get(data::Dict{String, Any}, keys::Vector{String}, default::Any=0)
    try
        for k in keys
            data = data[k]
        end
        return data
    catch e
        if e isa KeyError
            return default
        else
            rethrow(e)
        end
    end
end
#function to convert Scenario to Dictionary
function scenario_to_dict(scenario::Scenario)
    field_dict = Dict()
    for field in fieldnames(typeof(scenario))
        field_value = getfield(scenario, field)
        field_dict[Symbol(field)] = field_value
    end
    return field_dict
end
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

#println(GCR_eq(41.000, 1))
"""
Below is the basic tilt function.
    According to PVWatts, the tilt is set at the Site's latitude for 
    fixed rack on roof or ground mount.
"""

function basic_tilt(latitude::Real, array_type::Integer)
    if array_type == 1
        tilt = 15
    elseif array_type == 0 || array_type == 2
        tilt = round(latitude, digits=0)
    else
        tilt = 15
    end
    return tilt
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

function power_density_converter_for_REopt(;
    value::Real, #a decimal
    incoming::String, #options are 'acres/kw', 'kw/ft2... whatever was not chosen for incoming, will be the outgoing unit
    )
    if incoming == "acres/kw"
        step_1 = 1 / value #converts acres/kw to kw/acres
        new_val = step_1 / 43560 #conversion to get kw/ft2
        return println("Conversion of ", value, incoming, " is ", new_val, " kw/ft2")
    elseif incoming == "kw/ft2"
        step_1 = value * 43560 #converts to kw/acres
        new_val = 1 / step_1 #converts to acres/kw
        return println("Conversion of ", value, incoming, " is ", new_val, " acres/kw")
    else 
        println("Function can not handle that unit... only acres/kw and kw/ft2")
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
        electric_estimated_or_random = []
        #find the MatchID in the electric loads file 
        find_match_id_row = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_e)
        if isempty(find_match_id_row)
            println("Got stuck trying to find a match file", counter)
            counter += 1
            continue
        else #so if a match was found 
            simu_id = find_match_id_row.simulationID[1]
            simu_id = lpad(simu_id, 6, '0')
            estimated_or_real = find_match_id_row.energy_estimated_or_random[i]
            append!(electric_estimated_or_random, estimated_or_real)
            println("Got the simu_id = ", simu_id, " and load is ", estimated_or_real)
        end
        ng_load = []
        ng_estimated_or_random = []
        #find the MatchID in the natural gas file
        find_match_id_row_ng = filter(row -> !ismissing(row.MatchID) && row.MatchID == match_id, file_ng)
        if isempty(find_match_id_row_ng)
            println("Did not find a MatchID for the natural gas consumption.")
            append!(ng_load, 0)
            append!(ng_estimated_or_random, "NaN")
        else #so if a match was found 
            annual_ng_mmbtu = find_match_id_row_ng.Estimated_Annual_Natural_Gas_MMBtu
            estimated_or_real_ng = find_match_id_row_ng.Natural_Gas_Estimated_E_or_Random_R[i]
            append!(ng_load, annual_ng_mmbtu)
            append!(ng_estimated_or_random, estimated_or_real_ng)
            println("Got the simu_id = ", simu_id, " and load is ", estimated_or_real)
        end

        #get the file path for the set file that contains the hourly loads 
        file_path_s = joinpath(folder_path, "$set_file$(string(counter)).csv")
        file_s = CSV.read(file_path_s, DataFrame)
        electric_hourly_load = []
        electric_estimated_or_random = []
        append!(electric_hourly_load, file_s[!, simu_id])
        return electric_hourly_load, electric_estimated_or_random, ng_load, ng_estimated_or_random
    end
end
#columns to select from csv file for parcel data
cols = [:MatchID, :naicsCode, :place_name, :latitude, :longitude, :state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, :critical_habs_2_int, :wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :rooftop_area_m2, :solarPV_ground_area, :wind_ground_area, :under_1_acre]
function read_csv_parcel_file(file_path::String)
    # Read the CSV into a DataFrame
    initial_df = CSV.read(file_path, DataFrame)
    
    #get the selected cols from the df
    df = initial_df[:, cols]
    df = df[df.state .!= "HI", :]
    df = df[df.state .!= "AK", :]
    df = df[.!ismissing.(df.latitude), :]
    df = df[.!isnan.(df.latitude), :]
    return df 
end

#Set-up inputs file for PV runs 
data_file = "solar_runs_v2.json"
input_data = JSON.parsefile("./Input Resources/$data_file")

#parcel file path in IEDO Teams 
parcel_file = "C:/Users/dbernal/OneDrive - NREL/General - IEDO Onsite Energy/Data/PNNL Parcel Land Coverage Analysis/updated_11_27_2024/LC_facility_parcels_NREL_11_27.csv"
#get data from CSV file for parcel data 
data = read_csv_parcel_file(parcel_file)

#get load data from IEDO Teams
electric_load_folder_path = "C:/Users/dbernal/OneDrive - NREL/General - IEDO Onsite Energy/Data/Facility Load Profiles/"
load_traits_text = "Load Facility Traits Set"
load_data_text = "Load Facility Set"
ng_load_folder_path = "C:/Users/dbernal/OneDrive - NREL/General - IEDO Onsite Energy/Data/Industrial Facility Annual Energy Consumption Estimates/"

#establish number of runs 
number_of_runs = collect(1:50)

#store results
analysis_runs = DataFrame()
input_data_dic = [] #to store the input_data_site
input_REopt_dic = [] #to store the inputs = REoptinputs(s)
"""
emissions = DataFrame(
    #identifier information
    MatchID = "",
    name = "",
    #emissions information
    BAU_emissions_aer_gen_co2e_c_wo_tech = 0.0, #store BAU emissions before any tech is deployed
    emissions_srmer_from_tech = 0.0, #store emission difference when tech is deployed using srmer, net load
    RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech = 0.0, #store BAU emissions calculated with aer_gen_co2e_c - srmer emissions from net load of tech deployment
    PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech = 0.0, #percent 
    RESULT_BAU_emissions_aer_minus_aer_emissions_w_tech = 0.0, #store BAU emissions calculated with aer_gen_co2e_c - aer_gen_co2e_c emissions from net load of tech deployment 
    PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech = 0.0, #percent 
    BAU_emissions_srmer_co2e_wo_tech = 0.0, #store BAU emissions if not using aer but using srmer for the base case 
    RESULT_bau_emissions_srmer_minus_srmer_emissions_w_tech = 0.0, #emissions using srmer for base case and tech deployment case 
    PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech = 0.0 #percent 
) #to store emissions (average and short run), "srmer_co2e_c" is short run, average is "aer_gen_co2e_c" """
emissions = DataFrame() 

sites_iter = eachindex(number_of_runs)
for i in sites_iter
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
    land_acres = round(land_acres / 4046.86, digits=4) #conversion from m2 to acres
    #roofspace conversion
    roof_space = data[!, :rooftop_area_m2][i]
    #roof_space = String(roof_space)
    roof_space = convert_string_to_float(roof_space)
    roof_space = round(roof_space / 10.7639, digits=4) #conversion from m2 to ft2
    println("Got data from parcel file")
    
    #get MatchID in data DataFrame to start getting other data from loads 
    match_id = data[!, :MatchID][i]
    load_vector = get_load_data(match_id=match_id, folder_path_e=electric_load_folder_path, folder_path_ng=ng_load_folder_path)[1]
    println(typeof(load_vector))

    input_data_site["ElectricLoad"]["loads_kw"] = load_vector
    input_data_site["ElectricLoad"]["annual_kwh"] = sum(load_vector)
    println(input_data_site["ElectricLoad"]["annual_kwh"])
    #srmer_co2e_c <- for emissions reduction
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
        end
    end
    println("The power density for ground mounted PV is: ", input_data_site["PV"][1]["acres_per_kw"])
    
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
    println("Max size for PV on ground (kW) based on annual load for # $i is: ", PV_ground_max_size_based_on_load)
    println("Max size for PV on roof (kW) based on annual load for # $i is: ", PV_roof_max_size_based_on_load)
    #now getting max size based on space
    PV_ground_max_size_based_on_space = land_acres * (1/PV_ground_power_density) #inputs1.max_sizes["ground_mount"] #max size based on space 
    PV_roof_max_size_based_n_space = roof_space * (PV_roof_power_density) #inputs1.max_sizes["roof_fixed"] #max size based on space 
    println("Max size for PV on ground (kW) based on ground space for # $i is: ", PV_ground_max_size_based_on_space)
    println("Max size for PV on roof (kW) based on roof space for # $i is: ", PV_roof_max_size_based_n_space)
    inputs1.max_sizes["ground_mount"] = PV_ground_max_size_based_on_space
    inputs1.max_sizes["roof_fixed"] = PV_roof_max_size_based_n_space
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
        println("PV roof max size based on load for # $i is greater than PV max size based on space.")
        #fix the minimum kW for the roof 
        input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_n_space * 0.99
        #identify the max kW for ground remainder after knowing the total load 
        roof_PV_size = PV_roof_max_size_based_n_space
        input_data_site["PV"][2]["max_kw"] = roof_PV_size
        ground_remainder = (PV_ground_max_size_based_on_load - roof_PV_size) >= PV_ground_max_size_based_on_space ? PV_ground_max_size_based_on_space : (PV_ground_max_size_based_on_load - roof_PV_size) 
        println("The max size for ground: ", ground_remainder)
        input_data_site["PV"][1]["max_kw"] = ground_remainder
        ground_PV_size = input_data_site["PV"][1]["max_kw"]
        println("Assinged max ground PV size (kW): ", ground_PV_size)
        input_data_site["PV"][1]["min_kw"] = ground_PV_size * 0.99

        input_data_site["PV"][1]["array_type"] = ground_PV_size <= 1356 ? 0 : 2 #if over 1 MW -> 2, one-axis tracking
        array_type_i = input_data_site["PV"][1]["array_type"]
        println("The array type is: ", array_type_i)

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
    println("The object type of PV_ground_prod_factor_series for # $i is: ", typeof(PV_ground_prod_factor_series))
    PV_roof_prod_factor_series = inputs2.production_factor["roof_fixed", :].data #this gets the production factor series for the roof PV 
    println("The set PV ground size (kW) for # $i is: ", ground_PV_size)
    println("The set PV roof size (kW) for # $i is: ", roof_PV_size)
    PV_ground_prod_kwh_series = PV_ground_prod_factor_series * ground_PV_size
    PV_roof_prod_kwh_series = PV_roof_prod_factor_series * roof_PV_size
    PV_production_total_kwh_series = PV_ground_prod_kwh_series + PV_roof_prod_kwh_series
    println("The object type of PV_ground_prod_kwh_series for # $i is a ", typeof(PV_ground_prod_kwh_series))
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
    append!(analysis_runs, DataFrame(
        #identifier information
        MatchID = data[sites_iter, :MatchID][i],
        place_name = data[sites_iter, :place_name][i],
        naicsCode = data[sites_iter, :naicsCode][i],
        state = data[sites_iter, :state][i],
        DAC = data[sites_iter, :cejst_dac_int][i],
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
        PV_serving_load_kwh_total = sum(PV_serving_load_total_series)
    ))

    bau_inputs1 = REopt.BAUInputs(inputs1)

    BAU_emissions_aer_total = bau_inputs1.s.site.bau_emissions_lb_CO2_per_year #this is an aggregate total 
    println("The total BAU emissions for # $i is: ", BAU_emissions_aer_total)
    BAU_emissions_aer_series = bau_inputs1.s.electric_utility.emissions_factor_series_lb_CO2_per_kwh #this is an 8760 series 

    bau_inputs2 = REopt.BAUInputs(inputs2)

    BAU_emissions_srmer_total = bau_inputs2.s.site.bau_emissions_lb_CO2_per_year #this is an aggregate total 
    BAU_emissions_srmer_series = bau_inputs2.s.electric_utility.emissions_factor_series_lb_CO2_per_kwh #this is an 8760 series

    #now get the total aer emissions based on PV ground and roof sizes (kW) and prod factor series and emissions series
    #aer_emissions_delta = sum((BAU_emissions_aer_series * PV_ground_prod_kwh_series) + (BAU_emissions_aer_series * PV_roof_prod_kwh_series))
    srmer_emissions_delta = sum((BAU_emissions_srmer_series .* PV_ground_prod_kwh_series) + (BAU_emissions_srmer_series .* PV_roof_prod_kwh_series))
    aer_minus_srmer_w_tech_emissions = BAU_emissions_aer_total - srmer_emissions_delta #append 
    aer_minus_srmer_w_tech_emissions_percent_change = 1 - (aer_minus_srmer_w_tech_emissions ./ BAU_emissions_aer_total) #append

    #get net load 
    net_load = load_vector - (PV_ground_prod_kwh_series + PV_roof_prod_kwh_series)

    #now get srmer metrics
    srmer_emissions_w_tech = sum(net_load .* BAU_emissions_srmer_series) #this is appended to resultant df 
    srmer_emissions_w_tech_percent_change = 1 - (srmer_emissions_w_tech / BAU_emissions_srmer_total)  #appended 

    #now get aer_gen_co2e_c metrics
    aer_emissions_w_tech = sum(net_load .* BAU_emissions_aer_series) #append 
    aer_emissions_w_tech_percent_change = 1 - (aer_emissions_w_tech / BAU_emissions_aer_total) #append 
    
    """
    Emissions 1 is BAU aer emissions - srmer tech (multiplying tech production w srmer rates, emissions delta variable)
    Emissions 2 is BAU aer emissions - aer net load w tech 
    Emissions 3 is BAU srmer emissions - srmer net load w tech
    """

    append!(emissions, DataFrame(
        #identifier information
        MatchID = data[sites_iter, :MatchID][i],
        place_name = data[sites_iter, :place_name][i],
        #emissions information
        BAU_emissions_aer_gen_co2e_c_wo_tech = BAU_emissions_aer_total, #store BAU emissions before any tech is deployed
        emissions_srmer_from_tech = srmer_emissions_delta, #store emission difference when tech is deployed using srmer, net load
        RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech = aer_minus_srmer_w_tech_emissions, #store BAU emissions calculated with aer_gen_co2e_c - srmer emissions from net load of tech deployment
        PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech = aer_minus_srmer_w_tech_emissions_percent_change, #percent 
        RESULT_BAU_emissions_aer_minus_aer_emissions_w_tech = aer_emissions_w_tech, #store BAU emissions calculated with aer_gen_co2e_c - aer_gen_co2e_c emissions from net load of tech deployment 
        PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech = aer_emissions_w_tech_percent_change, #percent 
        BAU_emissions_srmer_co2e_wo_tech = BAU_emissions_srmer_total, #store BAU emissions if not using aer but using srmer for the base case 
        RESULT_bau_emissions_srmer_minus_srmer_emissions_w_tech = srmer_emissions_w_tech, #emissions using srmer for base case and tech deployment case 
        PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech = srmer_emissions_w_tech_percent_change #percent 
    ))

    println("Complated runs number $i")
    #store inputs = REoptInputs(s) into the dictionary
    append!(input_REopt_dic, [deepcopy(inputs2)])
    #store results
    #append!(analysis_runs, [(input_data_site, results)])
    #store inputs 
    append!(input_data_dic, [deepcopy(input_data_site)])
end
println("Completed optimization")

#write results onto JSON file
write.("./results/PV_results.json", JSON.json(analysis_runs))
println("Successfully printed results on JSON file")
#write inputs onto JSON file
write.("./results/inputs_PV_onsite_analysis.json", JSON.json(input_data_dic))
println("Successfully printed data site inputs onto JSON dictionary file")
#write inputs onto JSON file
write.("./results/inputs_REopt.json", JSON.json(input_REopt_dic))
println("Successfully printed REopt inputs onto JSON dictionary file")

# Set up extractable json file with all inputs to put onto DataFrame
inputs_file = "inputs_REopt.json"
inputs_all = JSON.parsefile("results/$inputs_file")
println("Successfuly parsed $inputs_file from JSON file")

inputs_file2 = "inputs_PV_onsite_analysis.json"
input_data_dic = JSON.parsefile("results/$inputs_file2")
println("Successful parsed $inputs_file2 from JSON file")

df = DataFrame(
    MatchID = analysis_runs[sites_iter, :MatchID],
    NAICS = analysis_runs[sites_iter, :naicsCode],
    name = analysis_runs[sites_iter, :place_name],
    state = analysis_runs[sites_iter, :state],
    DAC = analysis_runs[sites_iter, :DAC],
    input_Latitude = [round(inputs_all[i]["s"]["site"]["latitude"], digits=4) for i in sites_iter],
    input_Longitude = [round(inputs_all[i]["s"]["site"]["longitude"], digits=4) for i in sites_iter],
    input_roofsqft = [round(inputs_all[i]["s"]["site"]["roof_squarefeet"], digits=4) for i in sites_iter],
    input_landacres = [round(inputs_all[i]["s"]["site"]["land_acres"], digits=4) for i in sites_iter],
    input_annual_electric_load_kWh = [sum(inputs_all[i]["s"]["electric_load"]["loads_kw"]) for i in sites_iter],
    #ground PV inputs
    input_PV_ground_location = [input_data_dic[i]["PV"][1]["location"] for i in sites_iter],
    input_PV_array_type_ground = [input_data_dic[i]["PV"][1]["array_type"] for i in sites_iter],
    input_PV_ground_gcr = [input_data_dic[i]["PV"][1]["gcr"] for i in sites_iter],
    input_PV_ground_tilt = [input_data_dic[i]["PV"][1]["tilt"] for i in sites_iter],
    input_PV_ground_power_density = [input_data_dic[i]["PV"][1]["acres_per_kw"] for i in sites_iter],
    PV_size_kw_ground = [round(input_data_dic[i]["PV"][1]["max_kw"], digits=4) for i in sites_iter],
    PV_annual_kwh_energy_production_avg_ground = (round.(analysis_runs[sites_iter, :PV_ground_production_total], digits=0)),
    PV_max_possible_size_kw_ground = [round(inputs_all[i]["max_sizes"]["ground_mount"], digits=4) for i in sites_iter],
    #roof PV inputs
    input_PV_roof_location = [input_data_dic[i]["PV"][2]["location"] for i in sites_iter],
    input_PV_array_type_roof = [input_data_dic[i]["PV"][2]["array_type"] for i in sites_iter],
    input_PV_roof_gcr = [input_data_dic[i]["PV"][2]["gcr"] for i in sites_iter],
    input_PV_roof_tilt = [input_data_dic[i]["PV"][2]["tilt"] for i in sites_iter],
    input_PV_roof_power_density = [input_data_dic[i]["PV"][2]["kw_per_square_foot"] for i in sites_iter],
    PV_size_kw_roof = [round(input_data_dic[i]["PV"][2]["max_kw"], digits=4) for i in sites_iter],
    PV_annual_kwh_energy_production_avg_roof = (round.(analysis_runs[sites_iter, :PV_roof_production_total], digits=0)),
    PV_max_possible_size_kw_roof = [round(inputs_all[i]["max_sizes"]["roof_fixed"], digits=4) for i in sites_iter],
    #total metrics 
    PV_energy_exported_total = (round.(analysis_runs[sites_iter, :PV_export_total_kwh], digits=0)),
    PV_annual_kwh_energy_production_avg_total = (round.(analysis_runs[sites_iter, :PV_production_total], digits=0)),
    PV_serving_total_load_kwh = (round.(analysis_runs[sites_iter, :PV_serving_load_kwh_total], digits=0)),
    #all other metrics
    Grid_Electricity_Supplied_kWh_annual = (round.(analysis_runs[sites_iter, :annual_grid_supplied_kwh], digits=0)),
    BAU_Total_Annual_Emissions_lbs_CO2e_aer_gen = (round.(emissions[sites_iter, :BAU_emissions_aer_gen_co2e_c_wo_tech], digits=4)),
    emissions_srmer_delta_lbs_CO2e_w_tech = (round.(emissions[sites_iter, :emissions_srmer_from_tech], digits=4)),
    Total_Annual_Emissions_lbs_CO2_aer_gen_minus_srmer_w_tech = (round.(emissions[sites_iter, :RESULT_bau_emissions_aer_minus_srmer_emissions_w_tech], digits=4)),
    LifeCycle_Emission_Reduction_Fraction_aer_srmer = (round.(emissions[sites_iter, :PERCENT_CHANGE_from_bau_emissions_aer_srmer_emissions_w_tech], digits=2)),
    Total_Annual_Emissions_lbs_CO2_aer_gen_w_tech = (round.(emissions[sites_iter, :RESULT_BAU_emissions_aer_minus_aer_emissions_w_tech], digits=4)),
    LifeCycle_Emission_Reduction_Fraction_aer_aer = (round.(emissions[sites_iter, :PERCENT_CHANGE_from_bau_emissions_aer_aer_emissions_w_tech], digits=2)),
    BAU_Total_Annual_Emissions_lbs_CO2e_srmer = (round.(emissions[sites_iter, :BAU_emissions_srmer_co2e_wo_tech], digits=4)),
    Total_Annual_Emissions_lbs_CO2_srmer_w_tech = (round.(emissions[sites_iter, :RESULT_bau_emissions_srmer_minus_srmer_emissions_w_tech ], digits=4)),
    LifeCycle_Emission_Reduction_Fraction_srmer_srmer = (round.(emissions[sites_iter, :PERCENT_CHANGE_from_bau_emissions_srmer_srmer_emissions_w_tech], digits=2))
)
println(df)

# Define path to xlsx file
file_storage_location = "./results/validation_run_results.xlsx"

# Check if the Excel file already exists
if isfile(file_storage_location)
    # Open the Excel file in read-write mode
    XLSX.openxlsx(file_storage_location, mode="rw") do xf
        counter = 0
        while true
            sheet_name = "Results_" * string(counter)
            try
                sheet = xf[sheet_name]
                counter += 1
            catch
                break
            end
        end
        sheet_name = "Results_" * string(counter)
        # Add new sheet
        XLSX.addsheet!(xf, sheet_name)
        # Write DataFrame to the new sheet
        XLSX.writetable!(xf[sheet_name], df)
    end
else
    # Write DataFrame to a new Excel file
    XLSX.writetable!(file_storage_location, df)
end

println("Successful write into XLSX file: file_storage_location")

"""
Below we will plot the graphs of energy production over load.
"""

month_days = [
    31, #January
    28, #February
    31, #March
    30, #April 
    31, #May
    30, #June
    31, #July 
    31, #August 
    30, #September
    31, #October 
    30, #November 
    31 #December
]
function closest_first_hour_of_sunday(hour::Int)
    #check that the hour is valid
    if hour < 1 || hour > 8760
        error("Hour must be between 1 and 8760.")
    end
    
    #calculate the day of the year and the current hour of the day
    day_of_year = div(hour - 1, 24) + 1
    hour_of_day = mod(hour - 1, 24) + 1
    
    #get the current day of the week (1 = Sunday, 2 = Monday, ..., 7 = Saturday)
    day_of_week = ((day_of_year - 1) % 7) + 1
    
    #calculate the number of days to the closest Sunday
    days_to_next_sunday = mod(8 - day_of_week, 7)
    days_to_prev_sunday = mod(day_of_week - 1, 7)
    
    #calcualte the hour intervals for the nearest Sunday before and after
    next_sunday_hour = hour + days_to_next_sunday * 24 - (hour_of_day - 1)
    prev_sunday_hour = hour - days_to_prev_sunday * 24 - (hour_of_day - 1)
    
    #determine the closest Sunday
    if next_sunday_hour > 8760
        next_sunday_hour = prev_sunday_hour  # Ensure within range
    elseif prev_sunday_hour < 1
        prev_sunday_hour = next_sunday_hour  # Ensure within range
    end
    
    #output the closest Sunday hour
    if abs(hour - next_sunday_hour) < abs(hour - prev_sunday_hour)
        return next_sunday_hour
    else
        return prev_sunday_hour
    end
end
#create function to get interval hour in year given month and day 
function get_hour_interval(month::Integer, day::Integer)
    #h1 = month #gets the month number 
    #d1 = day #gets day 
    h1 = month - 1 #this gets me the previous month number  
    #println(h1)
    h2 = collect(1: h1)
    #println(h2)
    days_so_far = 0
    for i in h2
        days_per_month = month_days[i]
        days_so_far += days_per_month
        #println(days_so_far)
    end
    days = days_so_far + day
    #println(days)
    hours = (days - 1) * 24
    return hours
end
month_days
#below we have the figures output directory
meets_50_summer = "./results/Plots/Meets More Than 50%/summer/"
meets_50_winter = "./results/Plots/Meets More Than 50%/winter/"
meets_10_summer = "./results/Plots/Meets Less Than 10%/summer week/"
meets_10_winter = "./results/Plots/Meets Less Than 10%/winter week/"
meets_90_summer = "./results/Plots/Meets More Than 90%/summer/"
meets_90_winter = "./results/Plots/Meets More Than 90%/winter/"

#sets the hours to index at 
summer_solstice_hours = get_hour_interval(6, 20)
summer_solstice_hour = closest_first_hour_of_sunday(summer_solstice_hours)

winter_solstice_hours = get_hour_interval(12, 21)
winter_solstice_hour = closest_first_hour_of_sunday(winter_solstice_hours)

iters = length(analysis_runs[!, :MatchID])

num_runs = collect(1:iters)
plots_iter = eachindex(num_runs)
for i in plots_iter

    place_name = data[!, :place_name][i]
    naics_code = Int(data[!, :naicsCode][i])

    summer_week_start = summer_solstice_hour
    summer_week_end = summer_solstice_hour + 167
    x_summer = collect(summer_week_start:summer_week_end)

    winter_week_start = winter_solstice_hour
    winter_week_end = winter_solstice_hour + 167
    x_winter = collect(winter_week_start:winter_week_end)

    consumption = sum(inputs_all[i]["s"]["electric_load"]["loads_kw"])
    production_aggregate = analysis_runs[!, :PV_serving_load_kwh_total][i]
    condition = consumption/production_aggregate #across the entire year, not necessarily just the summer solstice or winter solstice week
    if condition <= 0.10
        plots_dir_summer = meets_10_summer
        plots_dir_winter = meets_10_winter
    elseif condition >= 0.50
        plots_dir_summer = meets_50_summer
        plots_dir_winter = meets_50_winter
    elseif condition >= 0.90
        plots_dir_summer = meets_90_summer
        plots_dir_winter = meets_90_winter
    else
        return nothing
    end

    consumption_summer_week = inputs_all[i]["s"]["electric_load"]["loads_kw"][summer_week_start:summer_week_end]
    #println("The type of object is ", typeof(consumption_summer_week))
    production = analysis_runs[!, :PV_serving_load_kwh_series][i]
    #println("The type of object is ", typeof(production))
    production_summer_week = production[summer_week_start:summer_week_end]
    #println("The type of object is ", typeof(production_summer_week))

    plot(
        x_summer,
        production_summer_week,
        label="PV Serving Load",
        fillalpha=0.3,
        seriestype=:shape,
        color=:blue, #PV roof and ground production
        xlabel="Hour Interval",
        ylabel="Electricity (kW)",
        title="$place_name PV Production over Consumption Summer Week",
        titlefontsize=:10,
        grid=:true
    )
    plot!(
        x_summer,
        consumption_summer_week,
        label="Consumption",
        color=:black, #load profile
        linestyle=:dash
    )
    savefig(plots_dir_summer * "$(place_name)_naics_$(naics_code)_summer" )
    println("Completed summer plotting figure #$i")

    consumption_winter_week = inputs_all[i]["s"]["electric_load"]["loads_kw"][winter_week_start:winter_week_end]
    #println("The type of object is ", typeof(consumption_summer_week))
    production = analysis_runs[!, :PV_serving_load_kwh_series][i]
    #println("The type of object is ", typeof(production))
    production_winter_week = production[winter_week_start:winter_week_end]
    #println("The type of object is ", typeof(production_summer_week))

    plot(
        x_winter,
        production_winter_week,
        label="PV Serving Load",
        fillalpha=0.3, #transparency for the fill 
        seriestype=:shape, #fills the area below the curve
        color=:blue, #PV roof and ground production
        xlabel="Hour Interval",
        ylabel="Electricity (kW)",
        title="$place_name PV Production over Consumption Winter Week",
        titlefontsize=:10,
        grid=:true
    )
    plot!(
        x_winter,
        consumption_winter_week,
        label="Consumption",
        color=:black, #load profile
        linestyle=:dash
    )
    savefig(plots_dir_winter * "$(place_name)_naics_$(naics_code)_winter" )
    println("Completed winter plotting figure #$i")

end