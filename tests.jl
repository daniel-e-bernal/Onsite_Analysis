"""Ready for use"""

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
    value::Real,
    incoming::String,
    outgoing::String,
    array_type::Integer #(0: Ground Mount Fixed (Open Rack); 1: Rooftop, Fixed; 2: Ground Mount 1-Axis Tracking)
    )
    
    
end

#columns to select from csv file for parcel data
cols = [:MatchID, :naicsCode, :place_name, :latitude, :longitude, :state, :airport_5km_int, :urban_areas_int, :DOD_airspace_int, :critical_habs_1_int, :critical_habs_2_int, :wind_exclusion, :cejst_dac_int, :FLD_PFS, :WFR_PFS, :rooftop_area_m2, :solarPV_ground_area, :wind_ground_area, :under_1_acre]
function read_csv_parcel_file(file_path::String)
    # Read the CSV into a DataFrame
    initial_df = CSV.read(file_path, DataFrame)
    
    #get the selected cols from the df
    df = initial_df[:, cols]
    return df 
end

#Set-up inputs file for PV runs 
data_file = "solar_runs_v2.json"
input_data = JSON.parsefile("./Input Resources/$data_file")

#parcel file path of internal
#file_name = "C:/Users/dbernal/OneDrive - NREL/Non-shared files/IEDO/Onsite Energy Program/Analysis Team/Input Resources/LC_facility_parcels_NREL_9_27.csv"
#parcel file path in IEDO Teams 
file_name_ = "C:/Users/dbernal/OneDrive - NREL/General - IEDO Onsite Energy/Data/PNNL Parcel Land Coverage Analysis/updated_9_27_2024/LC_facility_parcels_NREL_9_27.csv"
#get data from CSV file for parcel data 
data = read_csv_parcel_file(file_name_)

#establish number of runs 
number_of_runs = collect(1:1)

#store results
analysis_runs = []
input_data_dic = [] #to store the input_data_site
input_REopt_dic = [] #to store the inputs = REoptinputs(s)

sites_iter = eachindex(number_of_runs)
for i in sites_iter
    input_data_site = copy(input_data)

    #get the standard inputs 
    #conversion for latitude
    latitude = data[!, :latitude][i]
    latitude = String(latitude)
    latitude = convert_string_to_float(latitude)
    #converstion for longitude
    longitude = data[!, :longitude][i]
    longitude = String(longitude)
    longitude = convert_string_to_float(longitude)
    #conversion for land_acres
    land_acres = data[!, :solarPV_ground_area][i]
    land_acres = round(land_acres / 4046.86, digits=4) #conversion from m2 to acres
    #roofspace conversion
    roof_space = data[!, :rooftop_area_m2][i]
    roof_space = String(roof_space)
    roof_space = convert_string_to_float(roof_space)
    roof_space = round(roof_space / 10.7639, digits=4) #conversion from m2 to ft2 
    
    """ Work in progress
    #if there is greater than 1 acre of PV land available AND roofspace available then have 2 different types of PVs
    if land_acres > 0 && roof_space > 0
        input_data_site["PV"]["name"] = ["PV_ground", "PV_roof"]
    else 
        input_data_site["PV"]["name"] = "PV"
    """
    #input into the site specific data
    input_data_site["Site"]["latitude"] = latitude
    input_data_site["Site"]["longitude"] = longitude
    input_data_site["Site"]["land_acres"] = land_acres
    input_data_site["Site"]["roof_squarefeet"] = roof_space
    for name in [1, 2] #["ground_fixed", "ground_one_axis", "roof_fixed"]
        if input_data_site["PV"][name]["name"] == "ground_mount"
            #assign over 1 MW threshold of changing from ground mount FIXED to ground mount one-axis tracking 
            input_data_site["PV"][name]["array_type"] = land_acres < 6 ? 0 : 2 #if over 6 acres one-axis tracking
            array_type_i = input_data_site["PV"][name]["array_type"]
            #println(array_type_i)
            input_data_site["PV"][name]["acres_per_kw"] = power_density(array_type=array_type_i)
            input_data_site["PV"][name]["kw_per_square_foot"] = power_density(array_type=array_type_i)
            #assign gcr 
            gcr_pv = GCR_eq(latitude=latitude, array_type=array_type_i)
            input_data_site["PV"][name]["gcr"] = gcr_pv
            input_data_site["PV"][name]["tilt"] = tilt_pv(latitude=latitude, GCR=gcr_pv, array_type=array_type_i)
        elseif input_data_site["PV"][name]["name"] == "roof_fixed"
            input_data_site["PV"][name]["array_type"] = 1 #roof-fixed rack
            array_type_i = input_data_site["PV"][name]["array_type"]
            #println(array_type_i)
            input_data_site["PV"][name]["acres_per_kw"] = power_density(array_type=array_type_i)
            input_data_site["PV"][name]["kw_per_square_foot"] = power_density(array_type=array_type_i)
            #assign gcr 
            gcr_pv = GCR_eq(latitude=latitude, array_type=array_type_i)
            input_data_site["PV"][name]["gcr"] = gcr_pv
            input_data_site["PV"][name]["tilt"] = tilt_pv(latitude=latitude, GCR=gcr_pv, array_type=array_type_i)
        end
    end
    
    s = Scenario(input_data_site)
    inputs = REoptInputs(s)
    #testing
    #PV_prod = sum(inputs.s.pvs.production_factor_series) / 8760
    println("The pv production factor series is: ", inputs.s.pvs[1].production_factor_series)
    m1 = Model(optimizer_with_attributes(Xpress.Optimizer, "OUTPUTLOG" => 0))
    m2 = Model(optimizer_with_attributes(Xpress.Optimizer, "OUTPUTLOG" => 0))
    results = run_reopt([m1,m2], inputs)
    
    #store inputs = REoptInputs(s) into the dictionary
    append!(input_REopt_dic, [deepcopy(inputs)])
    s = scenario_to_dict(s)
    #store results
    append!(analysis_runs, [(input_data_site, results)])
    #store inputs 
    append!(input_data_dic, [deepcopy(input_data_site)])
end
println("Completed optimization")

#write onto JSON file
write.("./results/PV_results.json", JSON.json(analysis_runs))
println("Successfully printed results on JSON file")
#write inputs onto JSON file
write.("./results/inputs_PV_onsite_analysis", JSON.json(input_data_dic))
println("Successfully printed inputs onto JSON dictionary file")
#write inputs onto JSON file
write.("./results/inputs_REopt.json", JSON.json(input_REopt_dic))
println("Successfully printed REopt inputs onto JSON dictionary file")

# Set up extractable json file with all inputs to put onto DataFrame
inputs_file = "inputs_REopt.json"
inputs_all = JSON.parsefile("results/$inputs_file")
println("Successfuly parsed $inputs_file from JSON file")
"""
df = DataFrame(
    MatchID = data[sites_iter, :MatchID],
    NAICS = data[sites_iter, :naicsCode],
    name = data[sites_iter, :place_name],
    state = data[sites_iter, :state],
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
    PV_size_kw_ground = [round(analysis_runs[i][2]["PV"][1]["size_kw"], digits=0) for i in sites_iter],
    PV_year1_production_ground = [round(analysis_runs[i][2]["PV"][1]["year_one_energy_produced_kwh"], digits=0) for i in sites_iter],
    PV_annual_energy_production_avg_ground = [round(analysis_runs[i][2]["PV"][1]["annual_energy_produced_kwh"], digits=0) for i in sites_iter],
    PV_energy_lcoe_ground = [round(analysis_runs[i][2]["PV"][1]["lcoe_per_kwh"], digits=0) for i in sites_iter],
    PV_energy_exported_ground = [round(analysis_runs[i][2]["PV"][1]["annual_energy_exported_kwh"], digits=0) for i in sites_iter],
    #roof PV inputs
    input_PV_roof_location = [input_data_dic[i]["PV"][2]["location"] for i in sites_iter],
    input_PV_array_type_roof = [input_data_dic[i]["PV"][2]["array_type"] for i in sites_iter],
    input_PV_roof_gcr = [input_data_dic[i]["PV"][2]["gcr"] for i in sites_iter],
    input_PV_roof_tilt = [input_data_dic[i]["PV"][2]["tilt"] for i in sites_iter],
    input_PV_roof_power_density = [input_data_dic[i]["PV"][2]["acres_per_kw"] for i in sites_iter],
    PV_size_kw_roof = [round(analysis_runs[i][2]["PV"][2]["size_kw"], digits=0) for i in sites_iter],
    PV_year1_production_roof = [round(analysis_runs[i][2]["PV"][2]["year_one_energy_produced_kwh"], digits=0) for i in sites_iter],
    PV_annual_energy_production_avg_roof = [round(analysis_runs[i][2]["PV"][2]["annual_energy_produced_kwh"], digits=0) for i in sites_iter],
    PV_energy_lcoe_roof = [round(analysis_runs[i][2]["PV"][2]["lcoe_per_kwh"], digits=0) for i in sites_iter],
    PV_energy_exported_roof = [round(analysis_runs[i][2]["PV"][2]["annual_energy_exported_kwh"], digits=0) for i in sites_iter],
    #PV_energy_curtailed_roof = [sum(safe_get(analysis_runs[i][2], ["PV", "electric_curtailed_series_kw"], 0)) for i in sites_iter],
    Grid_Electricity_Supplied_kWh_annual = [round(analysis_runs[i][2]["ElectricUtility"]["annual_energy_supplied_kwh"], digits=0) for i in sites_iter],
    Total_Annual_Emissions_CO2 = [round(analysis_runs[i][2]["Site"]["annual_emissions_tonnes_CO2"], digits=4) for i in sites_iter],
    ElecUtility_Annual_Emissions_CO2 = [round(analysis_runs[i][2]["ElectricUtility"]["annual_emissions_tonnes_CO2"], digits=4) for i in sites_iter],
    BAU_Total_Annual_Emissions_CO2 = [round(analysis_runs[i][2]["Site"]["annual_emissions_tonnes_CO2_bau"], digits=4) for i in sites_iter],
    LifeCycle_Emissions_CO2 = [round(analysis_runs[i][2]["Site"]["lifecycle_emissions_tonnes_CO2"], digits=2) for i in sites_iter],
    BAU_LifeCycle_Emissions_CO2 = [round(analysis_runs[i][2]["Site"]["lifecycle_emissions_tonnes_CO2_bau"], digits=2) for i in sites_iter],
    LifeCycle_Emission_Reduction_Fraction = [round(analysis_runs[i][2]["Site"]["lifecycle_emissions_reduction_CO2_fraction"], digits=2) for i in sites_iter],
    LifeCycle_capex_costs_for_generation_techs = [round(analysis_runs[i][2]["Financial"]["lifecycle_generation_tech_capital_costs"], digits=2) for i in sites_iter],
    Initial_upfront_capex_wo_incentives = [round(analysis_runs[i][2]["Financial"]["initial_capital_costs"], digits=2) for i in sites_iter],
    Initial_upfront_capex_w_incentives = [round(analysis_runs[i][2]["Financial"]["initial_capital_costs_after_incentives"], digits=2) for i in sites_iter],
    Yr1_energy_cost_after_tax = [round(analysis_runs[i][2]["ElectricTariff"]["year_one_energy_cost_before_tax"], digits=2) for i in sites_iter],
    Yr1_demand_cost_after_tax = [round(analysis_runs[i][2]["ElectricTariff"]["year_one_demand_cost_before_tax"], digits=2) for i in sites_iter],
    Yr1_total_energy_bill_before_tax = [round(analysis_runs[i][2]["ElectricTariff"]["year_one_bill_before_tax"], digits=2) for i in sites_iter],
    Yr1_export_benefit_before_tax = [round(analysis_runs[i][2]["ElectricTariff"]["year_one_export_benefit_before_tax"], digits=2) for i in sites_iter],
    Annual_renewable_electricity_kwh = [round(analysis_runs[i][2]["Site"]["annual_renewable_electricity_kwh"], digits=2) for i in sites_iter],
    Annual_renewable_electricity_kwh_fraction = [round(analysis_runs[i][2]["Site"]["renewable_electricity_fraction"], digits=2) for i in sites_iter],
    npv = [round(analysis_runs[i][2]["Financial"]["npv"], digits=2) for i in sites_iter],
    lcc = [round(analysis_runs[i][2]["Financial"]["lcc"], digits=2) for i in sites_iter]
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

println("Successful write into XLSX file: file_storage_location")"""