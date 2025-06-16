"""
The functions below are used to calculate emissions for the CSP model.
"""
#using ArchGDAL

using ArchGDAL
using HTTP
using DelimitedFiles
using DataFrames
using FilePaths
using CSV
using XLSX
using Dates
using JSON


function cambium_profile(; scenario::String = "Mid-case", 
                        location_type::String = "GEA Regions 2023", 
                        latitude::Real, 
                        longitude::Real,
                        start_year::Int = 2025,
                        lifetime::Int = 1,
                        metric_col::String = "lrmer_co2e_c",
                        grid_level::String = "enduse",
                        time_steps_per_hour::Int = 1,
                        load_year::Int = 2017,
                        profile_year::Int = 2017
                        )

    url = "https://scenarioviewer.nrel.gov/api/get-levelized/" # Production 
    project_uuid = "0f92fe57-3365-428a-8fe8-0afc326b3b43" # Cambium 2023 
    

    payload=Dict(
            "project_uuid" => project_uuid,
            "scenario" => scenario,
            "location_type" => location_type,  
            "latitude" => string(round(latitude, digits=3)),
            "longitude" => string(round(longitude, digits=3)), 
            "start_year" => string(start_year), # data year covers nominal year and 4 years proceeding; e.g., 2040 values cover time range starting in 2036
            "lifetime" => string(lifetime), # Integer 1 or greater (Default 25 yrs)
            "discount_rate" => "0.0", # Zero = simple average (a pwf with discount rate gets applied to projected CO2 costs, but not quantity.)
            "time_type" => "hourly", # hourly or annual
            "metric_col" => metric_col, # lrmer_co2e
            "smoothing_method" => "rolling", # rolling or none (only applicable to hourly queries). "rolling" best with TMY data; "none" best if 2012 weather data used.
            "gwp" => "100yrAR6", # Global warming potential values. Default: "100yrAR6". Options: "100yrAR5", "20yrAR5", "100yrAR6", "20yrAR6" or a custom tuple [1,10.0,100] with GWP values for [CO2, CH4, N2O]
            "grid_level" => grid_level, # enduse or busbar 
            "ems_mass_units" => "lb" # lb or kg
    )

    try
        r = HTTP.get(url; query=payload) 
        response = JSON.parse(String(r.body)) # contains response["status"]
        output = response["message"]
        # Convert from [lb/MWh] to [lb/kWh] if the metric is emissions-related
        data_series = occursin("co2", metric_col) ? output["values"] ./ 1000 : convert(Array{Float64,1}, output["values"])
        # Align day of week of emissions or clean energy and load profiles (Cambium data starts on Sundays so assuming profile_year=2017)
        data_series = align_profile_with_load_year(load_year=load_year, profile_year=profile_year, profile_data=data_series)
        
        if time_steps_per_hour > 1
            data_series = repeat(data_series, inner=time_steps_per_hour)
        end
        
        response_dict = Dict{String, Any}(
            "description" => "Hourly CO2 (or CO2e) grid emissions factors or clean energy fraction for applicable Cambium location and location_type, adjusted to align with load year $(load_year).",
            "units" => occursin("co2", metric_col) ? "Pounds emissions per kWh" : "Fraction of clean energy",
            "location" => output["location"],
            "metric_col" => output["metric_col"], 
            "data_series" => data_series 
        )
        return response_dict
    catch
        return Dict{String, Any}(
            "error" => "Could not look up Cambium $(metric_col) profile from point ($(latitude), $(longitude)). 
             Location is likely outside contiguous US or something went wrong with the Cambium API request."
        )
    end
end

function align_profile_with_load_year(; load_year::Int, profile_year::Int, profile_data::Array{<:Real,1})
    
    ef_start_day = dayofweek(Date(profile_year,1,1)) # Monday = 1; Sunday = 7
    load_start_day = dayofweek(Date(load_year,1,1)) 
    
    if ef_start_day == load_start_day
        profile_data_adj = profile_data
    else
        # Example: Emissions year = 2017; ef_start_day = 7 (Sunday). Load year = 2021; load_start_day = 5 (Fri)
        cut_days = 7+(load_start_day-ef_start_day) # Ex: = 7+(5-7) = 5 --> cut Sun, Mon, Tues, Wed, Thurs
        wrap_ts = profile_data[25:24+24*cut_days] # Ex: = profile_data[25:144] wrap Mon-Fri to end
        profile_data_adj = append!(profile_data[24*cut_days+1:end],wrap_ts) # Ex: now starts on Fri and end Fri to align with 2021 cal
    end

    return profile_data_adj
end

"""
AVERT emissions functions.
"""
function avert_emissions_profiles(; avert_region_abbr::String="", latitude::Real, longitude::Real, time_steps_per_hour::Int=1, load_year::Int=2017, avert_data_year::Int=2023)
    if avert_region_abbr == "" # Region not supplied
        avert_region_abbr, avert_meters_to_region = avert_region_abbreviation(latitude, longitude)
    else
        avert_meters_to_region = 0.0
    end
    avert_emissions_region = region_abbr_to_name(avert_region_abbr)
    if isnothing(avert_region_abbr)
        return Dict{String, Any}(
                "error"=>
                "Could not look up AVERT emissions region within 5 miles from point ($(latitude), $(longitude)).
                Location is likely invalid or well outside continental US, AK, and HI."
            )
    end

    response_dict = Dict{String, Any}(
        "avert_region_abbr" => avert_region_abbr,
        "avert_region" => avert_emissions_region,
        "units" => "Pounds emissions per kWh",
        "description" => "Regional hourly grid emissions factors for applicable EPA AVERT region, adjusted to align days of week with load year $(load_year).",
        "avert_meters_to_region" => avert_meters_to_region
    )
    for ekey in ["CO2", "NOx", "SO2", "PM25"]
        # Columns 1 and 2 do not contain AVERT region information, so skip them.
        file_df = joinpath(@__DIR__, "..", "..", "data", "emissions", "AVERT_Data", "AVERT_$(avert_data_year)_$(ekey)_lb_per_kwh.csv")
        if !(isfile(file_df))
            file_df = joinpath(@__DIR__, "tests", "emissions", "AVERT_$(avert_data_year)_$(ekey)_lb_per_kwh.csv")
        end
        avert_df = readdlm(file_df, ',')[:, 3:end]
        # Find col index for region. Row 1 does not contain AVERT data so skip that.
        emissions_profile_unadjusted = round.(avert_df[2:end,findfirst(x -> x == avert_region_abbr, avert_df[1,:])], digits=6)
        # Adjust for day of week alignment with load
        ef_profile_adjusted = align_profile_with_load_year(load_year=load_year, profile_year=avert_data_year, profile_data=emissions_profile_unadjusted) 
        # Adjust for non-hourly timesteps 
        if time_steps_per_hour > 1
            ef_profile_adjusted = repeat(ef_profile_adjusted,inner=time_steps_per_hour)
        end
        response_dict["emissions_factor_series_lb_"*ekey*"_per_kwh"] = ef_profile_adjusted
    end
    return response_dict
end

function avert_region_abbreviation(latitude, longitude)
    
    file_path = joinpath(@__DIR__, "..", "..", "data", "emissions", "AVERT_Data", "avert_4326.shp")
    if !(isfile(file_path))
        file_path = joinpath(@__DIR__, "tests", "emissions", "avert_4326.shp")
    end

    abbr = nothing
    meters_to_region = nothing

    shpfile = ArchGDAL.read(file_path)
	avert_layer = ArchGDAL.getlayer(shpfile, 0)

	point = ArchGDAL.fromWKT(string("POINT (",longitude," ",latitude,")"))
    
	for i in 1:ArchGDAL.nfeature(avert_layer)
		ArchGDAL.getfeature(avert_layer,i-1) do feature # 0 indexed
			if ArchGDAL.contains(ArchGDAL.getgeom(feature), point)
				abbr = ArchGDAL.getfield(feature,"AVERT")
                meters_to_region = 0.0;
			end
		end
	end
    if isnothing(abbr)
        @warn "Could not find AVERT region containing site latitude/longitude. Checking site proximity to AVERT regions."
    else
        return abbr, meters_to_region
    end

    another_file = joinpath(@__DIR__, "..", "..", "data", "emissions", "AVERT_Data", "avert_102008.shp")
    if !(isfile(another_file))
        another_file = joinpath(@__DIR__, "tests", "emissions", "avert_102008.shp")
    end

    shpfile = ArchGDAL.read(another_file)
    avert_102008 = ArchGDAL.getlayer(shpfile, 0)

    pt = ArchGDAL.createpoint(latitude, longitude)

    try
        fromProj = ArchGDAL.importEPSG(4326)
        toProj = ArchGDAL.importPROJ4("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
        ArchGDAL.createcoordtrans(fromProj, toProj) do transform
            ArchGDAL.transform!(pt, transform)
        end
    catch
        @warn "Could not look up AVERT emissions region closest to point ($(latitude), $(longitude)). Location is
        likely invalid or well outside continental US, AK and HI. Grid emissions assumed to be zero."
        return abbr, meters_to_region #nothing, nothing
    end

    distances = []
    for i in 1:ArchGDAL.nfeature(avert_102008)
        ArchGDAL.getfeature(avert_102008,i-1) do f # 0 indexed
            push!(distances, ArchGDAL.distance(ArchGDAL.getgeom(f), pt))
        end
    end
    
    ArchGDAL.getfeature(avert_102008,argmin(distances)-1) do feature	# 0 indexed
        meters_to_region = distances[argmin(distances)]

        if meters_to_region > 8046
            @warn "Your site location ($(latitude), $(longitude)) is more than 5 miles from the nearest AVERT region. Cannot calculate emissions."
            return abbr, meters_to_region #nothing, #
        else
            return ArchGDAL.getfield(feature,"AVERT"), meters_to_region
        end
    end
end

function region_abbr_to_name(region_abbr)
    lookup = Dict(
        "CA" => "California",
        "CENT" => "Central",
        "FL" => "Florida",
        "MIDA" => "Mid-Atlantic",
        "MIDW" => "Midwest",
        "NCSC" => "Carolinas",
        "NE" => "New England",
        "NW" => "Northwest",
        "NY" => "New York",
        "RM" => "Rocky Mountains",
        "SE" => "Southeast",
        "SW" => "Southwest",
        "TN" => "Tennessee",
        "TE" => "Texas",
        "AKGD" => "Alaska",
        "HIMS" => "Hawaii (except Oahu)",
        "HIOA" => "Hawaii (Oahu)"
    )
    return get(lookup, region_abbr, "")
end

function region_name_to_abbr(region_name)
    lookup = Dict(
        "California" => "CA",
        "Central" => "CENT",
        "Florida" => "FL",
        "Mid-Atlantic" => "MIDA",
        "Midwest" => "MIDW",
        "Carolinas" => "NCSC",
        "New England" => "NE",
        "Northwest" => "NW",
        "New York" => "NY",
        "Rocky Mountains" => "RM",
        "Southeast" => "SE",
        "Southwest" => "SW",
        "Tennessee" => "TN",
        "Texas" => "TE",
        "Alaska" => "AKGD",
        "Hawaii (except Oahu)" => "HIMS",
        "Hawaii (Oahu)" => "HIOA"
    )
    return get(lookup, region_name, "")
end

""" Emissions calculations flow """

#calculate emissions
function emissions_calc(latitude::Real, longitude::Real, avert_emissions_region::String="")
    #attempt to get cambium profile
    cambium_emissions = cambium_profile(
                                        latitude = latitude, 
                                        longitude = longitude)
    if haskey(cambium_emissions, "error")
        # Get AVERT emissions region    
        if avert_emissions_region == ""
            region_abbr, meters_to_region = avert_region_abbreviation(latitude, longitude)
            avert_emissions_region = region_abbr_to_name(region_abbr)
        else
            region_abbr = region_name_to_abbr(avert_emissions_region)
            meters_to_region = 0
        end
        #attempt to get avert profile if cambium fails
        emissisons_avert = avert_emissions_profiles(
                                latitude = latitude,
                                longitude = longitude)
        
        return emissisons_avert
    else
        return cambium_emissions
    end
end