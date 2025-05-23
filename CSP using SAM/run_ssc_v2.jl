### Script for running any type of SSC module
using JSON
import HTTP
using DelimitedFiles
using DataFrames
using CSV
using Base.Iterators

function get_weatherdata(lat::Float64,lon::Float64,debug::Bool)
    ### Call NSRDB
    # api_jgifford = "wKt35uq0aWoNHnzuwbcUxElPhVuo0K18YPSgZ9Ph"
    api_jgifford = "IOEFhaZoaOxkB2l9XvkaFswXgTsJxRqob1hBbPdv" #as of 9/13/2024 2pm
    attributes_tmy = "ghi,dhi,dni,wind_speed,wind_direction,air_temperature,surface_pressure,dew_point"
    url = string("http://developer.nrel.gov/api/nsrdb/v2/solar/psm3-2-2-tmy-download.csv?api_key=",api_jgifford,
        "&wkt=POINT(",lon,"%20",lat,")&attributes=",attributes_tmy,
        "&names=tmy&utc=false&leap_day=true&interval=60&email=jeffrey.gifford@nrel.gov")
    # r = HTTP.request("GET", url)
    
    df = DataFrame(CSV.File(HTTP.get(url).body,delim=",",silencewarnings=true))
    
    ### Write csv file for checking (can comment out/delete when not debugging)
    debug = false
    if debug
        weatherfile_string = string("weatherfile_",lat,"_",lon,"_wdir.csv")
        CSV.write(weatherfile_string,df)
    end
    
    ### Create weather data dataframe for SAM
    weatherdata = Dict()
    weatherdata["tz"] = parse(Int64,df."Time Zone"[1])
    weatherdata["elev"] = parse(Float64,df."Elevation"[1])
    weatherdata["lat"] = parse(Float64,df."Latitude"[1])
    weatherdata["lon"] = parse(Float64,df."Longitude"[1])
    new_df = vcat(df[3:end, :])
    weatherdata["year"] = parse.(Int64,new_df."Source") # Source --> year 
    weatherdata["month"] = parse.(Int64,new_df."Location ID") # Location ID --> month
    weatherdata["day"] = parse.(Int64,new_df."City") # City --> day 
    weatherdata["hour"] = parse.(Int64,new_df."State") # State --> hour
    weatherdata["minute"] = parse.(Int64,new_df."Country") # Country --> minute
    weatherdata["dn"] = parse.(Float64,new_df."Time Zone") # Time Zone --> dn (DNI)
    weatherdata["df"] = parse.(Float64,new_df."Longitude") # Longitude --> df (DHI)
    weatherdata["gh"] = parse.(Float64,new_df."Latitude") # Latitude --> gh (GHI)
    weatherdata["wspd"] = parse.(Float64,new_df."Elevation") # Elevation --> wspd
    weatherdata["wdir"] = parse.(Int64,new_df."Local Time Zone") # Local Time Zone --> wdir
    weatherdata["tdry"] = parse.(Float64,new_df."Dew Point Units") # Dew Point Units --> tdry
    weatherdata["tdew"] = parse.(Float64,new_df."DNI Units") # Clearsky DNI Units --> rhum (RH)
    weatherdata["pres"] = parse.(Float64,new_df."DHI Units") # Clearsky DHI Units --> pres
    ### Full list of weather data types (not all are required)
    # (numbers): lat, lon, tz, elev, 
    # (arrays): year, month, day, hour, minute, gh, dn, df, poa, wspd, wdir, tdry, twet, tdew, rhum, pres, snow, alb, aod
    
    return weatherdata
end
function set_ssc_data_from_dict(D,model,data)
    j = 0
    for (key, value) in D
        # if key == "solar_resource_file"
        #     print(typeof(value))
        #     continue
        # elseif typeof(value) == String
        #     @ccall hdl.ssc_data_set_string(data::Ptr{Cvoid},key::Cstring,D[key]::Cstring)::Cvoid
        #     j += 1
        if typeof(value) == String
            @ccall hdl.ssc_data_set_string(data::Ptr{Cvoid},key::Cstring,D[key]::Cstring)::Cvoid
            j += 1
        elseif typeof(D[key]) in [Int64,Float64]
            @ccall hdl.ssc_data_set_number(data::Ptr{Cvoid},key::Cstring,D[key]::Cdouble)::Cvoid
            j += 1
        elseif typeof(D[key]) == Vector{Any} || typeof(D[key]) == Vector{Float64} || typeof(D[key]) == Vector{Int64}
            nrows, ncols = length(D[key]), length(D[key][1])
            c_matrix = []
            for k in 1:nrows
                for l in 1:ncols
                    push!(c_matrix,D[key][k][l])
                end
            end
            if ncols == 1 && (nrows > 2 || model == "mst") || (model == "fresnel" && key == "ppa_price_input")
                # if key == "ppa_price_input"
                #     println("I am assigning ppa_prive_input as an array.")
                # end
                c_matrix = convert(Array{Float64},c_matrix)
                @ccall hdl.ssc_data_set_array(data::Ptr{Cvoid},key::Cstring,c_matrix::Ptr{Cdouble},length(D[key])::Cint)::Cvoid
                j += 1
            else
                # if key == "ppa_price_input"
                #     println("I am assigning ppa_prive_input as a matrix.")
                # end
                c_matrix = convert(Array{Float64},c_matrix)
                @ccall hdl.ssc_data_set_matrix(data::Ptr{Cvoid},key::Cstring,c_matrix::Ptr{Cdouble},Cint(nrows)::Cint,Cint(ncols)::Cint)::Cvoid
                j += 1
            end
        elseif typeof(D[key]) == Dict{Any,Any}
            table = @ccall hdl.ssc_data_create()::Ptr{Cvoid}  # data pointer
            set_ssc_data_from_dict(D[key],model,table)
            @ccall hdl.ssc_data_set_table(data::Ptr{Cvoid}, key::Cstring, table::Ptr{Cvoid})::Cvoid
            @ccall hdl.ssc_data_free(table::Ptr{Cvoid})::Cvoid
        else
            print("Could not assign variable " * key)
        end
        
    end
end

function run_ssc(cst_type::String,facility_id::Int,rated_power_electric::Float64)
    model = cst_type
    # relates internal names to specific models in SAM
    # TODO: Update load profile/dispatch to be equal to 8760 of rated power electric
    # TODO: Molten salt --> ideally we want optimize tower height and heliostat layout for each facility
    model_ssc = Dict(
        "trough" => "trough_physical",
        "fresnel" => "fresnel_physical",
        "mst" => "tcsmolten_salt"
    ) 
    # Assign all defaults
    defaults_file = joinpath(@__DIR__,"sam","defaults","defaults_" * model_ssc[model] * ".json") 
    defaults = JSON.parsefile(defaults_file)

    # Set inputs (only updating weather file location)
    defaults["solar_resource_file"] = joinpath(@__DIR__,"weatherfiles","weatherfile_" * string(facility_id) * ".csv") 

    # TODO: Add in design electric power assuming a SM = 3 based on given land
    

    R = Dict()
    error = ""
    
    if !(model in collect(keys(model_ssc)))
        error =  error * "Model is not available at this time. \n"
    else
        ### Setup SSC
        global hdl = nothing
        libfile = "ssc_2025.dll"
        global hdl = joinpath(@__DIR__, "sam", libfile)
        ssc_module = @ccall hdl.ssc_module_create(model_ssc[model]::Cstring)::Ptr{Cvoid}
        data = @ccall hdl.ssc_data_create()::Ptr{Cvoid}  # data pointer
        @ccall hdl.ssc_module_exec_set_print(1::Cint)::Cvoid # change to 1 to print outputs/errors (for debugging)

        ### Set inputs/defaults
        set_ssc_data_from_dict(defaults,model,data)

        ### Execute simulation
        test = @ccall hdl.ssc_module_exec(ssc_module::Ptr{Cvoid}, data::Ptr{Cvoid})::Cint
        #println(test)
        if test != 1
            println(string("Test = ", test,". SAM Simulation Failed."))
        else
            println(string("Test = ", test,". SAM Simulation Successful."))
        end
        ### Retrieve results
        
        outputs_dict = Dict(
            "trough" => ["P_out_net","total_land_area"], # Units: [MWe, acre]
            "fresnel" => ["P_out_net","total_land_area"], # Units: [MWe, acre]
            "mst" => ["P_out_net","total_land_area"] # Units: [MWe, acre]
        )
        outputs = outputs_dict[model]
        annual_energy = []
        total_area = 0.0
        for k in outputs
            if k in ["P_out_net"]
                len = 0
                len_ref = Ref(len)
                c_response = @ccall hdl.ssc_data_get_array(data::Ptr{Cvoid}, k::Cstring, len_ref::Ptr{Cvoid})::Ptr{Float64}
                
                for i in 1:8760
                    push!(annual_energy,unsafe_load(c_response,i))  # For array type outputs
                end
                println(string("Annual Electricity Production [MWhe]: ",round(Int,sum(annual_energy))))
                #annual_energy_string = joinpath(@__DIR__, "results", string(model), string("annual_energy_MW_",facility_id,".csv"))
                #df_annual_energy = DataFrame("Electricity Produced [MW]" => annual_energy)
                #CSV.write(annual_energy_string,df_annual_energy)
                #R[k] = sum(annual_energy)

            elseif k in ["total_land_area"]
                val = convert(Cdouble, 0.0)
                ref = Ref(val)
                @ccall hdl.ssc_data_get_number(data::Ptr{Cvoid}, k::Cstring, ref::Ptr{Cdouble})::Cvoid
                total_area = Float64(ref[])
                println(string("Total Land Area [acre]: ",round(Int,total_area)))
            end
            ### Can we get rated power (electric and thermal) as well as total area
            
        end
        df_summary = DataFrame("Annual Electricity, Net [MWh]" => sum(annual_energy),
                               "Electricity Produced [MW]" => annual_energy,
                               "Total Area [acre]" => total_area)
        df_summary_string = joinpath(@__DIR__, "results", string(model), string("results_",model,"_",facility_id,".csv"))
        CSV.write(df_summary_string,df_summary)
        ### Free SSC
        @ccall hdl.ssc_module_free(ssc_module::Ptr{Cvoid})::Cvoid   
        @ccall hdl.ssc_data_free(data::Ptr{Cvoid})::Cvoid

    end
    
    ### Check for errors
    if error == ""
        error = "No errors found."
    end
    R["error"] = error
    #return R
    return R
end

result = run_ssc("mst",100,200.0)
