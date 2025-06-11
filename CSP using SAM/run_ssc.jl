### Script for running any type of SSC module
using JSON
import HTTP
using DelimitedFiles
using DataFrames
using CSV
using Base.Iterators

function get_weatherdata(;lat::Float64, lon::Float64, debug::Bool, api_dv::String, facility_id::String)
    ### Call NSRDB
    # api_jgifford = "wKt35uq0aWoNHnzuwbcUxElPhVuo0K18YPSgZ9Ph"
    api_jgifford = "IOEFhaZoaOxkB2l9XvkaFswXgTsJxRqob1hBbPdv" #as of 9/13/2024 2pm
    if api_dv == "" #if it is empty, use the default api key
        api_dv = api_jgifford  #as of 9/13/2024 2pm
    end
    attributes_tmy = "ghi,dhi,dni,wind_speed,wind_direction,air_temperature,surface_pressure,dew_point"
    url = string("http://developer.nrel.gov/api/nsrdb/v2/solar/psm3-2-2-tmy-download.csv?api_key=",api_dv,
        "&wkt=POINT(",lon,"%20",lat,")&attributes=",attributes_tmy,
        "&names=tmy&utc=false&leap_day=true&interval=60&email=jeffrey.gifford@nrel.gov")
    # r = HTTP.request("GET", url)
    
    df = DataFrame(CSV.File(HTTP.get(url).body,delim=",",silencewarnings=true))
    
    ### Write csv file for checking (can comment out/delete when not debugging)
    #debug = false
    if debug
        weatherfile_path = joinpath(@__DIR__,"weatherfiles", "weatherfile_$(facility_id)_wdir.csv")
        CSV.write(weatherfile_path, df)
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
                c_matrix = convert(Array{Float64},c_matrix)
                @ccall hdl.ssc_data_set_array(data::Ptr{Cvoid},key::Cstring,c_matrix::Ptr{Cdouble},length(D[key])::Cint)::Cvoid
                j += 1
            else
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

function run_ssc(csp_type::String, facility_id::String, rated_power::Float64) 
    ### Function input definitions
    # csp_type (String) = type of csp technology
    # facility_id (String) = String assigned to that facility, will be used for referencing weather files, saving results, etc.
    # rated_power (Float) = Peak electrical demand of the facility [MW]
    
    model = csp_type
    # relates internal names to specific models in SAM
    model_ssc = Dict(
        "trough" => "trough_physical",
        "fresnel" => "fresnel_physical",
        "mst" => "tcsmolten_salt"
    ) 
    
    R = Dict()
    error = ""
    
    if !(model in collect(keys(model_ssc)))
        println("Model is not available at this time")
    else
        ### Setup SSC
        global hdl = nothing
        libfile = "ssc_2025.dll"
        global hdl = joinpath(@__DIR__, "sam", libfile)
        ssc_module = @ccall hdl.ssc_module_create(model_ssc[model]::Cstring)::Ptr{Cvoid}
        data = @ccall hdl.ssc_data_create()::Ptr{Cvoid}  # data pointer
        @ccall hdl.ssc_module_exec_set_print(0::Cint)::Cvoid # 1 to print outputs/errors (for debugging)

        ### Set defaults
        defaults_file = joinpath(@__DIR__,"sam","defaults","defaults_" * model_ssc[model] * ".json") 
        defaults = JSON.parsefile(defaults_file)

        # Set inputs (only updating weather file location)
        defaults["solar_resource_file"] = joinpath(@__DIR__, "weatherfiles", "weatherfile_$(facility_id)_wdir.csv") 

        sm = 3.0 # This is going to be a constant, external function sets rated power
        if model in ["trough"]
            defaults["P_ref"] = rated_power
            defaults["use_solar_mult_or_aperture_area"] = 0 # 0 = use solar multiple
            defaults["specified_solar_multiple"] = sm
            defaults["tshours"] = 12.0
        elseif model in ["fresnel"]
            defaults["P_ref"] = rated_power 
            defaults["use_solar_mult_or_aperture_area"] = 0 # 0 = use solar multiple
            defaults["specified_solar_multiple"] = sm
            defaults["tshours"] = 12.0
        else #model in ["mst"]
            defaults["P_ref"] = rated_power 
            defaults["field_model_type"] = 0 # 0 = optimize tower height and heliostat field
            defaults["solarm"] = sm
            defaults["tshours"] = 12.0
        end

        set_ssc_data_from_dict(defaults, model, data)

        ### Execute simulation
        test = @ccall hdl.ssc_module_exec(ssc_module::Ptr{Cvoid}, data::Ptr{Cvoid})::Cint
        #println(test)
        if test != 1
            println("SAM Simulation Failed.")
        else
            println("SAM Simulation Successful.")
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

            elseif k in ["total_land_area"]
                val = convert(Cdouble, 0.0)
                ref = Ref(val)
                @ccall hdl.ssc_data_get_number(data::Ptr{Cvoid}, k::Cstring, ref::Ptr{Cdouble})::Cvoid
                total_area = Float64(ref[])
                println(string("Total Land Area [acre]: ",round(Int,total_area)))
            end
        end
        ### ADD System size parameters that 
        results_summary = Dict("Annual Electricity, Net [MWh]" => sum(annual_energy),
                               "Electricity Produced [MW]" => annual_energy,
                               "Total Area [acre]" => total_area,
                               "Solar Multiple [-]" => sm,
                               "Rated Power [MW]" => rated_power,
                               "Capacity Factor [-]" => sum(annual_energy)/(8760*rated_power))
        
        ### Free SSC
        @ccall hdl.ssc_module_free(ssc_module::Ptr{Cvoid})::Cvoid   
        @ccall hdl.ssc_data_free(data::Ptr{Cvoid})::Cvoid

    end
    
    return results_summary 
end

# models = ["trough"] #,"fresnel","mst"]
# for m in models
#     println(m)
#     result = run_ssc(m,100,200.0)
# end

function store_ssc(result::Dict, csp_type::String, facility_id::String, option::String)
    # Determine the max length of the array values
    maxlen = maximum([isa(v, AbstractArray) ? length(v) : 1 for v in values(result)])

    # Normalize: expand scalars to [value, "", "", ...]
    normalized = Dict{String, Vector}()

    for (k, v) in result
        if isa(v, AbstractArray)
            normalized[k] = v
        else
            normalized[k] = [string(v); fill("", maxlen - 1)]
        end
    end

    # Convert to DataFrame and write CSV
    result_df = DataFrame(normalized)
    result_string = joinpath(@__DIR__, "results", csp_type, string("option ",option), "result_$(facility_id)_$(csp_type).csv") 
    CSV.write(result_string, result_df )
end


function run_ssc_options(csp_type::String, facility_id::String, option::String, peak_power::Float64, annual_demand::Float64, available_area::Float64; tol=0.01, max_iter=100, damping_factor=0.75)
    ### Function input definitions
    # csp_type (String) : type of csp technology
    # facility_id (Int) : Integer assigned to that facility, will be used for referencing weather files, saving results, etc.
    # option (String) : A = no export; B = annual load; C = maximum area
    # peak_power (Float) : Facility's peak electrical demand [MW]
    # annual_demand (Float) : Annual total electricity demanded at the facility [MWhe]
    # available_area (Float) : Total available land at facility [acre]
    # tol (Float) : tolerance of convergence
    # max_iter (Int) : maximum number of iterations to converge
    # damping_factor (float) : iteration term to prevent unstable oscillations
    if option in ["B"]
        rated_power = peak_power
        println("Beginning Option B.")
        for i in 1:max_iter
            result = run_ssc(csp_type, facility_id, rated_power)
            annual_generation = result["Annual Electricity, Net [MWh]"]
            #print(annual_generation)
            error_ratio = (annual_generation - annual_demand) / annual_demand
            println(string("Error ratio: ",string(round(error_ratio,digits=4))))
            if abs(error_ratio) <= tol
                println("Converged after $i iterations.")
                store_ssc(result, csp_type, facility_id, option)
                println("Final Rated Power: ", rated_power)
                println("Successfully completed Option B.")
                return result
            end

            # Adjust rated power proportionally
            rated_power *= 1 / (1 + damping_factor*error_ratio)
        end

        println("Did not converge after $max_iter iterations.")
        return result
    elseif option in ["C"]
        rated_power = peak_power
        println("Beginning Option C.")
        for i in 1:max_iter
            result = run_ssc(csp_type, facility_id, rated_power)
            total_area = result["Total Area [acre]"]
            #print(annual_generation)
            error_ratio = (total_area - available_area) / available_area
            println(string("Error ratio: ",string(round(error_ratio,digits=4))))
            if abs(error_ratio) <= tol
                println("Converged after $i iterations.")
                store_ssc(result, csp_type, facility_id, option)
                println("Successfully completed Option C.")
                return result
            end

            # Adjust rated power proportionally
            rated_power *= 1 / (1 + damping_factor*error_ratio)
        end

        println("Did not converge after $max_iter iterations.")
        return result
    end
    
end

#=
csp_type = "trough"
facility_id = 100
option = "C"
peak_power = 150.0
annual_demand = 100.0*8760
println(string("Annual Energy Demand [Mwh]: ", string(round(Int,annual_demand))))
available_area = 300.0
rated_power = run_ssc_options(csp_type, facility_id, option, peak_power, annual_demand, available_area)
print(rated_power)
=#
