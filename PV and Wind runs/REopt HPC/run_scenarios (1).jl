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
    
    using JSON, Xpress, REopt, JuMP, DelimitedFiles, Dates

    # File imports
    loads = readdlm("load_profiles_all_zones.csv",',')
    loads[1,1] = "1A"

    monthly_kwh = readdlm("usage_corrected.csv", ',')
    pv_prod_factors = readdir("pv_prod_factors")

    # results_superset = readdir("results")
    
    # Add common site information
    function run_site(i::Vector)

        lat = i[24]
        lon = i[25]
        location_number = i[2]
        urdb_rate = i[end-2]
        blended_rate = i[end-4]

        filename = string("pv10_bess_", location_number, ".json")

        pvprf_1 = string(location_number, ".json")
        pvprf_2 = string("bess_", location_number, ".json")

        request = Dict()
        request["Site"] = Dict()
        request["Site"]["latitude"] = lat
        request["Site"]["longitude"] = lon

        request["Financial"] = Dict()
        request["Financial"]["om_cost_escalation_rate_fraction"] = 0.023 #2.3% per NIST Office of Applied Economics
        request["Financial"]["elec_cost_escalation_rate_fraction"] = 0.019 # 1.9% REopt default
        request["Financial"]["offtaker_discount_rate_fraction"] = 0.047 # 4.7% per vz bond rate of return (cbonds)
        request["Financial"]["third_party_ownership"] = false
        request["Financial"]["offtaker_tax_rate_fraction"] = 0.21
        request["Financial"]["analysis_years"] = 25

        request["ElectricLoad"] = Dict()
        request["ElectricLoad"]["year"] = 2021 # Change?

        # zone is second element of tuple.
        zone = REopt.find_ashrae_zone_city(lat, lon; get_zone=true)[2]
        
        # load column by climate zone
        load_profile_column = findall(x -> x == String(zone), loads[1,:])[1] # returns a vector
        # normally this is a matrix, so index into it.
        loads_kw = loads[2:end,load_profile_column][:,1]
        site_idx = findall(x -> x == location_number, monthly_kwh)[1][1]
        # scale loads by month.
    
        v = []
        dr = DateTime(2021):Dates.Hour(1):DateTime(2021,12,31,23,45)
        
        for i in 1:12
            s = sum(loads_kw[findall(x -> month(x) == i, dr)])
            scaling_factor = monthly_kwh[site_idx, i+2] / s
            push!(v, scaling_factor.*loads_kw[findall(x -> month(x) == i, dr)])
        end
        request["ElectricLoad"]["loads_kw"] = vcat(v...)

        request["ElectricTariff"] = Dict()
        if urdb_rate == ""
            request["ElectricTariff"]["blended_annual_energy_rate"] = blended_rate # $/kWh rate
        else
            request["ElectricTariff"]["urdb_response"] = JSON.parsefile(joinpath("urdb_rates", string(urdb_rate, ".json")))
        end

        ## PV only scenarios
        request["PV"] = Dict()
        request["PV"]["array_type"] = 1 # fixed rack rooftop
        request["PV"]["module_type"] = 1 # fixed rack rooftop
        request["PV"]["inv_eff"] = 0.95 # 5% efficiency hit due to charge controller
        request["PV"]["tilt"] = 5 # not default 10 for rooftop PV to prevent obstruction
        request["PV"]["azimuth"] = 180
        request["PV"]["dc_ac_ratio"] = 1.0 # no DC -> AC conversion losses
        request["PV"]["installed_cost_per_kw"] = 5000 # includes all costs per VZ
        request["PV"]["om_cost_per_kw"] = 17 # includes all costs per VZ
        request["PV"]["macrs_option_years"] = 5
        request["PV"]["macrs_bonus_fraction"] = 0.8
        request["PV"]["macrs_itc_reduction"] = 0.5
        request["PV"]["federal_itc_fraction"] = 0.3
        request["PV"]["location"] = "roof"
        try
            if pvprf_1 in pv_prod_factors
                request["PV"]["production_factor_series"] = readdlm(joinpath("pv_prod_factors", pvprf_1))[:,1]
            elseif pvprf_2 in pv_prod_factors
                request["PV"]["production_factor_series"] = readdlm(joinpath("pv_prod_factors", pvprf_2))[:,1]
            else
                nothing
            end
        catch e
            nothing
        end

        ## PV BESS scenarios
        request["ElectricStorage"] = Dict()
        request["ElectricStorage"]["charge_efficiency"] = 1.0 # PV to BESS losses accounted for in PV
        request["ElectricStorage"]["discharge_efficiency"] = 0.95 # Account for internal BESS losses
        request["ElectricStorage"]["grid_charge_efficiency"] = 0.969 # 96.9% efficiency of rectifier

        request["ElectricStorage"]["soc_min_fraction"] = 0.2 # Min SOC %
        request["ElectricStorage"]["can_grid_charge"] = true
        request["ElectricStorage"]["inverter_replacement_year"] = 12
        request["ElectricStorage"]["battery_replacement_year"] = 12

        request["ElectricStorage"]["installed_cost_per_kw"] = 113.0 # rectifier costs
        request["ElectricStorage"]["installed_cost_per_kwh"] = 500.0
        request["ElectricStorage"]["replace_cost_per_kw"] = 113.0 # rectifier costs
        request["ElectricStorage"]["replace_cost_per_kwh"] = 400.0
        
        request["ElectricStorage"]["macrs_option_years"] = 5
        request["ElectricStorage"]["macrs_bonus_fraction"] = 0.0
        request["ElectricStorage"]["macrs_itc_reduction"] = 0.0
        request["ElectricStorage"]["total_itc_fraction"] = 0.3

        # PV 5 kW and optimal BESS
        request["PV"]["min_kw"] = 10.0
        request["PV"]["max_kw"] = 10.0
        
        m1 = Model(optimizer_with_attributes(Xpress.Optimizer, "OUTPUTLOG" => 0))
        m2 = Model(optimizer_with_attributes(Xpress.Optimizer, "OUTPUTLOG" => 0))
        results = run_reopt([m1,m2], request)
        empty!(m1)
        empty!(m2)
        results["Request"] = request

        # outage simulator starts here.
        time_steps = 8760

        # Generator only
        bess_kwh = 0.0
        bess_kw = 0.0
        initsoc = zeros(time_steps)

        microgrid_only = false
        pv_kw_ac_hourly = zeros(time_steps)

        o1 = simulate_outages(;
            batt_kwh = bess_kwh, 
            batt_kw = bess_kw, 
            pv_kw_ac_hourly = pv_kw_ac_hourly,
            init_soc = initsoc, 
            critical_loads_kw = request["ElectricLoad"]["loads_kw"],
            wind_kw_ac_hourly = zeros(time_steps),
            batt_roundtrip_efficiency = 0.95,
            diesel_kw = 60.0, # generator size and fuel tank size per Verizon
            fuel_available = 210.0,
            b = 0.0, # y intercept of fuel curve slope
            m = 0.076, # slope of fuel curve slope
            diesel_min_turndown = 0.0
        )
        
        results["Gen_only_osim"] = o1
        
        # DER with generators.
        if haskey(results, "ElectricStorage")
            bess_kwh = results["ElectricStorage"]["size_kwh"]
            bess_kw = results["ElectricStorage"]["size_kw"]
            initsoc = results["ElectricStorage"]["soc_series_fraction"]
        else
            bess_kwh = 0.0 
            bess_kw = 0.0
            initsoc = zeros(time_steps)
        end;

        microgrid_only = false
        pv_kw_ac_hourly = (
            get(results["PV"], "electric_to_storage_series_kw", zeros(time_steps))
            + get(results["PV"], "electric_curtailed_series_kw", zeros(time_steps))
            + get(results["PV"], "electric_to_load_series_kw", zeros(time_steps))
            + get(results["PV"], "electric_to_grid_series_kw", zeros(time_steps))
        );

        o2 = simulate_outages(;
            batt_kwh = bess_kwh, 
            batt_kw = bess_kw, 
            pv_kw_ac_hourly = pv_kw_ac_hourly,
            init_soc = initsoc, 
            critical_loads_kw = request["ElectricLoad"]["loads_kw"],
            wind_kw_ac_hourly = zeros(time_steps),
            batt_roundtrip_efficiency = 0.95,
            diesel_kw = 60,
            fuel_available = 210.0,
            b = 0.0,
            m = 0.076, 
            diesel_min_turndown = 0.0
        )
        
        results["Gen_DERs_osim"] = o2

        write(joinpath("results", filename), JSON.json(results))
        sleep(1.0)
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
scenarios = readdlm("inputs_with_rates.csv", ',');

@info size(scenarios)

sta_idx = 2
end_idx = size(scenarios)[1]

@time pmap(25000:end_idx) do i
    try
        # Pass a vector of values for each site.
        run_site(scenarios[i,:])
    catch e
    	@info e
        @warn "Error for " scenarios[i, 2]
    end
end

# remove the workers
for i in workers()
    rmprocs(i)
end