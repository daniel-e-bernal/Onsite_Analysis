# From https://researchcomputing.princeton.edu/support/knowledge-base/julia#distributed
# https://github.com/Arpeggeo/julia-distributed-computing

using Distributed

if Sys.islinux()
    # Create worker information
    hostnames = split(chomp(read(`scontrol show hostnames`, String)), "\n")
    n_procs_per_node = parse(Int, chomp(read(`nproc`, String)))
    n_nodes = length(hostnames)
    if n_nodes == 1
        num_workers = n_procs_per_node - 1
    else
        num_workers = []
        for (i, hostname) in enumerate(hostnames)
            if i == 1
                push!(num_workers, (hostname, n_procs_per_node-1))
            else
                push!(num_workers, (hostname, n_procs_per_node))
            end
        end
    end
elseif Sys.isapple()
    n_procs = parse(Int, chomp(read(`sysctl -n hw.physicalcpu`, String)))
    num_workers = n_procs - 1
else
    n_procs = 4
    num_workers = n_procs - 1
end

Distributed.addprocs(num_workers; topology=:master_worker)

# Env instantiation
@everywhere begin
    using Pkg
    Pkg.activate(@__DIR__)
end

# Function declarations, module imports
@everywhere begin
    
    using JSON, REopt, JuMP, DelimitedFiles, Dates
    
    # Add common site information
    function run_simulator(d::Dict, fname::String)

        d["PV_5kW"]["osim_genonly"] = run_outage_simulator(d["Request"]["ElectricLoad"]["loads_kw"])
        d["PV_5kW"]["osim_dersgen"] = run_outage_simulator(d["PV_5kW"])

        d["PV_10kW"]["osim_genonly"] = d["PV_5kW"]["osim_genonly"]
        d["PV_10kW"]["osim_dersgen"] = run_outage_simulator(d["PV_10kW"])

        d["BESS_only"]["osim_genonly"] = d["PV_5kW"]["osim_genonly"]
        d["BESS_only"]["osim_dersgen"] = run_outage_simulator(d["BESS_only"])

        d["PV_5kW_BESS"]["osim_genonly"] = d["PV_5kW"]["osim_genonly"]
        d["PV_5kW_BESS"]["osim_dersgen"] = run_outage_simulator(d["PV_5kW_BESS"])

        d["PV_10kW_BESS"]["osim_genonly"] = d["PV_5kW"]["osim_genonly"]
        d["PV_10kW_BESS"]["osim_dersgen"] = run_outage_simulator(d["PV_10kW_BESS"])

        write(joinpath("/lustre/eaglefs/projects/bulkreopt/v2_verizon/results_with_osim", string(fname, "")), JSON.json(d))

        GC.gc()
    end

    function run_outage_simulator(v::Vector)
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
            critical_loads_kw = v,
            wind_kw_ac_hourly = zeros(time_steps),
            batt_roundtrip_efficiency = 0.95,
            diesel_kw = 60.0, # generator size and fuel tank size per Verizon
            fuel_available = 210.0,
            b = 0.0, # y intercept of fuel curve slope
            m = 0.076, # slope of fuel curve slope
            diesel_min_turndown = 0.0
        )
        
        return o1
    end
    
    function run_outage_simulator(results::Dict)
        # outage simulator starts here.
        time_steps = 8760

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
            critical_loads_kw = results["ElectricUtility"]["electric_to_load_series_kw_bau"],
            wind_kw_ac_hourly = zeros(time_steps),
            batt_roundtrip_efficiency = 0.95,
            diesel_kw = 60,
            fuel_available = 210.0,
            b = 0.0,
            m = 0.076, 
            diesel_min_turndown = 0.0
        )
        
        return o2
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
executed_files = readdir("v2_results")

"""
MODIFY start and end index limits here:

    test: 2:20

    node1: 2:20,000
    node2: 20,001:40,000
    node3: 40,001:size(scenarios)[1]
"""

@info size(executed_files)

sta_idx = 1
end_idx = length(executed_files)

@time pmap(sta_idx:end_idx) do i
    try
        j = JSON.parsefile(joinpath("v2_results",executed_files[i]))
        if !haskey(j["PV_10kW_BESS"], "osim_genonly")
            run_simulator(j, executed_files[i])
            @info "executed scenarios for " executed_files[i]
        else
            write(joinpath("/lustre/eaglefs/projects/bulkreopt/v2_verizon/results_with_osim", executed_files[i]), JSON.json(j))
        end
    catch e
    	@info e
        @warn "Error for " executed_files[i,2]
    end
end

# remove the workers
for i in workers()
    rmprocs(i)
end
