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
        n_procs = num_workers
    else
        num_workers = []
        for (i, hostname) in enumerate(hostnames)
            if i == 1
                push!(num_workers, (hostname, n_procs_per_node-1))
            else
                push!(num_workers, (hostname, n_procs_per_node))
            end
        end
        n_procs = Int(n_procs_per_node*n_nodes)
    end
elseif Sys.isapple()
    n_procs = parse(Int, chomp(read(`sysctl -n hw.physicalcpu`, String)))
    num_workers = n_procs - 1
else
    n_procs = 11
    num_workers = n_procs - 1
end

Distributed.addprocs(num_workers)

# Env instantiation
@everywhere begin
    using Pkg
    Pkg.activate(@__DIR__)
end

# Function declarations, module imports
@everywhere begin
    
    using JSON, DelimitedFiles, CSV, DataFrames

    ## Read the files
    results_files = readdir("v2_results"; join=true);
    
    function extract_results(staidx::Int,endidx::Int)
        
        @info staidx,",", endidx

        df = DataFrame(
            LocationNumber=[],
            Scenario=[],
            Status=[],
            Lat=[],
            Lon=[],
            solver_seconds=[],
            URDB_Rate=[],
            Blended_Rate=[],
            PV_kW=[],
            BESS_kW=[],
            BESS_kwh=[],
            LCC_bau=[],
            LCC=[],
            NPV=[],
            SPP=[],
            IRR=[],
            Breakeven_cost_co2pertonne = [],
            Year1_demand_ch_bau=[],
            Year1_demand_ch=[],
            Year1_energy_ch_bau=[],
            Year1_energy_ch=[],
            Total_demand_ch_bau=[],
            Total_demand_ch=[],
            Total_energy_ch_bau=[],
            Total_energy_ch=[],
            Year_one_co2_em_bau=[],
            Year_one_co2_em = [],
            LCC_CO2_em_bau = [],
            LCC_CO2_em = [],
            Pct_Renewable_Electricity_bau = [],
            Pct_Renewable_Electricity = [],
            Energy_supplied_kwh_bau=[],
            Energy_supplied_kwh=[],
            PV_to_site_kwh=[],
            Curtailed_solar_kwh=[],
            Grid_to_bess_kwh=[],
            Bess_to_site_kwh=[]
        )

        for idx in staidx:endidx

            sleep(0.5)

            f = JSON.parsefile(results_files[idx])

            urdbrate = "na"
            blndrate = 0.0
            try
                urdbrate = f["Request"]["ElectricTariff"]["urdb_response"]["label"]
            catch
                blndrate = f["Request"]["ElectricTariff"]["blended_annual_energy_rate"]
            end

            lat = f["Request"]["Site"]["latitude"]
            lon = f["Request"]["Site"]["longitude"]

            results_keys = ["PV_10kW_BESS", "BESS_only", "PV_10kW", "PV_5kW_BESS", "PV_5kW"]

            for k in results_keys

                vec = []

                try

                    if haskey(f[k], "PV")
                        pv_size_kw =  f[k]["PV"]["size_kw"]

                        pv_to_load = sum(f[k]["PV"]["electric_to_load_series_kw"])
                        pv_to_curt = sum(f[k]["PV"]["electric_curtailed_series_kw"])
                    else
                        pv_size_kw = 0.0
                        pv_to_load = 0.0
                        pv_to_curt = 0.0
                    end

                    if haskey(f[k], "ElectricStorage")
                        bess_size_kw =  f[k]["ElectricStorage"]["size_kw"]
                        bess_size_kwh = f[k]["ElectricStorage"]["size_kwh"]

                        if bess_size_kwh > 0
                            grid_to_bess = sum(f[k]["ElectricUtility"]["electric_to_storage_series_kw"])
                            stor_to_load = sum(f[k]["ElectricStorage"]["storage_to_load_series_kw"])
                        else
                            grid_to_bess = 0.0
                            stor_to_load = 0.0
                        end
                    else
                        bess_size_kw = 0.0
                        bess_size_kwh = 0.0

                        grid_to_bess = 0.0
                        stor_to_load = 0.0
                    end

                    Breakeven_cost_co2pertonne = nothing
                    try
                        Breakeven_cost_co2pertonne = f[k]["Financial"]["breakeven_cost_of_emissions_reduction_per_tonne_CO2"]
                    catch
                        Breakeven_cost_co2pertonne = 0.0
                    end

                    if haskey(f[k], "ElectricUtility")
                        lcc_co2_em_bau = f[k]["ElectricUtility"]["lifecycle_emissions_tonnes_CO2_bau"]
                        lcc_co2_em = f[k]["ElectricUtility"]["lifecycle_emissions_tonnes_CO2"]

                        co2_em_annual_bau = f[k]["ElectricUtility"]["annual_emissions_tonnes_CO2_bau"]
                        co2_em_annual = f[k]["ElectricUtility"]["annual_emissions_tonnes_CO2"]

                        grid_to_site_bau = f[k]["ElectricUtility"]["annual_energy_supplied_kwh_bau"]
                        grid_to_site = f[k]["ElectricUtility"]["annual_energy_supplied_kwh"]
                    else
                        lcc_co2_em_bau = 0.0
                        lcc_co2_em = 0.0

                        co2_em_annual_bau = 0.0
                        co2_em_annual = 0.0

                        grid_to_site_bau = 0.0
                        grid_to_site = 0.0
                    end

                    if haskey(f[k], "Site")
                        pctre_elec_bau = f[k]["Site"]["renewable_electricity_fraction_bau"]
                        pctre_elec = f[k]["Site"]["renewable_electricity_fraction"]
                        pctre_ener_bau = f[k]["Site"]["total_renewable_energy_fraction_bau"]
                        pctre_ener = f[k]["Site"]["total_renewable_energy_fraction"]
                    else
                        pctre_elec_bau = 0.0
                        pctre_elec = 0.0
                        pctre_ener_bau = 0.0
                        pctre_ener = 0.0
                    end
                    
                    push!(vec, [
                        results_files[idx],
                        k,
                        f[k]["status"],
                        lat,
                        lon,
                        f[k]["solver_seconds"],
                        urdbrate,
                        blndrate,
                        pv_size_kw,
                        bess_size_kw,
                        bess_size_kwh,
                        f[k]["Financial"]["lcc_bau"],
                        f[k]["Financial"]["lcc"],
                        f[k]["Financial"]["npv"],
                        f[k]["Financial"]["simple_payback_years"],
                        f[k]["Financial"]["internal_rate_of_return"],
                        Breakeven_cost_co2pertonne,
                        f[k]["ElectricTariff"]["year_one_demand_cost_before_tax_bau"],
                        f[k]["ElectricTariff"]["year_one_demand_cost_before_tax"],
                        f[k]["ElectricTariff"]["year_one_energy_cost_before_tax_bau"],
                        f[k]["ElectricTariff"]["year_one_energy_cost_before_tax"],
                        f[k]["ElectricTariff"]["lifecycle_demand_cost_after_tax_bau"],
                        f[k]["ElectricTariff"]["lifecycle_demand_cost_after_tax"],
                        f[k]["ElectricTariff"]["lifecycle_energy_cost_after_tax_bau"],
                        f[k]["ElectricTariff"]["lifecycle_energy_cost_after_tax"],
                        co2_em_annual_bau,
                        co2_em_annual,
                        lcc_co2_em_bau,
                        lcc_co2_em,
                        pctre_elec_bau,
                        pctre_elec,
                        grid_to_site_bau,
                        grid_to_site,
                        pv_to_load,
                        pv_to_curt,
                        grid_to_bess,
                        stor_to_load
                    ])
                catch e
                    push!(vec, vcat(
                        [
                            results_files[idx],
                            k,
                            "err"
                        ],
                        repeat([0], 34)
                        )
                    )
                    @info results_files[idx]
                end
                push!(df, vec...)
            end
        end
        # LN = replace(results_files[idx], "v2_results/"=>"", ".json"=>"")
        fname = string("results_", getpid(), ".csv")
        CSV.write(joinpath("/lustre/eaglefs/projects/bulkreopt/v2_verizon/extracted_results", fname), df)
        # writedlm(joinpath("/lustre/eaglefs/projects/bulkreopt/v2_verizon/extracted_results", fname), vec, ',')
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
results_files = length(readdir("v2_results"; join=true));

# create range
splits = [round(Int, s) for s in range(0, 10000; length=n_procs)]

println("\n*********************************************\n")

pmap(1:length(splits)-1) do s
    try
        extract_results(splits[s]+1,splits[s+1])
    catch e
        @info e
    end
end

# remove the workers
for i in workers()
    rmprocs(i)
end

