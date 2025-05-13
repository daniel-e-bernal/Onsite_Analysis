"""
This file contains the script to run REopt on the HPC.

This is not meant to run the full optimized techo-economic capabilities of REopt.
We're only doing technical potential so the only real inputs are 
    A) space availability, (constraint)
    B) load profile, (hopefully scaled to the specific industrial sector)
        i. electrical,
        ii. thermal.
    C) lat/long.
"""

using REopt, DataFrames, XLSX, JSON, DeliminatedFiles, Xpress

#also potentially using PyCall to call on the datasets that are within the HPC storage such as NSRDBX and WINDToolkit
#using PyCall

# Function to safely extract values from JSON with default value if key is missing
function safe_get(data::Dict{String, Any}, keys::Vector{String}, default=0)
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

#read CSV file with coordinates
csv_land_file = #input location of file, supposed to contain the lat/long and space availability of all/representative sites 
coord_data = CSV.read( ,DataFrame)
lat = coord_data.Latitude 
long = coord_data.Longitude 

#read XLSX file with site information
xlsx_site_file = #input location of file, supposed to contain the lat/long and space availability of all/representative sites 


# make a new logic for selecting the wind size_class
function wind_size_class()
    """ First Constraint - Migration Paths
    the first priority in determing size class is within migration paths 
    if within migration paths then we need to cap the hub height to XX meters
    """
    # if else statement
    """ Second Constraint - Military Installations
    The next is determining if the site is within XX miles radius of a military installation.
    """
    #if else statement
    """ Third Constraint - Land availability
    This constraint is going to be determined by ["Site"]["land_acres"]
    """
    #if else statement
end

#Setup inputs for Solar run 
data_file = "solar_runs.json"
input_data = JSON.parsefile("Input Resources/$data_file")

# Electric Loads and Thermal Loads
"""
I will need to figure out about the electric and thermal load because it all depends on how
Blake and Indraneel output the different loads based on sector and climate zone.
"""

# site analysis to store results
sites_iter = eachindex() #will contain the number of sites, likely ~50,000
for i in sites_iter
    input_data_site = copy(input_data)

    # site specific inputs
    input_data_site["Site"]["latitude"] = lat[i]
    input_data_site["Site"]["longitude"] = long[i]
    input_data_site["ElectricLoad"]["loads_kw"] = #will input the loads here, likely a function to select type of load depending on sector and climate zone
    input_data_site["Wind"]["size_class"] = wind_size_class()


    # set up scenario and inputs
    s = Scenario(input_data_site)
    inputs = REoptInputs(s)

    append!(site_analysis, [(input_data_site)])
end
print("Completed run")

# Populate the DataFrame with the results produced and inputs
df = DataFrame(
    Site = sites,
    PV_size = [round(safe_get(site_analysis[i][2], ["PV", "size_kw"]), digits=0) for i in sites_iter],
    PV_year1_production = [round(safe_get(site_analysis[i][2], ["PV", "year_one_energy_produced_kwh"]), digits=0) for i in sites_iter],
    PV_annual_energy_production_avg = [round(safe_get(site_analysis[i][2], ["PV", "annual_energy_produced_kwh"]), digits=0) for i in sites_iter],
    PV_energy_lcoe = [round(safe_get(site_analysis[i][2], ["PV", "lcoe_per_kwh"]), digits=0) for i in sites_iter],
    PV_energy_exported = [round(safe_get(site_analysis[i][2], ["PV", "annual_energy_exported_kwh"]), digits=0) for i in sites_iter],
    PV_energy_curtailed = [sum(safe_get(site_analysis[i][2], ["PV", "electric_curtailed_series_kw"], 0)) for i in sites_iter],
    PV_energy_to_Battery_year1 = [sum(safe_get(site_analysis[i][2], ["PV", "electric_to_storage_series_kw"], 0)) for i in sites_iter],
    Grid_Electricity_Supplied_kWh_annual = [round(safe_get(site_analysis[i][2], ["ElectricUtility", "annual_energy_supplied_kwh"]), digits=0) for i in sites_iter],
    Total_Annual_Emissions_CO2 = [round(safe_get(site_analysis[i][2], ["Site", "annual_emissions_tonnes_CO2"]), digits=4) for i in sites_iter],
    ElecUtility_Annual_Emissions_CO2 = [round(safe_get(site_analysis[i][2], ["ElectricUtility", "annual_emissions_tonnes_CO2"]), digits=4) for i in sites_iter],
    BAU_Total_Annual_Emissions_CO2 = [round(safe_get(site_analysis[i][2], ["Site", "annual_emissions_tonnes_CO2_bau"]), digits=4) for i in sites_iter],
    LifeCycle_Emissions_CO2 = [round(safe_get(site_analysis[i][2], ["Site", "lifecycle_emissions_tonnes_CO2"]), digits=2) for i in sites_iter],
    BAU_LifeCycle_Emissions_CO2 = [round(safe_get(site_analysis[i][2], ["Site", "lifecycle_emissions_tonnes_CO2_bau"]), digits=2) for i in sites_iter],
    LifeCycle_Emission_Reduction_Fraction = [round(safe_get(site_analysis[i][2], ["Site", "lifecycle_emissions_reduction_CO2_fraction"]), digits=2) for i in sites_iter],
    LifeCycle_capex_costs_for_generation_techs = [round(safe_get(site_analysis[i][2], ["Financial", "lifecycle_generation_tech_capital_costs"]), digits=2) for i in sites_iter],
    LifeCycle_capex_costs_for_battery = [round(safe_get(site_analysis[i][2], ["Financial", "lifecycle_storage_capital_costs"]), digits=2) for i in sites_iter],
    Initial_upfront_capex_wo_incentives = [round(safe_get(site_analysis[i][2], ["Financial", "initial_capital_costs"]), digits=2) for i in sites_iter],
    Initial_upfront_capex_w_incentives = [round(safe_get(site_analysis[i][2], ["Financial", "initial_capital_costs_after_incentives"]), digits=2) for i in sites_iter],
    Yr1_energy_cost_after_tax = [round(safe_get(site_analysis[i][2], ["ElectricTariff", "year_one_energy_cost_before_tax"]), digits=2) for i in sites_iter],
    Yr1_demand_cost_after_tax = [round(safe_get(site_analysis[i][2], ["ElectricTariff", "year_one_demand_cost_before_tax"]), digits=2) for i in sites_iter],
    Yr1_total_energy_bill_before_tax = [round(safe_get(site_analysis[i][2], ["ElectricTariff", "year_one_bill_before_tax"]), digits=2) for i in sites_iter],
    Yr1_export_benefit_before_tax = [round(safe_get(site_analysis[i][2], ["ElectricTariff", "year_one_export_benefit_before_tax"]), digits=2) for i in sites_iter],
    Annual_renewable_electricity_kwh = [round(safe_get(site_analysis[i][2], ["Site", "annual_renewable_electricity_kwh"]), digits=2) for i in sites_iter],
    Annual_renewable_electricity_kwh_fraction = [round(safe_get(site_analysis[i][2], ["Site", "renewable_electricity_fraction"]), digits=2) for i in sites_iter],
    npv = [round(safe_get(site_analysis[i][2], ["Financial", "npv"]), digits=2) for i in sites_iter],
    lcc = [round(safe_get(site_analysis[i][2], ["Financial", "lcc"]), digits=2) for i in sites_iter]
    )
println(df)

write.("./results/validation_run_results.json", JSON.json(site_analysis))
println("Successfully printed results on JSON file")

# where to store the results
file_storage_location = "./results/validation_run_results.xlsx" # later to change

# Check if the Excel file already exists
if isfile(file_storage_location)
    # then open the excel file in read-write mode 'rw'
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
        # add new sheet
        XLSX.addsheet!(xf, sheet_name)
        # write dataframe to the new sheet
        XLSX.writetable!(xf[sheet_name], df)
    end
else
    # write dataframe to new excel file
    XLSX.writetable!(file_storage_location, df)
end
