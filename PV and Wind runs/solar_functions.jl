#setting up PV functions
using DataFrames
using Statistics

#=create a function that determines the sizing option for the PV system:
    Option A: Size the PV system to a baseload. The user will input a value (kW) for the baseload and it'll be sized to that value IF 
        there is enough space (land acres + roof space) to accommodate the system. If there isn't enough space, then the system will be sized to the maximum
        that can fit on the land acres + roof space.
    Option B: Size the PV system to offset the entire load. The user will input a value (kW) for the load and it'll be sized to that value IF 
        there is enough space (land acres + roof space) to accommodate the system. If there isn't enough space, then the system will be sized to the maximum
        that can fit on the land acres + roof space.
    Option C: Size the PV system to max out the roof + land space available.
    =#
function pv_size_option(;
    sizing_option::String, # "A", "B", or "C"
    roof_prod_series::Vector, #production factor series for fixed roof-mounted PV system
    ground_fixed_prod_series::Vector, #production factor series for fixed ground-mounted PV system
    ground_axis_prod_series::Vector, #production factor series for ground-mounted PV system with one-axis tracking
    baseload::Real, #kW
    load_kw::Vector, #series of kw
    land_avail::Real, #land acres available
    roof_avail::Real, #roof space available in ft2
    data_dict::Dict{String, Any} #dictionary with data for the site
)
    #check if the sizing option that was entered is valid 
    if !(sizing_option in ["A", "B", "C"])
        error("Invalid sizing option. Must be 'A', 'B', or 'C'.")
    end

    #check if the baseload and load_kwh are non-negative
    if baseload < 0 || land_avail < 0 || roof_avail < 0
        error("Baseload, land acres, and roof space must be non-negative.")
    end

    # setting up local variables
    roof_fixed_prod_factor_series = roof_prod_series
    ground_fixed_prod_factor_series = ground_fixed_prod_series
    ground_axis_prod_factor_series = ground_axis_prod_series
    load_vector = load_kw
    land_acres = land_avail
    roof_space = roof_avail
    one_axis_deployment = false
    input_data_site = data_dict

    #calculating capacity factors 
    PV_roof_capacity_factor = sum(roof_fixed_prod_factor_series) / 8760 #actual capacity factor
    PV_fixed_ground_cap_factor = sum(ground_fixed_prod_factor_series) / 8760 #actual capacity factor for fixed rack ground PV
    #if axis prod factor series is 0 then cap factor will be 0, otherwise, calculate actual capacity factor
    PV_axis_ground_cap_factor = ground_axis_prod_factor_series == 0 ? 0 : sum(ground_axis_prod_factor_series) / 8760
    
    #calculate max roof PV size that would be necessary to satisfy site's load
    PV_roof_max_size_based_on_load = sum(load_vector) / (PV_roof_capacity_factor * 8760)
    #calcualte max roof PV size possible based on roof space
    PV_roof_max_size_based_n_space = roof_space * (0.017) #power density is kw per sqft 

    #calculate max ground PV size that would be necessary to satisfy site's load
    PV_fixed_ground_max_size_based_on_load = sum(load_vector) / (PV_fixed_ground_cap_factor * 8760)
    PV_axis_ground_max_size_based_on_load = PV_axis_ground_cap_factor == 0 ? 0 : (sum(load_vector) / (PV_axis_ground_cap_factor * 8760))

    #calculate max ground PV sizes possible based on ground space
    PV_fixed_ground_max_size_based_on_space = land_acres * (1/0.0031) #power density is acres per kW 
    PV_axis_ground_max_size_based_on_space = land_acres * (1/0.0042) #power density is acres per kW
    PV_ground_max_size_based_on_space = PV_fixed_ground_max_size_based_on_space

    #If option = A - a set baseload (kW)
    if sizing_option == "A"
        #check if the baseload is less than or equal to the max roof size based on space 
        if baseload <= PV_roof_max_size_based_n_space
            println("Baseload $baseload (kW) for $match_id is less than or equal to PV roof max size based on space $PV_roof_max_size_based_n_space.")
            input_data_site["PV"][2]["min_kw"] = baseload * 0.999 #set the minimum kW for the roof PV
            input_data_site["PV"][2]["max_kw"] = baseload #set the maximum kW for the roof PV
            input_data_site["PV"][1]["max_kw"] = 0 #set the maximum kW for the ground PV to 0 
            return baseload, 0, one_axis_deployment, input_data_site
        elseif baseload > PV_roof_max_size_based_n_space
            println("Baseload $baseload (kW) for $match_id is greater than PV roof max size based on space $PV_roof_max_size_based_n_space.")
            input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_n_space * 0.999 #set the minimum kW for the roof PV
            input_data_site["PV"][2]["max_kw"] = PV_roof_max_size_based_n_space #set the maximum kW for the roof PV
            remaining_kw = baseload - PV_roof_max_size_based_n_space # to determine remaining kw 
            if remaining_kw < 1000 #if less than than 1 MW then stick to fixed ground mount 
                ground_PV_size = PV_fixed_ground_max_size_based_on_space < remaining_kw ? PV_fixed_ground_max_size_based_on_space : remaining_kw
                input_data_site["PV"][1]["min_kw"] = ground_PV_size * 0.999 #set the minimum kW for the ground PV
                input_data_site["PV"][1]["max_kw"] = ground_PV_size #set the maximum kW for the ground PV
            elseif remaining_kw >= 1000 && PV_axis_ground_max_size_based_on_space >= 1000 #if at or greater than 1 MW AND enough space for axis (>= 4.2 acres) then go to axis ground mount
                println("Switching to ground 1-axis PV for $match_id")
                ground_PV_size = PV_axis_ground_max_size_based_on_space >= remaining_kw ? remaining_kw : PV_axis_ground_max_size_based_on_space
                input_data_site["PV"][1]["min_kw"] = ground_PV_size * 0.999 #set the minimum kW for the ground PV
                input_data_site["PV"][1]["max_kw"] = ground_PV_size #set the maximum kW for the ground PV
                one_axis_deployment = true
            end
            return PV_roof_max_size_based_n_space, ground_PV_size, one_axis_deployment, input_data_site
        end
    #If option = B - attempt to offset entire load 
    elseif sizing_option == "B"
        #now re-set the sizes for PV on the roof and ground, with roof as the priority
        if PV_roof_max_size_based_on_load <= PV_roof_max_size_based_n_space
            println("PV roof max size based on load for $match_id is less than or equal to PV max size based on space.")
            input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_on_load * 0.999
            roof_PV_size = PV_roof_max_size_based_on_load
            input_data_site["PV"][2]["max_kw"] = roof_PV_size
            #since the roof size for the total load was calculated above and it is smaller than the total based on roof size, then we do not have any more PV to deploy and the ground PV is 0
            ground_PV_size = 0 
            input_data_site["PV"][1]["max_kw"] = ground_PV_size
        elseif PV_roof_max_size_based_on_load > PV_roof_max_size_based_n_space
            println("PV roof max size based on load for $match_id is greater than PV max size based on space.")
            #fix the minimum kW for the roof 
            input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_n_space * 0.999
            #identify the max kW for ground remainder after knowing the total load 
            roof_PV_size = PV_roof_max_size_based_n_space
            println("The PV roof size is: ", roof_PV_size)
            input_data_site["PV"][2]["max_kw"] = roof_PV_size
            
            #get remaining kWh that need to be supplied by subtracting total load - load produced from roof PV
            annual_production_from_roof_pv = sum(roof_PV_size * roof_fixed_prod_factor_series)
            production_defecit = sum(load_vector) - annual_production_from_roof_pv

            #calculate ground size using fixed rack production factor
            ground_PV_size = production_defecit / (PV_fixed_ground_cap_factor * 8760)
            println("The FIRST pass to assign ground PV size for $match_id is: ", ground_PV_size)
            ground_PV_size = ground_PV_size > PV_fixed_ground_max_size_based_on_space ? PV_fixed_ground_max_size_based_on_space : ground_PV_size
            println("The SECOND pass to assign ground PV size for $match_id is: ", ground_PV_size)
            #inputs1.max_sizes["ground_mount"] = PV_fixed_ground_max_size_based_on_space
            #sleep(10)

            #check if the remainder ground size is over 1354 kW, if so, then we need to simulate 1-axis which would equal 1000 kW
            if ground_PV_size > 1354 && PV_axis_ground_cap_factor != 0
                println("Switching to ground 1-axis PV")
                ground_PV_size = production_defecit / (PV_axis_ground_cap_factor * 8760)
                println("The First pass to assign ground axis PV size for $match_id is: ", ground_PV_size)
                ground_PV_size = ground_PV_size > PV_axis_ground_max_size_based_on_space ? PV_axis_ground_max_size_based_on_space : ground_PV_size
                println("The Second pass to assign ground axis PV size for $match_id is: ", ground_PV_size)
                #if ground_PV_size > PV_axis_ground_max_size_based_on_space
                    #ground_PV_size = PV_axis_ground_max_size_based_on_space
                #end
                #inputs1.max_sizes["ground_mount"] = PV_axis_ground_max_size_based_on_space
                PV_ground_max_size_based_on_space = PV_axis_ground_max_size_based_on_space
                println("The PV ground max size based on space for $match_id is (acres): ", PV_axis_ground_max_size_based_on_space)
                one_axis_deployment = true
                #sleep(10)
            end
        end
        return roof_PV_size, ground_PV_size, one_axis_deployment, input_data_site
    #If option = C - max out the roof + land space available
    elseif sizing_option == "C"
        #decide on axis or fixed for ground mount 
        if land_acres < 4.2 #if land acres is less than 4.2 then use fixed ground mount
            input_data_site["PV"][1]["max_kw"] = PV_fixed_ground_max_size_based_on_space
            input_data_site["PV"][1]["min_kw"] = PV_fixed_ground_max_size_based_on_space * 0.999
            input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_n_space * 0.999
            input_data_site["PV"][2]["max_kw"] = PV_roof_max_size_based_n_space
            return PV_roof_max_size_based_n_space, PV_fixed_ground_max_size_based_on_space, one_axis_deployment, input_data_site
        else #if land acres is greater than 4.2 then use axis ground mount
            input_data_site["PV"][1]["max_kw"] = PV_axis_ground_max_size_based_on_space
            input_data_site["PV"][1]["min_kw"] = PV_axis_ground_max_size_based_on_space * 0.999
            input_data_site["PV"][2]["min_kw"] = PV_roof_max_size_based_n_space * 0.999
            input_data_site["PV"][2]["max_kw"] = PV_roof_max_size_based_n_space
            one_axis_deployment = true
            return PV_roof_max_size_based_n_space, PV_axis_ground_max_size_based_on_space, one_axis_deployment, input_data_site
        end
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

#select production factor series based on match_id
function select_prod_factor(; 
    match_id::Any,
    roof_prod_factor_folder::String = pv_roof_prod_factors,
    ground_prod_fixed_factor_folder::String = pv_ground_fixed_prod_factors,
    ground_prod_axis_factor_folder::String = pv_ground_axis_prod_factors)
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

    get_roof_file = joinpath(roof_prod_factor_folder, "$match_id.csv")
    get_ground_fixed_file = joinpath(ground_prod_fixed_factor_folder, "$match_id.csv")
    get_ground_axis_file = joinpath(ground_prod_axis_factor_folder, "$match_id.csv")

    # Check if either file exists
    roof_exists = isfile(get_roof_file)
    ground_fixed_exists = isfile(get_ground_fixed_file)
    ground_axis_exists = isfile(get_ground_axis_file)

    if !roof_exists
        error("No production factor file found for match_id: $match_id in roof_prod_factor_folder.")
        @info "Downloading now...."
    end
    if !ground_fixed_exists && !ground_axis_exists
        error("No production factor file found for match_id: $match_id in ground_prod_fixed_factor_folder OR ground_prod_axis_factor_folder")
    end

    #read csv file and convert it to array, each file only has 1 column no header
    try
        data_ground_fixed = readdlm(get_ground_fixed_file, ',', Float64)  # Read only the first column
        data_ground_fixed = vec(data_ground_fixed)
        data_ground_axis = 0
        if ground_axis_exists == true
            data_ground_axis = readdlm(get_ground_axis_file, ',', Float64)  # Read only the first column
            data_ground_axis = vec(data_ground_axis)
        end
        #data_ground = collect(data_ground)  # Convert DataFrame column to an array
        data_roof = readdlm(get_roof_file, ',', Float64)
        data_roof = vec(data_roof)
        #data_roof = CSV.read(get_roof_file, DataFrame)[:, 1]  # Read only the first column
        #data_roof = collect(data_roof)  # Convert DataFrame column to an array
        return data_roof, data_ground_fixed, data_ground_axis
    catch e
        error("Error reading file $get_ground_fixed_file, $get_ground_axis_file and $get_roof_file: $e")
    end
end