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

iters = length(data[!, :latitude])
#iters = 1
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
    production_aggregate = sum(analysis_runs[i][2]["PV"][1]["electric_to_load_series_kw"] + analysis_runs[i][2]["PV"][2]["electric_to_load_series_kw"])
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
    production = analysis_runs[i][2]["PV"][1]["electric_to_load_series_kw"] + analysis_runs[i][2]["PV"][2]["electric_to_load_series_kw"]
    #println("The type of object is ", typeof(production))
    production_summer_week = production[summer_week_start:summer_week_end]
    #println("The type of object is ", typeof(production_summer_week))

    plot(
        x_summer,
        production_summer_week,
        label="Production",
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
    production = analysis_runs[i][2]["PV"][1]["electric_to_load_series_kw"] + analysis_runs[i][2]["PV"][2]["electric_to_load_series_kw"]
    #println("The type of object is ", typeof(production))
    production_winter_week = production[winter_week_start:winter_week_end]
    #println("The type of object is ", typeof(production_summer_week))

    plot(
        x_winter,
        production_winter_week,
        label="Production",
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