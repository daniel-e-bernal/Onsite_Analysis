"""
The code below is just psuedo code for deciding on wind and solar for the Analysis runs.
"""

function solar_array_type(PV_placement::String, land_available::Real)
    """
    In REopt the array types are as follows:
    0: Ground Mount Fixed (Open Rack); 
    1: Rooftop, Fixed; 
    2: Ground Mount 1-Axis Tracking; 
    3 : 1-Axis Backtracking; 
    4: Ground Mount, 2-Axis Tracking... #We will likely only use 0 and 1
    #maybe 2 once SAM runs are done
    """
    
    tilt_formula =

    # attain location for solar from [PV][location] output
    if PV_placement = "roof"
        array_type = 1
        tilt = #formula of tilt depending on latitude or keep as 20 or 15 degrees?
    elseif PV_placement = "ground"
        array_type = 0
        tilt = #formula of tilt depending on latitude or keep as 20 or 15 degrees?
    elseif PV_placement = "ground" && land_available > 3
        array_type = 2
        tilt = #formula of tilt depending on latitude or keep as 20 or 15 degrees?
    elseif PV_placement = "both"
        array_type = 1
        tilt = #formula of tilt depending on latitude or keep as 20 or 15 degrees?