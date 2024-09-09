from rex import WindX
from rex.sam_resource import SAMResource
import numpy as np 
# import pandas as pd

# NOTE: creating these variables as placeholders, these will need to be updated to self.X attributes of the WindResource class object
year = 2012
# flatirons_site.yaml
lat = 35.2018863
lon = -101.945027
hub_height_meters = 97.0
allowed_hub_height_meters = [10, 40, 60, 80, 100, 120, 140, 160, 200]

# NOTE: Current setup of files on HPC v1.0.0 = 2007-2013, v1.1.0 = 2014
if year < 2014:
    wtk_file = '/datasets/WIND/conus/v1.0.0/wtk_conus_{year}.h5'.format(year=year)
elif year == 2014:
    wtk_file = '/datasets/WIND/conus/v1.1.0/wtk_conus_{year}.h5'.format(year=year)

def calculate_heights_to_download():
    """
    Given the system hub height, and the available hubheights from WindToolkit,
    determine which heights to download to bracket the hub height
    """

    # evaluate hub height, determine what heights to download

    heights = [hub_height_meters]
    if hub_height_meters not in allowed_hub_height_meters:
        height_low = allowed_hub_height_meters[0]
        height_high = allowed_hub_height_meters[-1]
        for h in allowed_hub_height_meters:
            if h < hub_height_meters:
                height_low = h
            elif h > hub_height_meters:
                height_high = h
                break
        heights[0] = height_low
        heights.append(height_high)

    return heights

heights = calculate_heights_to_download()

with WindX(wtk_file, hsds=False) as f:
    # get gid of location closest to given lat/lon coordinates and timezone offset
    site_gid = f.lat_lon_gid((lat,lon))
    time_zone = f.meta['timezone'].iloc[site_gid]

    # instantiate dictionary to hold datasets
    wind_dict = {}
    # loop through hub heights to download, capture dataset and simultaneously shift the data by the timezone offset
    # NOTE: pressure datasets unit = Pa, convert to atm via division by 101325
    for h in heights:
        wind_dict['temperature_{height}m_arr'.format(height=h)] = SAMResource.roll_timeseries((f['temperature_{height}m'.format(height=h), :, site_gid]), time_zone, 1)
        wind_dict['pressure_{height}m_arr'.format(height=h)] = SAMResource.roll_timeseries((f['pressure_{height}m'.format(height=h), :, site_gid]/101325), time_zone, 1)
        wind_dict['windspeed_{height}m_arr'.format(height=h)] = SAMResource.roll_timeseries((f['windspeed_{height}m'.format(height=h), :, site_gid]), time_zone, 1)
        wind_dict['winddirection_{height}m_arr'.format(height=h)] = SAMResource.roll_timeseries((f['winddirection_{height}m'.format(height=h), :, site_gid]), time_zone, 1)    
    
# Remove data from feb29 on leap years
if (year % 4) == 0:
    feb29 = np.arange(1416,1440)
    for key, value in wind_dict.items():
        wind_dict[key] = np.delete(value, feb29)
    
# round to desired precision and concatenate data into format needed for data dictionary
if len(heights) == 2:
    wind_dict['temperature_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['temperature_{h}m_arr'.format(h=heights[0])]), decimals=1)
    wind_dict['pressure_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['pressure_{h}m_arr'.format(h=heights[0])]), decimals=2)
    wind_dict['windspeed_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['windspeed_{h}m_arr'.format(h=heights[0])]), decimals=3)
    wind_dict['winddirection_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['winddirection_{h}m_arr'.format(h=heights[0])]), decimals=1)
    wind_dict['temperature_{h}m_arr'.format(h=heights[1])] = np.round((wind_dict['temperature_{h}m_arr'.format(h=heights[1])]), decimals=1)
    wind_dict['pressure_{h}m_arr'.format(h=heights[1])] = np.round((wind_dict['pressure_{h}m_arr'.format(h=heights[1])]), decimals=2)
    wind_dict['windspeed_{h}m_arr'.format(h=heights[1])] = np.round((wind_dict['windspeed_{h}m_arr'.format(h=heights[1])]), decimals=3)
    wind_dict['winddirection_{h}m_arr'.format(h=heights[1])] = np.round((wind_dict['winddirection_{h}m_arr'.format(h=heights[1])]), decimals=1)
    combined_data = [list(a) for a in zip(wind_dict['temperature_{h}m_arr'.format(h=heights[0])],
                                            wind_dict['pressure_{h}m_arr'.format(h=heights[0])],
                                            wind_dict['windspeed_{h}m_arr'.format(h=heights[0])],
                                            wind_dict['winddirection_{h}m_arr'.format(h=heights[0])],
                                            wind_dict['temperature_{h}m_arr'.format(h=heights[1])],
                                            wind_dict['pressure_{h}m_arr'.format(h=heights[1])],
                                            wind_dict['windspeed_{h}m_arr'.format(h=heights[1])],
                                            wind_dict['winddirection_{h}m_arr'.format(h=heights[1])])]

elif len(heights) == 1:
    wind_dict['temperature_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['temperature_{h}m_arr'.format(h=heights[0])]), decimals=1)
    wind_dict['pressure_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['pressure_{h}m_arr'.format(h=heights[0])]), decimals=2)
    wind_dict['windspeed_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['windspeed_{h}m_arr'.format(h=heights[0])]), decimals=3)
    wind_dict['winddirection_{h}m_arr'.format(h=heights[0])] = np.round((wind_dict['winddirection_{h}m_arr'.format(h=heights[0])]), decimals=1)
    combined_data = [list(a) for a in zip(wind_dict['temperature_{h}m_arr'.format(h=heights[0])],
                                            wind_dict['pressure_{h}m_arr'.format(h=heights[0])],
                                            wind_dict['windspeed_{h}m_arr'.format(h=heights[0])],
                                            wind_dict['winddirection_{h}m_arr'.format(h=heights[0])])]

# NOTE: possible faster concatenation with np.concat however precision / dtype of all elements changes (use if concerned about speed and HOPP/NED can handle additional float precision)
# if len(heights) == 2:
#     combined_data = np.concatenate((wind_dict['temperature_{h}m_arr'.format(h=heights[0])].reshape(8760,1),
#                                     wind_dict['pressure_{h}m_arr'.format(h=heights[0])].reshape(8760,1),
#                                     wind_dict['windspeed_{h}m_arr'.format(h=heights[0])].reshape(8760,1),
#                                     wind_dict['winddirection_{h}m_arr'.format(h=heights[0])].reshape(8760,1),
#                                     wind_dict['temperature_{h}m_arr'.format(h=heights[1])].reshape(8760,1),
#                                     wind_dict['pressure_{h}m_arr'.format(h=heights[1])].reshape(8760,1),
#                                     wind_dict['windspeed_{h}m_arr'.format(h=heights[1])].reshape(8760,1),
#                                     wind_dict['winddirection_{h}m_arr'.format(h=heights[1])].reshape(8760,1)
#                                     ), axis=1).tolist()
# elif len(heights) == 1:
#     combined_data = np.concatenate((wind_dict['temperature_{h}m_arr'.format(h=heights[0])].reshape(8760,1),
#                                     wind_dict['pressure_{h}m_arr'.format(h=heights[0])].reshape(8760,1),
#                                     wind_dict['windspeed_{h}m_arr'.format(h=heights[0])].reshape(8760,1),
#                                     wind_dict['winddirection_{h}m_arr'.format(h=heights[0])].reshape(8760,1)
#                                     ), axis=1).tolist()


# Final dictionary for use in wind_resource.py
_data = {'heights': [float(h) for h in heights for i in range(4)],
         'fields':  [1, 2, 3, 4] * len(heights),
         'data':    combined_data
        }

# validation
# print(type(_data['data']))
# print(type(_data['data'][0]))
# print(_data['data'][0])
