from rex import NSRDBX
from rex.sam_resource import SAMResource
import numpy as np 
import pandas as pd

# NOTE: creating these variables as placeholders, these will need to be updated to self.X attributes of the SolarResource class object
year = 2012
# flatirons_site.yaml
lat = 35.2018863
lon = -101.945027
#other site
# lat = 27.18624
# lon = -96.9516

# NOTE: Old deprecated version of PSM v3.1 which corresponds to /api/nsrdb/v2/solar/psm3-download (what HOPP is currently calling)
nsrdb_file = '/datasets/NSRDB/deprecated_v3/nsrdb_{year}.h5'.format(year=year)

# NOTE: Current version of PSM v3.2.2 which corresponds to /api/nsrdb/v2/solar/psm3-2-2-download 
# nsrdb_file = '/datasets/NSRDB/current/nsrdb_{year}.h5'.format(year=year)

with NSRDBX(nsrdb_file, hsds=False) as f:
    # get gid of location closest to given lat/lon coordinates
    site_gid = f.lat_lon_gid((lat,lon))

    # extract timezone, elevation, latitude and longitude from meta dataset with gid
    time_zone = f.meta['timezone'].iloc[site_gid]
    elevation = f.meta['elevation'].iloc[site_gid]
    latitude = f.meta['latitude'].iloc[site_gid]
    longitude = f.meta['longitude'].iloc[site_gid]
    
    # extract remaining datapoints: year, month, day, hour, minute, dn, df, gh, wspd,tdry, pres, tdew
    # NOTE: datasets have readings at 0 and 30 minutes each hour, HOPP/SAM workflow requires only 30 minute reading values -> therefore filter 0 minute readings with [1::2]
    # NOTE: datasets are not auto shifted by timezone offset -> therefore wrap extraction in SAMResource.roll_timeseries(input_array, timezone, #steps in an hour)
    # NOTE: solar_resource.py code references solar_zenith_angle and RH = relative_humidity but I couldn't find them being used in the code. Captured them below just in case.
    year_arr = f.time_index.year.values[1::2]
    month_arr = f.time_index.month.values[1::2]
    day_arr = f.time_index.day.values[1::2]
    hour_arr = f.time_index.hour.values[1::2]
    minute_arr = f.time_index.minute.values[1::2]
    dni_arr = SAMResource.roll_timeseries((f['dni', :, site_gid][1::2]), time_zone, 1)
    dhi_arr = SAMResource.roll_timeseries((f['dhi', :, site_gid][1::2]), time_zone, 1)
    ghi_arr = SAMResource.roll_timeseries((f['ghi', :, site_gid][1::2]), time_zone, 1)
    wspd_arr = SAMResource.roll_timeseries((f['wind_speed', :, site_gid][1::2]), time_zone, 1)
    tdry_arr = SAMResource.roll_timeseries((f['air_temperature', :, site_gid][1::2]), time_zone, 1)
    # relative_humidity_arr = SAMResource.roll_timeseries((f['relative_humidity', :, site_gid][1::2]), time_zone, 1)
    # solar_zenith_arr = SAMResource.roll_timeseries((f['solar_zenith_angle', :, site_gid][1::2]), time_zone, 1)
    pres_arr = SAMResource.roll_timeseries((f['surface_pressure', :, site_gid][1::2]), time_zone, 1)
    tdew_arr = SAMResource.roll_timeseries((f['dew_point', :, site_gid][1::2]), time_zone, 1)
    
# Remove data from feb29 on leap years
if (year % 4) == 0:
    feb29 = np.arange(1416,1440)
    year_arr = np.delete(year_arr, feb29)
    month_arr = np.delete(month_arr, feb29)
    day_arr = np.delete(day_arr, feb29)
    hour_arr = np.delete(hour_arr, feb29)
    minute_arr = np.delete(minute_arr, feb29)
    dni_arr = np.delete(dni_arr, feb29)
    dhi_arr = np.delete(dhi_arr, feb29)
    ghi_arr = np.delete(ghi_arr, feb29)
    wspd_arr = np.delete(wspd_arr, feb29)
    tdry_arr = np.delete(tdry_arr, feb29)
    # relative_humidity_arr = np.delete(relative_humidity_arr, feb29)
    # solar_zenith_arr = np.delete(solar_zenith_arr, feb29)
    pres_arr = np.delete(pres_arr, feb29)
    tdew_arr = np.delete(tdew_arr, feb29)

# round to desired precision and convert to desired data type
time_zone = float(time_zone)
elevation = round(float(elevation), 0)
latitude = round(float(latitude), 2)
longitude = round(float(longitude),2)
year_arr = list(year_arr.astype(float, copy=False))
month_arr = list(month_arr.astype(float, copy=False))
day_arr = list(day_arr.astype(float, copy=False))
hour_arr = list(hour_arr.astype(float, copy=False))
minute_arr = list(minute_arr.astype(float, copy=False))
dni_arr = list(dni_arr.astype(float, copy=False))
dhi_arr = list(dhi_arr.astype(float, copy=False))
ghi_arr = list(ghi_arr.astype(float, copy=False))
wspd_arr = list(wspd_arr.astype(float, copy=False))
tdry_arr = list(tdry_arr.astype(float, copy=False))
# relative_humidity_arr = list(np.round(relative_humidity_arr, decimals=1))
# solar_zenith_angle_arr = list(np.round(solar_zenith_angle_arr, decimals=1))
pres_arr = list(pres_arr.astype(float, copy=False))
tdew_arr = list(tdew_arr.astype(float, copy=False))




# Final dictionary for use in solar_resource.py
_data = {'tz' :     time_zone,
         'elev' :   elevation,
         'lat' :    latitude,
         'lon' :    longitude,
         'year' :   year_arr,
         'month' :  month_arr,
         'day' :    day_arr,
         'hour' :   hour_arr,
         'minute' : minute_arr,
         'dn' :     dni_arr,
         'df' :     dhi_arr,
         'gh' :     ghi_arr,
         'wspd' :   wspd_arr,
         'tdry' :   tdry_arr,
         'pres' :   pres_arr,
         'tdew' :   tdew_arr
        }

# Validation against HOPP/examples/01-wind-solar.ipynb > hi.system.site.solar_resource._data 
# test_df = pd.DataFrame.from_dict(_data)
# print(test_df.tail(15))
# print(type(_data['tdew']))
# print(type(_data['tdew'][10]))
# print(_data['tdew'][10])
# print(_data['lon'])

