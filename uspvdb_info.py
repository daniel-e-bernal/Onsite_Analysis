"""
This file is only to read the uspvb file which is in the USGS folder.
"""

import pandas as pd
import glob
import matplotlib.pyplot as plt

file_uspvdb = "C:\\Users\\dbernal\\Documents\\Github (esque)\\DOE Reference Buildings Comparison\\USGS\\uspvdbCSV\\uspvdb_v2_0_20240801.csv"
#C:\Users\dbernal\Documents\Github (esque)\DOE Reference Buildings Comparison\USGS\uspvdbCSV
df = pd.read_csv(file_uspvdb)
print("Successfully read ", file_uspvdb)

fixed_tilt = []
single_axis = []

filtered_df = df[df['p_agrivolt'] == 'non-agrivoltaic']
filtered_df = filtered_df[filtered_df['p_year'] >= 2014]

for index, row in filtered_df.iterrows():
    if row['p_axis'] == 'single-axis':
        single_axis.append(row)
    elif row['p_axis'] == 'fixed-tilt':
        fixed_tilt.append(row)
        
fixed_tilt_df = pd.DataFrame(fixed_tilt)
single_axis_df = pd.DataFrame(single_axis)

print(fixed_tilt_df.columns)
print(single_axis_df.columns)

fixed_tilt_df['p_power_density_dc'] = fixed_tilt_df['p_cap_dc'] / fixed_tilt_df['p_area'] * 1000000
single_axis_df['p_power_density_dc'] = single_axis_df['p_cap_dc'] / single_axis_df['p_area'] * 1000000

print(fixed_tilt_df.columns)
print(single_axis_df.columns)

#Prepare dataframes for plotting
fixed_tilt_plot = fixed_tilt_df.groupby('p_year')['p_power_density_dc'].agg(['mean', 'median'])
single_axis_plot = single_axis_df.groupby('p_year')['p_power_density_dc'].agg(['mean', 'median'])

# Graph 1: fixed tilt
plt.figure(figsize=(10, 6))
plt.plot(fixed_tilt_plot.index, fixed_tilt_plot['mean'], marker='o', linestyle='-', color='blue', label='Mean')
plt.plot(fixed_tilt_plot.index, fixed_tilt_plot['median'], marker='o', linestyle='--', color='red', label='Median')
plt.xlabel('Year')
plt.ylabel('Power Density DC (MW/km^2)')
plt.title('Power Density DC Over the Years (Fixed Tilt)')
plt.grid(True)
plt.legend()

#show every year in the x axis
plt.xticks(fixed_tilt_plot.index)

# Save the figure
output_directory = "C:\\Users\\dbernal\\Documents\\Github (esque)\\DOE Reference Buildings Comparison\\pictures\\"
file_name = "fixed_tilt_power_density"
#plt.savefig(f"{output_directory}{file_name}.png", format='png', dpi=600)  # You can change the filename and format as needed
plt.close()  # Close the figure if you're done with it

plt.show()

# Graph 2: single axis
plt.figure(figsize=(10, 6))
plt.plot(fixed_tilt_plot.index, single_axis_plot['mean'], marker='o', linestyle='-', color='blue', label='Mean')
plt.plot(fixed_tilt_plot.index, single_axis_plot['median'], marker='o', linestyle='--', color='red', label='Median')
plt.xlabel('Year')
plt.ylabel('Power Density DC (MW/km^2)')
plt.title('Power Density DC Over the Years (Single-Axis)')
plt.grid(True)
plt.legend()

#show every year in the x axis
plt.xticks(fixed_tilt_plot.index)

#save new fig under new name to reflect type of rack system
#file_name = "single_axis_power_density"
#plt.savefig(f"{output_directory}{file_name}.png", format='png', dpi=600)  # You can change the filename and format as needed
plt.close()  # Close the figure if you're done with it

plt.show()