import os
import pandas as pd

def combine_csv_files(directory, output_file):
    # Initialize a list to hold DataFrames
    dataframes = []
    
    # Loop through each file in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            # Construct the full file path
            file_path = os.path.join(directory, filename)
            
            try:
                # Read the CSV file in chunks to avoid memory issues for large files
                df = pd.read_csv(file_path)
                
                # Append the DataFrame to the list (avoid appending one by one to the final file)
                dataframes.append(df)
            
            except Exception as e:
                print(f"Error reading {filename}: {e}")
                continue
    
    if dataframes:
        # Concatenate all DataFrames in the list into one DataFrame
        combined_df = pd.concat(dataframes, ignore_index=True)

        # Save the combined DataFrame to a new CSV file
        combined_df.to_csv(output_file, index=False)
        print(f"Combined data saved to {output_file}")
    else:
        print("No valid CSV files found or there were errors reading the files.")

if __name__ == "__main__":
    # Example usage
    directory = 'C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/results_turbine_0/'  # Replace with your directory path
    output_file = 'C:/Users/dbernal/Documents/GitHub/Onsite_Analysis/results/Wind/test_newlogic_2.csv'  # Replace with your desired output file name
    combine_csv_files(directory, output_file)
