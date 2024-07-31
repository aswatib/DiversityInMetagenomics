import os
import pandas as pd

# Directory containing Bracken output files
directory = '/home/people/s220868/project/sewage_data/kraken_run/bracken_temp_out'

# Initialize an empty DataFrame to store combined results
combined_df = pd.DataFrame()

# Iterate through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('_bracken.report'):

        # Get the full file path
        filepath = os.path.join(directory, filename)

        # Read the Bracken report file
        df = pd.read_csv(filepath, sep='\t')

        # Extract the file identifier from the filename
        file_id = filename.replace('_bracken.report', '')

        # Select the species name and new_est_reads columns, add the file identifier
#       df = df[['name', 'new_est_reads']]
        df = df[['name', 'kraken_assigned_reads']]
        df.columns = ['species', file_id]
        df.set_index('species', inplace=True)

        # Merge with the combined DataFrame
        if combined_df.empty:
            combined_df = df
        else:
            combined_df = combined_df.join(df, how='outer')

# Replace NA values with 0
combined_df.fillna(0, inplace=True)

# Reset the index to make species names as columns
combined_df.reset_index(inplace=True)

# Rename the index column to 'species'
combined_df.rename(columns={'index': 'species'}, inplace=True)

# Save the combined DataFrame to a CSV file
#combined_df.to_csv('/home/people/s220868/project/sewage_data/kraken_run/bracken_temp_out/combined_bracken_output.csv', index=True)
combined_df.to_csv('/home/people/s220868/project/sewage_data/kraken_run/bracken_temp_out/combined_kraken2_output.csv', index=True)
