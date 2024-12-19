import pandas as pd

# Read the CSV file into a DataFrame
df = pd.read_csv('/rcfs/projects/SoilSFA_r0/Data_by_thrust_area/3_Field_Metaphenomes/A3.4_Interkingdom_interactions/metaT_metaG_virome_analyses/iPHoP/out_iPHoP_3.4_interkingdom_DNA_RNA_iPHoP/parsed_Host_prediction_to_genome_m90.csv')

# Define a function to concatenate values with a comma
def concatenate_values(series):
    return ', '.join(series.dropna().astype(str))

# Group by the 'Virus' column
grouped = df.groupby('Virus').agg({
    'Host genome': concatenate_values,
    'Host taxonomy': concatenate_values,
    'Main method': concatenate_values,
    'Confidence score': concatenate_values,
    'Additional methods': concatenate_values
}).reset_index()

# Count occurrences of each 'Virus' ID
virus_counts = df['Virus'].value_counts()

# Create a new DataFrame where each 'Virus' ID appears only once
final_df = grouped.groupby('Virus').filter(lambda x: virus_counts[x.name] > 1)

# Add the unique 'Virus' entries without any modifications
unique_df = df.drop_duplicates('Virus', keep=False)

# Combine the processed duplicates and unique entries
combined_df = pd.concat([final_df, unique_df], ignore_index=True)

# Save the resulting DataFrame to a new CSV file
combined_df.to_csv('parsed_Host_prediction_to_genome_m90_duplicates_fixed.csv', index=False)
