import pandas as pd
import glob
import os

def extract_sample_name(file_path):
    # Extract the filename without the path
    filename = os.path.basename(file_path)
    # Extract sample name part before the first underscore (e.g., SM1003_MT_Bulk)
    sample_name = filename.split('_gea.tab')[0]
    return sample_name

def transform_and_rename(file_path):
    # Extract sample name from the filename
    sample_name = extract_sample_name(file_path)
    
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep='\t')
    
    # Apply transformations
    df.loc[df['mean_plus'] < 1, 'reads_cnt_plus'] = 0
    df.loc[df['covered_fraction_plus'] < 0.75, 'reads_cnt_plus'] = 0
    df.loc[df['mean_minus'] < 1, 'reads_cnt_minus'] = 0
    df.loc[df['covered_fraction_minus'] < 0.75, 'reads_cnt_minus'] = 0
    
    # Select required columns
    df = df[['gene_id', 'scaffold_id', 'strand', 'length', 'reads_cnt_plus', 'reads_cnt_minus', 'product']]
    
    # Rename columns to include sample name
    df.rename(columns={
        'reads_cnt_plus': f'{sample_name}_reads_cnt_plus',
        'reads_cnt_minus': f'{sample_name}_reads_cnt_minus'
    }, inplace=True)
    
    return df

# Get all file paths
file_paths = glob.glob("/rcfs/projects/SoilSFA_r0/Data_by_thrust_area/3_Field_Metaphenomes/A3.4_Interkingdom_interactions/metaT_metaG_virome_analyses/95_85_clustered_RNA_viruses/gea_parser/gea_tab_out/*.tab")

# Transform and rename columns for each file
transformed_dfs = [transform_and_rename(file_path) for file_path in file_paths]

# Start merging DataFrames on common columns
for i, df in enumerate(transformed_dfs):
    if i == 0:
        final_df = df
    else:
        final_df = pd.merge(final_df, df, on=['gene_id', 'scaffold_id', 'strand', 'length', 'product'], how='inner')

# Drop the 'product' column if it is no longer needed
final_df.drop(columns=['product'], inplace=True)

# Save the final merged DataFrame
final_df.to_csv('vOTUs_MT_bulk_output.tab', index=False, sep='\t')
