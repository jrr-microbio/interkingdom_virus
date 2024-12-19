import pandas as pd

# Read the TSV file
df = pd.read_csv('/rcfs/projects/SoilSFA_r0/Data_by_thrust_area/3_Field_Metaphenomes/A3.4_Interkingdom_interactions/metaT_metaG_virome_analyses/95_85_clustered_RNA_viruses/gea_parser/vOTUs_MT_bulk_output.tab', sep='\t')  # Replace with your actual TSV file path

# Initialize a dictionary to hold data for the new table
combined_data = {}

# Iterate through the rows of the dataframe
for index, row in df.iterrows():
    gene_id = row['gene_id']
    strand = row['strand']
    
    # Initialize dictionary for gene_id if not already present
    if gene_id not in combined_data:
        combined_data[gene_id] = {'scaffold_id': row['scaffold_id'], 'length': row['length'], 'flag': ''}
    
    # Iterate through all columns and assign values to their respective abundance and expression columns
    for col in df.columns:
        if 'reads_cnt' in col:
            sample = col.split('_reads_cnt')[0]
            if f'{sample}_abundance' not in combined_data[gene_id]:
                combined_data[gene_id][f'{sample}_abundance'] = 0
                combined_data[gene_id][f'{sample}_expression'] = 0
                
            if strand == '+':
                if 'plus' in col:
                    combined_data[gene_id][f'{sample}_abundance'] = row[col]
                elif 'minus' in col:
                    combined_data[gene_id][f'{sample}_expression'] = row[col]
            elif strand == '-':
                if 'minus' in col:
                    combined_data[gene_id][f'{sample}_abundance'] = row[col]
                elif 'plus' in col:
                    combined_data[gene_id][f'{sample}_expression'] = row[col]
    
    # Check for 'careful' flag condition
    abundant = any(combined_data[gene_id][f'{sample}_abundance'] > 0 for sample in set(col.split('_reads_cnt')[0] for col in df.columns if 'reads_cnt' in col))
    expressed = any(combined_data[gene_id][f'{sample}_expression'] > 0 for sample in set(col.split('_reads_cnt')[0] for col in df.columns if 'reads_cnt' in col))
    
    if expressed and not abundant:
        combined_data[gene_id]['flag'] = 'careful - expressed but not abundant'

# Convert the dictionary to a DataFrame
combined_df = pd.DataFrame.from_dict(combined_data, orient='index').reset_index()
combined_df.rename(columns={'index': 'gene_id'}, inplace=True)

# Write the dataframe to a tsv file
combined_df.to_csv('vOTU_output_abund_exp_parsed.tsv', sep='\t', index=False)

print("File created successfully: vOTUs_output_abund_exp_parsed.tsv")
