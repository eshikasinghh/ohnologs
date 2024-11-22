import pandas as pd
import os

def filter_and_merge_ohno_pairs(species_list):
    # Define column names based on the file format
    columns = ['Ohno1', 'Ohno2', 'Outgroup Support', 'Multiplication for P1(>=k)', 
               'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 
               'dmelanogaster', 'P-self(>=k)']
    
    # Define base directories
    base_dir1 = '/Volumes/Ohnolog/Ohnologs_v3.1/5_Run_synteny/5_CombineSelfOutgroup'
    base_dir2 = '/Volumes/Ohnolog/Ohnologs_v3.1/5_Run_synteny/5_CombineSelfOutgroup100-400'
    
    for species in species_list:
        # Create full paths for each species
        file1_path = os.path.join(base_dir1, species, f"{species}Ohno_Self+Outgp.txt")
        file2_path = os.path.join(base_dir2, species, f"{species}Ohno_Self+Outgp.txt")
        
        # Check if both files exist
        if not (os.path.exists(file1_path) and os.path.exists(file2_path)):
            print(f"Skipping {species}: One or both files not found")
            continue
        
        # Read files
        try:
            df1 = pd.read_csv(file1_path, delim_whitespace=True, names=columns, skiprows=1)
            df2 = pd.read_csv(file2_path, delim_whitespace=True, names=columns, skiprows=1)
        except Exception as e:
            print(f"Error reading files for {species}: {e}")
            continue
        
        # Normalize Ohno pairs to be unordered by sorting
        df1['OhnoPair'] = df1[['Ohno1', 'Ohno2']].apply(lambda x: tuple(sorted(x)), axis=1)
        df2['OhnoPair'] = df2[['Ohno1', 'Ohno2']].apply(lambda x: tuple(sorted(x)), axis=1)
        
        # Find matching Ohno pairs
        ohno_pairs1 = set(df1['OhnoPair'])
        ohno_pairs2 = set(df2['OhnoPair'])
        matching_pairs = ohno_pairs1.intersection(ohno_pairs2)
        
        # Filter rows with matching Ohno pairs in both files
        filtered_df1 = df1[df1['OhnoPair'].isin(matching_pairs)].copy()
        filtered_df2 = df2[df2['OhnoPair'].isin(matching_pairs)].copy()

        if filtered_df1.empty or filtered_df2.empty:
            print(f"No matching pairs found for {species}")
            continue

        # Merge the data by matching `OhnoPair`
        merged_df = pd.merge(
            filtered_df1[['Ohno1', 'Ohno2', 'OhnoPair', 'Multiplication for P1(>=k)']].rename(
                columns={'Multiplication for P1(>=k)': '100-500_Multiplication'}
            ),
            filtered_df2[['Ohno1', 'Ohno2', 'OhnoPair', 'Multiplication for P1(>=k)']].rename(
                columns={'Multiplication for P1(>=k)': '100-400_Multiplication'}
            ),
            on='OhnoPair',
            how='inner'
        )

        # Extract only the desired columns for output
        output_data = merged_df[['Ohno1_x', 'Ohno2_x', '100-500_Multiplication', '100-400_Multiplication']]
        output_data = output_data.rename(columns={'Ohno1_x': 'Ohno1', 'Ohno2_x': 'Ohno2'})

        # Output the filtered and merged data
        output_path = os.path.join(base_dir2, species, f"{species}.ohnopairs")
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        output_data.to_csv(output_path, sep='\t', index=False, header=True)

        print(f"Processed {species}: {len(output_data)} rows saved")

# List of species to process
species_list = ['cjacchus', 'cfamiliaris', 'hsapiens', 'rnorvegicus']

# Run the function
filter_and_merge_ohno_pairs(species_list)
