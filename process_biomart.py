import os
import subprocess

folder_path = '/Volumes/Ohnolog/Ohnologs_v3.1/1_All_Genes/1_BioMart_gene_attributes/'

output_file_path = 'all_unique_chromosomes.txt'


with open(output_file_path, 'w') as output_file:

    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):  # process only .txt files
            file_path = os.path.join(folder_path, filename)
            
            # run the `cut -f3` command to extract the 3rd column from the file
            result = subprocess.run(f"cut -f3 {file_path}", shell=True, capture_output=True, text=True)
            
            # Store unique values in a set
            unique_values = set()
            
            # Process each line in the output of the `cut -f3` command
            for line in result.stdout.splitlines():
                if line.strip() != "chromosome_name":  # Skip header if present
                    unique_values.add(line.strip())
            
            # Write the unique values for this species to the output file
            output_file.write(f"Unique values from {filename}:\n")
            output_file.write("\n".join(sorted(unique_values)))  # Sort for readability
            output_file.write("\n\n")  # Add a space between results from different files

print(f"All unique values have been written to {output_file_path}")
