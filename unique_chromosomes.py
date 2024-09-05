import pandas as pd

species_list = [
"acarolinensis_gene_ensembl",
"fcatus_gene_ensembl",
"ggallus_gene_ensembl",
"ptroglodytes_gene_ensembl",
"btaurus_gene_ensembl",
"clfamiliaris_gene_ensembl",
"ggorilla_gene_ensembl",
"ecaballus_gene_ensembl",
"hsapiens_gene_ensembl",
"mmulatta_gene_ensembl",
"cjacchus_gene_ensembl",
"mmusculus_gene_ensembl",
"panubis_gene_ensembl",
"mdomestica_gene_ensembl",
"pabelii_gene_ensembl",
"sscrofa_gene_ensembl",
"oanatinus_gene_ensembl",
"ocuniculus_gene_ensembl",
"rnorvegicus_gene_ensembl",
"oaries_gene_ensembl",
"mgallopavo_gene_ensembl",
"csabaeus_gene_ensembl",
"tguttata_gene_ensembl",
"loculatus_gene_ensembl",
"lchalumnae_gene_ensembl",
"gaculeatus_gene_ensembl",
"olatipes_gene_ensembl",
"tnigroviridis_gene_ensembl",
"drerio_gene_ensembl",
"trubripes_gene_ensembl",
"dmelanogaster_gene_ensembl",
"cintestinalis_gene_ensembl",
"csavignyi_gene_ensembl",
"celegans_gene_ensembl"
]

chromosome_dict = {}

for specie in species_list:
    file_path = f"/Volumes/Ohnolog/Ohnologs_v2/1_All_Genes_copy/new_1_BioMart_gene_attributes/{specie}_biomaRt_v112_20160703.txt"


    df = pd.read_csv(file_path, delimiter='\t')  # read the file using pandas, file is tab-delimited
    unique_chromosomes = df['chromosome_name'].unique()
    chromosome_dict[specie] = unique_chromosomes
#   print(unique_chromosomes)

df_chromosomes = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in chromosome_dict.items()]))
    # k = specie name
    # v = list of unique chromosomes
    #converts dictionary into datagrame

col_width = max(len(col) for col in df_chromosomes.columns) + 2
    #calculated length of the longest species name and adds 2 for extra space/readability

with open('unique_chromosomes.txt', 'w') as f: #opens file in write mode
    f.write(' '.join(f'{col:<{col_width}}' for col in df_chromosomes.columns) + '\n')
            # f-string formats each column name
            #< left-aligns the column name
            #.join joins the formatted columns into a single string separated by a space
            # \n is moves to the next line for the data rows to begin
    
    for row in df_chromosomes.itertuples(index=False, name=None):
            #itertuples with name = None mean we can access elements in each row by their position
        formatted_row = ' '.join(f'{str(cell):<{col_width}}' if pd.notna(cell) else ' ' * col_width for cell in row)
        f.write(formatted_row.rstrip() + '\n')