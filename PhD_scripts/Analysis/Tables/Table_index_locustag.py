import os

import pandas as pd
from MyPackage.GenbankEditing import GenbankEditing

# %% Extraction of information
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cya_clusterPipeline')

for file in os.listdir():
    if file.endswith('.gb'):
        outfile = open('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                       'Genome_Database/Locustag_indexes/' + file.replace('.gb', '_reference.fasta'), 'a')
        gbk = GenbankEditing(file)
        gbk.protein_fasta_creation(outfile)
        print(str(file.replace('.gb', '_reference.fasta')))
        outfile.close()

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'Genome_Database/FinalSelectionGenomes/Prokka_selection')

for file in os.listdir():
    if file.endswith('.gbk'):
        outfile = open('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                       'Genome_Database/Locustag_indexes/' + file.replace('.gbk', '_prokka.fasta'), 'a')
        gbk = GenbankEditing(file)
        gbk.protein_fasta_creation(outfile)
        print(str(file.replace('.gb', '_prokka.fasta')))
        outfile.close()

# --- Using the files is necessary to create individual blastp databases using the references and then query the prokka

# %% Creation of the SQL table with index and original locustags

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Locustag_indexes/Blastp_output')

Table_index_dict = {}

for file in os.listdir():
    genome_blast = pd.read_csv(file, sep='\t', header=None, index_col=0)
    for index, row in genome_blast.iterrows():
        Table_index_dict.update({index: [row[1], row[2], row[10]]})

Table_index = pd.DataFrame.from_dict(Table_index_dict, orient='index')

Table_index.to_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                   'MHC_Scripts/SQL_Tables/Locus_index_prokka.csv')

# %% Additional annotation from original genbank files

Table_index = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                          'MHC_Scripts/SQL_Tables/SQL_Locus_index_prokka.csv', index_col='NCBI_LT')

GenkeyOriginal = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                             'MHC_Scripts/Python_results/GeneKeyTable.csv', header=None, index_col=0)

subdf = GenkeyOriginal.loc[:, 2]

FinalTable = pd.merge(Table_index, subdf, how='left', left_index=True, right_index=True)

Table_index.to_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                   'MHC_Scripts/SQL_Tables/SQL_Locus_index_prokka_annotated.csv.csv')
