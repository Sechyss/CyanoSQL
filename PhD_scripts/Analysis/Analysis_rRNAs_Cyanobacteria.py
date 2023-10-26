import csv
import os

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# %% Extraction of genome sequence and creation of fasta file


os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Temporal')

for file in os.listdir():
    if file.endswith('.gb'):
        print(str(file.replace('.gb', '')))
        for genome in SeqIO.parse(file, 'genbank'):
            locus = str(genome.name)
            db = open(file.replace('.gb', '') + '_' + locus + '.fasta', 'w')
            cluster = str(file.replace('.gb', ''))
            id = str(cluster) + '_' + locus
            aaseq = genome.seq
            aa = SeqRecord(aaseq, id=id, description='')
            SeqIO.write(aa, db, "fasta")
            db.close()
# %%
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Temporal')

for file in os.listdir():
    if file.endswith('.gb'):
        print(str(file.replace('.gb', '')))
        db = open(file.replace('.gb', '') + '.fasta', 'a')
        for genome in SeqIO.parse(file, 'genbank'):
            id = genome.name
            aaseq = genome.seq
            aa = SeqRecord(aaseq, id=id, description='')
            SeqIO.write(aa, db, "fasta")

# %% Dataframe of the results from Barnapp

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/')

barnapp = pd.read_csv('Barnapp_cyano.txt', sep='\t', header=None)
b = []
df_dict = {}
for index, row in barnapp.iterrows():
    genome = str(row[0])
    print(genome)
    if '16S' in row[8]:
        if genome not in b:
            b.append(genome)
            df_dict.update({genome: [[row[3], row[4]]]})
        else:
            df_dict[genome].append([row[3], row[4]])

with open('Barnapp_simplify_16S.csv', 'w') as f:  # Just use 'w' mode in 3.x
    w = csv.DictWriter(f, df_dict.keys())
    w.writeheader()
    w.writerow(df_dict)
