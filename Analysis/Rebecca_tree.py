import os

import dendropy
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'MHC_Scripts/Phylogenetic_congruency/Rebecca/')

# %%
table_accessions = pd.read_excel('FinalPhageList_TaxonIDs.xlsx', header=None, sheet_name='Sheet1')
list_taxa = list(table_accessions[0])

fastas = SeqIO.parse('AlltRNAsCorrected20520.fasta', 'fasta')
with open('fasta_trnas_phage_reb.fasta', 'a') as file:
    for genome in fastas:
        if str(genome.id).split('_', 1)[0] in list_taxa:
            print(genome.id)
            aa = SeqRecord(genome.seq, id=genome.id, description='')
            SeqIO.write(aa, file, "fasta")

# %%
Tree = dendropy.Tree.get(path='updatedTreeAttempt.newick', schema='newick')

List_taxa = Tree.taxon_namespace.labels()

phage = pd.read_csv('AccessionNumbers.csv', header=None)

list_phage = list(phage[0])

phage_taxa = set(taxon for taxon in Tree.taxon_namespace
                 if taxon.label in list_phage)

Tree.retain_taxa_with_labels(list_phage)

with open('List_of_taxa_updatedTreeAttempt.txt', 'a') as file:
    for i in List_taxa:
        file.write(str(i) + '\n')
