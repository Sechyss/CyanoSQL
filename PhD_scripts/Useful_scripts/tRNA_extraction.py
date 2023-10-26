import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

os.chdir('/Users/u1866168/Downloads')

table = pd.read_csv('All_features_Raw_reads.csv')


list = np.array([])

for i in table['Locustag']:
    if '-' in i and 'Syn' in i:
        list = np.append(list, i)


genome = SeqIO.parse('/Users/u1866168/Documents/OneDrive - University of '
                     'Warwick/Genome_Database/Synechococcus_Cyanorak/Syn_WH7803.gbk', 'genbank')
file = open('tRNA_sequences.fasta', 'w')
for i in genome:
    for feature in i.features:
        if feature.type == 'tRNA':
            locustag = str(feature.qualifiers["locus_tag"][0]).replace(' ', '_')
            if locustag in list:
                ntseq = feature.extract(i.seq)
                print(ntseq)
                fasta_aa = SeqRecord(ntseq, locustag, description='')
                SeqIO.write(fasta_aa, file, 'fasta')