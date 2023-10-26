import os

import pandas as pd
from Bio import SeqIO

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/GET_HOMOLOGUES/GH_Prokka_run/OMCL_alg/TranslatorX/')

# %% Creation of a dictionary which will transform partition file into cluster
Counter = 1
DictionaryAlignCluster = {}
filelist = os.listdir()
for file in sorted(filelist):
    if file.endswith('.nt_ali.fasta'):
        f2 = SeqIO.parse(file, 'fasta')
        Gene = str(file.replace('.transx.nt_ali.fasta', ''))
        Length = 0
        for seq_record in f2:
            Length = len(seq_record.seq) - 1
            break
        DictionaryAlignCluster.update({Gene: [Counter, (Counter + Length)]})
        Counter = Counter + Length + 1

del Counter
# %%  Parsing the table obtained by RDP4

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/Phylogenetic_congruency/')
Table_clusters_Core = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                                  'MHC_Scripts/Phylogenetic_congruency/RDP_analysis.csv', index_col=None)

# Table needs to be treated to remove spaces in the columns and excess of rows at the beginning.
# Also is important to remove * : sed -i .backup 's/\*//g' RDP_analysis.csv

# --- Setting variables before parsing the table.
Index_storage = {}
Parental_List = []
Temporal_gene = []
Recombination_event = ['Start']
# --- Parsing table and extracting the important features
for index, row in Table_clusters_Core.iterrows():

    if str(row['Begin']) != 'nan':

        for i in DictionaryAlignCluster.keys():
            if row['Begin'] > DictionaryAlignCluster[i][0] and row['End'] < DictionaryAlignCluster[i][1] \
                    and (row['Begin'] < row['End']):
                Temporal_gene = i
                Recombination_event = row["RecombinationEventNumber"]
                Index_storage.update({row["RecombinationEventNumber"]: [row["RDP"], row["GENECONV"],
                                                                        row["Bootscan"], row["Maxchi"],
                                                                        row["Chimaera"], row["SiSscan"],
                                                                        row["3Seq"]]})
                Parental_List.append([int(Recombination_event),
                                      i,
                                      row["Recombinant_Sequence(s)"],
                                      row["Minor_Parental_Sequence(s)"],
                                      row["Major_Parental_Sequence(s)"]])

    elif Recombination_event == row["RecombinationEventNumber"] and str(row['Begin']) == 'nan':

        Parental_List.append([int(Recombination_event),
                              Temporal_gene,
                              row["Recombinant_Sequence(s)"],
                              row["Minor_Parental_Sequence(s)"],
                              row["Major_Parental_Sequence(s)"]])

# --- Creation of tables and csv files
Table_RDP = pd.DataFrame.from_dict(Index_storage, orient='index', columns=["RDP", "GENECONV",
                                                                           "Bootscan", "Maxchi", "Chimaera", "SiSscan",
                                                                           "3Seq"])
Table_Recom = pd.DataFrame(Parental_List, columns=['Recombination', 'Gene', 'Recombinantseq', 'Minor Parental seq',
                                                   'Major Parental Sequence'])

Table_RDP.to_csv('/Users/u1866168/Documents/OneDrive - University of '
                 'Warwick/MHC_Scripts/Phylogenetic_congruency/RDP4_TableSQL.csv')
Table_Recom.to_csv('/Users/u1866168/Documents/OneDrive - University of '
                   'Warwick/MHC_Scripts/Phylogenetic_congruency/RDP_ParentalSeqsSQL.csv')
