import glob
import os
import pickle
from collections import Counter

import pandas as pd
from Bio import SeqIO

# %% Creation of the table with genomic islands per genome
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/IslandViewer4/GenomicIslands/')

dictionary_index = {}
for file in glob.glob('*gff'):
    try:
        table = pd.read_csv(file, sep='\t', skiprows=1, header=None)
        dictionary_index.update({file.replace('.gff', ''): [list(table[3]), list(table[4])]})
    finally:
        continue

# %% Import of dictionary with gene clusters ID per locustag
Dictionary1 = open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                   "Experiments/Experiments_Python/Dict_Cluster.pickle", "rb")
Clusters_GH = pickle.load(Dictionary1)

# %% Extraction of genes in the regions of putative genomic islands

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'Genome_Database/FinalSelectionGenomes/Prokka_selection/')
list_to_df = []
for file in glob.glob('*gbk'):
    infile = file.replace('.gbk', '')
    if infile in dictionary_index.keys():
        for i in range(len(dictionary_index[infile][0])):
            GI_beginning = dictionary_index[infile][0][i]
            GI_end = dictionary_index[infile][1][i]
            genbank = SeqIO.parse(file, 'genbank')
            for genome in genbank:
                for feature in genome.features:
                    if (feature.type == 'CDS' or feature.type == 'tRNA') and \
                            (
                                    feature.location._start.position > GI_beginning and feature.location._end.position < GI_end):
                        gen_feature = feature.type
                        ID = str(feature.qualifiers["locus_tag"][0]).replace(' ', '')
                        product = str(feature.qualifiers["product"][0])
                        start_pos = feature.location._start.position
                        end_pos = feature.location._end.position
                        if gen_feature == 'CDS':
                            list_to_df.append([ID, infile, gen_feature, Clusters_GH[ID], product, start_pos, end_pos])
                        else:
                            list_to_df.append([ID, infile, gen_feature, 'tRNA', product, start_pos, end_pos])

dataframe = pd.DataFrame(list_to_df, columns=['Locustag', 'Strain', 'Gene_feature',
                                              'ClusterID', 'Product', 'Start', 'End'])

# %% Analysis of frequency of genomic island per gene

clusters = dataframe['ClusterID'].value_counts()

clusters = clusters.drop(['tRNA'])

Genes = []
for i in clusters.index:
    Genes.append(str(i).split('_', 1)[1])

Genes = Counter(Genes)
Genes.pop('hypothetical_protein')
