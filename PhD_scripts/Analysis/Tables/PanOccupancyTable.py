import os

import pandas as pd
from Bio import SeqIO

# Working directory
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cya_clusterPipeline')
# Opening of the pangenome matrix created from GET_HOMOLOGUES
Occupancy = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                        'Warwick/Experiments/GET_HOMOLOGUES/FinalDraftGetHomologues/Matrices/pangenome_matrix_t0.tr'
                        '.csv', index_col='Gene')

ListGenomes = list(Occupancy)
for i in range(len(ListGenomes)):  # Translation of file name from genbank into their organisms
    file = ListGenomes[i]
    if file in os.listdir('/Users/u1866168/Documents/OneDrive - University of '
                          'Warwick/Genome_Database/Cya_clusterPipeline'):
        if not file.startswith('.'):
            gbank = SeqIO.parse(file, "genbank")
            if file.endswith('.gb'):
                for genome in gbank:
                    strain = genome.annotations['source']
                    ListGenomes[i] = strain
                    continue
            else:
                for genome in gbank:
                    for feature in genome.features:
                        if feature.type == "source":
                            strain = str(feature.qualifiers["organism"][0])  # Storage of the variable strain
                            ListGenomes[i] = strain
                            continue
Occupancy.columns = ListGenomes  # Creation of the new csv with organisms instead of file names
Occupancy.to_csv('/Users/u1866168/Documents/OneDrive - University of '
                 'Warwick/Experiments/Experiments_Python/pangenome_matrix_occupancy_corrected.csv')
