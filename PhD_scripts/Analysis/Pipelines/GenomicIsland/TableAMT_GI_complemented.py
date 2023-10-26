import pandas as pd
from MyPackage import translate_cluster_locustag_id
from tqdm import tqdm

SQL_table_GI = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                           'MHC_Scripts/SQL_Tables/SQL_queries/GenomicIslandGenes_AMT.csv', index_col=0, low_memory=True)
SQL_table_full = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                             'Warwick/MHC_Scripts/SQL_Tables/SQL_Prokka_clusters_AMT_TPM.csv', index_col=0,
                             low_memory=True)

index_GI = SQL_table_GI.index

New_table = SQL_table_full
New_table['Cluster'] = ''
New_table['GenomicIsland'] = ''
New_table['Total_Expression'] = New_table.sum(axis=1)

for index, row in tqdm(New_table.iterrows()):
    New_table.loc[index, 'Cluster'] = translate_cluster_locustag_id(index, "/Users/u1866168/Documents/OneDrive - "
                                                                           "University of "
                                                                           "Warwick/Experiments/Experiments_Python"
                                                                           "/Dict_Cluster.pickle")
    if index in index_GI:
        New_table.loc[index, 'GenomicIsland'] = 'True'
    else:
        New_table.loc[index, 'GenomicIsland'] = 'False'

New_table.to_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                 'MHC_Scripts/SQL_Tables/SQL_queries/GenomicIslandGenes_AMT_complemented.csv')
