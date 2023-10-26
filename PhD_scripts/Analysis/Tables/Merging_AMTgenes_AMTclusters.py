import pandas as pd

from SQL_analysis import amt_merging_genes_clusters

SQL_table_full = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                             'Warwick/MHC_Scripts/SQL_Tables/SQL_Prokka_clusters_AMT_TPM.csv', index_col=0,
                             low_memory=True)
# --- Creation of cluster table for SQL database

Table_cluster = amt_merging_genes_clusters(SQL_table_full, "/Users/u1866168/Documents/OneDrive - "
                                                           "University of "
                                                           "Warwick/Experiments/Experiments_Python"
                                                           "/Dict_Cluster.pickle")
