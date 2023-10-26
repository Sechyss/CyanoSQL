#!/usr/bin/env python
import glob
import os
import pickle
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
import seaborn as sns
from MyPackage.PhyloDiversity import PhyloDiversity, mantel_pairwise
from SQL_analysis import tree_distancematrix_normalized
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from tqdm import tqdm

# %% Analysis Phylogeny of genes evolution. Creation of distances and collection
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/'
         'Phylogenetic_congruency/Prokka_Run/Divergence_evolution')
Dictionary_dataframes = {}
Dictionary_taxa = {}

for file1 in os.listdir():
    if file1.endswith('.contree'):
        # --- Creation of the first distance matrix for mantel
        pdm_gene, taxa_list = tree_distancematrix_normalized(file1)
        taxa_list.sort()
        print('Storaging the following gene ' + str(file1.rsplit(".", 4)[0]))
        Dictionary_dataframes.update({file1.rsplit(".", 4)[0]: pdm_gene})
        Dictionary_taxa.update({file1.rsplit(".", 4)[0]: taxa_list})

with open("../Dict_dataframe_dist.pickle", "wb") as wfp:
    pickle.dump(Dictionary_dataframes, wfp)

with open("../Dict_dataframe_taxa.pickle", "wb") as wfp:
    pickle.dump(Dictionary_taxa, wfp)

# %% --- Loading of the different databases from previous scripts

os.chdir('/Volumes/hcwebdav/M-Drive/HCSS6/Shared309/Alberto/MHC_Scripts/Phylogenetic_congruency/Prokka_phydiver_tmpfiles')
Dict1 = open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
             "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Dict_dataframe_dist.pickle", "rb")
df_database = pickle.load(Dict1)  # Input from the first part of the script
Dict2 = open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
             "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Dict_dataframe_taxa.pickle", "rb")
taxa_database = pickle.load(Dict2)  # Input from the second part of the script

# --- Creation of parameters for input and collection

keysdf = list(df_database.keys())
df_storage = np.empty([1, 3])
phylogeny = PhyloDiversity(keysdf, df_database, taxa_database, df_storage)

del Dict1, Dict2, df_database, taxa_database

# %% Creation of the multiplex array containing all matrix extracted

unique_list = []
print('Creating files with input')
for i in tqdm(keysdf):
    collector_lists = np.empty([1, 5], dtype=object)
    phylogeny.extraction_matrix(i, collector_lists, unique_list)
    collector_lists = np.delete(collector_lists, 0, 0)
    unique_list.append(i)

    with open('./tmpfiles/' + str(i) + '.data', 'ab') as collector_file:
        pickle.dump(collector_lists, collector_file)

del unique_list, collector_file, collector_lists, keysdf
# %% Mantel of the multiplex array containing all matrix extracted

# --- Mantel of the multiplex array containing all matrix extracted

if __name__ == '__main__':
    p = Pool(6)

    print('Running pairwise comparison of all the elements in our multiplex array')

    results = p.map(mantel_pairwise, tqdm(glob.glob('*data')))
    p.close()
    p.join()

    with open('../Recollection_object.pickle', 'wb') as outfile:
        pickle.dump(results, outfile)

# %% Analysis of the results for Mantel

Dict2 = open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
             "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Dict_dataframe_taxa.pickle", "rb")
df_database = pickle.load(Dict2)  # Input from the first part of the script
keysdf = list(df_database.keys())
os.chdir('/Volumes/HCSS6/Shared309/Alberto/MHC_Scripts/Phylogenetic_congruency/Test/')

raw_data = pickle.load(open('Recollection_object.pickle', 'rb'))

recollection_df = pd.DataFrame(columns=keysdf, index=keysdf)

for gene in tqdm(raw_data):
    for comparison in gene:
        recollection_df[comparison[0][0]][comparison[0][1]] = comparison[0][2]
        recollection_df[comparison[0][1]][comparison[0][0]] = comparison[0][2]

del raw_data, keysdf, Dict2, comparison, gene

for index, row in tqdm(recollection_df.iterrows()):
    for gene in row.index:
        if recollection_df.loc[index, gene] != float:
            taxa1 = df_database[index]
            len1 = len(taxa1)
            taxa2 = df_database[gene]
            len2 = len(taxa2)

            intersec = list(set(taxa1) & set(taxa2))

            recollection_df.loc[index, gene] = (1 - ((2 * len(intersec)) / (len1 + len2)))

# %% Upload of the pre-worked datasets/matrices and correction

mantel_matrix_raw = pickle.load(open('Mantel_matrix.pickle', 'rb'))
mantel_corrected = mantel_matrix_raw.fillna(1)

with open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
          "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Mantel_corrected.csv") as outfile:
    pickle.dump(mantel_corrected, outfile)

del mantel_matrix_raw, mantel_corrected
# %% Exploratory Tables and figures

mantel_corrected = pd.read_csv("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                               "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Mantel_corrected.csv", header=0, index_col=0)
mantel_corrected.index = mantel_corrected.index.map(lambda x: str(x).strip(" .\"\'"))
mantel_corrected.columns = mantel_corrected.columns.map(lambda x: str(x).strip(" .\"\'"))


# %% Creation of the model for HAC

def plot_dendrogram(linkage_model, **kwargs):  # from sklearn workbook
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(linkage_model.children_.shape[0])
    n_samples = len(linkage_model.labels_)
    for INDEX, merge in enumerate(linkage_model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[INDEX] = current_count

    linkage_matrix = np.column_stack([linkage_model.children_, linkage_model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


model = AgglomerativeClustering(distance_threshold=150, n_clusters=None, affinity='euclidean', linkage='ward')
model = model.fit(mantel_corrected)
# plot the top three levels of the dendrogram
plot_dendrogram(model, truncate_mode='level', p=4, count_sort=True)
# plt.axhline(y=200, color='r', linestyle='--')
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()

# %% Definition of the dictionaries to work
dictionary_index_gene = {}
for index in mantel_corrected.index:
    dictionary_index_gene.update({list(mantel_corrected).index(index): index})

table_fract = pd.read_csv("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                          "MHC_Scripts/SQL_Tables/SQL_ClusterGHTable.csv", header=None)
Dictionary_fraction = pd.Series(table_fract[3].values, index=table_fract[2]).to_dict()
Dictionary_fraction = {k.strip(" .\"\'"): v.strip(" .\"\'") for k, v in Dictionary_fraction.items()}

cluster_dict = {}
for term, clusterid in enumerate(model.labels_):
    if clusterid not in cluster_dict.keys():
        cluster_dict.update({clusterid: [dictionary_index_gene[term].strip(" .\"\'")]})
    else:
        cluster_dict[clusterid].append(dictionary_index_gene[term].strip(" .\"\'"))

# %% Heatmap and dendrogram

# Raw dendrogram to explore the distances
data = mantel_corrected.values
plt.figure(figsize=(10, 7))
dend = shc.dendrogram(shc.linkage(data, method='ward', optimal_ordering=True))
plt.xlabel('Genes')
plt.ylabel('Euclidean distances')
plt.axhline(y=150, color='r', linestyle='--')
plt.show()

# %% linkage method to use for calculating clusters: ward
lut = dict(zip(pd.Series(Dictionary_fraction.values()).unique(), ("#191970", "#C6E2FF", "#191970")))

genes = pd.Series(Dictionary_fraction).map(lut)

dendrogram_sns = sns.clustermap(mantel_corrected, metric="euclidean", method="ward", row_colors=genes, col_colors=genes)
plt.show()

# %% Rich ideas

list_reordered = list(
    x for x in Dictionary_fraction.keys() if 'Core' in Dictionary_fraction[x] and x in mantel_corrected.index)
list_reordered.extend(x for x in Dictionary_fraction.keys() if x not in list_reordered and x in mantel_corrected.index)

mantel_reordered = mantel_corrected.reindex(columns=list_reordered, index=list_reordered)

dend_reordered = shc.linkage(mantel_reordered, method='ward', optimal_ordering=True)
clustering_dend = sns.clustermap(mantel_reordered, metric="euclidean", method="ward", row_colors=genes,
                                 col_colors=genes,
                                 row_linkage=dend_reordered, col_linkage=dend_reordered,
                                 xticklabels=False, yticklabels=False,
                                 cmap="YlGnBu")

plt.show()

# %% Extraction of clusters

clusters = shc.fcluster(dend_reordered, 150, 'distance')

# clusters indices correspond to incites of original df
dict_clust_gene = {}
for i, cluster in enumerate(clusters):
    if cluster not in dict_clust_gene.keys():
        dict_clust_gene.update({cluster: [[mantel_reordered.index[i],
                                           Dictionary_fraction[mantel_reordered.index[i]]]]})
    else:
        dict_clust_gene[cluster].append([mantel_reordered.index[i],
                                         Dictionary_fraction[mantel_reordered.index[i]]])
    # print(mantel_reordered.index[i], cluster)
# %% Analysis of the presence of core and shell in the different clusters
for cluster in dict_clust_gene:
    core_counter = 0
    shell_counter = 0
    for element in dict_clust_gene[cluster]:
        if 'Core' in element[1]:
            core_counter += 1
        else:
            shell_counter += 1

    print('There is a ' + str((core_counter / 1154) * 100) + '% of core genes in ' + str(
        cluster))
    print('     The percentage core genes in the cluster is ' + str(
        (core_counter / (core_counter + shell_counter)) * 100))

for element in dict_clust_gene[1]:
    if 'Shell' in element[1]:
        print(element[0])

cluster1 = dict_clust_gene[1]
result = set(list([row[2], row[3], row[4]] for index, row in table_fract.iterrows() if row[2] in cluster1))


# %% Function for analysis of clusters

def taxa_genes(inputlist):
    df_key = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                         'MHC_Scripts/SQL_Tables/SQL_GeneKeyTable_Prokka.csv', header=None)
    df_key[2] = df_key[2].map(lambda x: x.strip(" .\"\'"))
    df_clade = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                           'MHC_Scripts/SQL_Tables/SQL_Clade_table_prokka.csv')
    emptydict = {}
    for index, row in df_key.iterrows():
        if row[2] in inputlist:
            if row[2] not in emptydict.keys():
                emptydict.update({row[2]: [row[6]]})
            else:
                emptydict[row[2]].append(row[6])

    taxag = list(emptydict.values())
    temporal = []
    for i in taxag:
        for y in i:
            temporal.append(y)

    heatmap = pd.DataFrame(index=(emptydict.keys()), columns=df_clade['Strain'])

    for key in emptydict:
        for element in emptydict[key]:
            heatmap.loc[key][element] = 1
    heatmap = heatmap.fillna(0)

    return heatmap


# %% Comparison G3 vs G4 relatively close congruency

# Remove any dot that may affect the index matching
df_groups = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                        'MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Diverge_clustegram_analysis/Nodes.csv')
df_clade = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                       'MHC_Scripts/SQL_Tables/SQL_Clade_table_prokka.csv')

for i in range(1, 17):
    df_groups[str(i)] = df_groups[str(i)].map(lambda x: str(x).strip(" .\"\'"))

    inputlist = list(df_groups[str(i)])

    heatmap = taxa_genes(inputlist)
    res = sns.heatmap(heatmap, yticklabels=False, xticklabels=False, cbar=False)
    plt.savefig('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                'MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Diverge_clustegram_analysis/NodeFigures/Node_'+str(i)+'.pdf'
                )

    plt.show()

list_reordered_nodes = [3, 2, 1, 14, 11, 13, 12, 8, 7, 9, 10, 4, 5, 6, 15, 16]
collection_df = pd.DataFrame(columns=df_clade['Strain'])

for i in list_reordered_nodes:
    df_groups[str(i)] = df_groups[str(i)].map(lambda x: str(x).strip(" .\"\'"))

    inputlist = list(df_groups[str(i)])

    heatmap = taxa_genes(inputlist)
    collection_df = collection_df.append(heatmap)

res= sns.heatmap(collection_df, yticklabels=False, xticklabels=False, cbar=False)