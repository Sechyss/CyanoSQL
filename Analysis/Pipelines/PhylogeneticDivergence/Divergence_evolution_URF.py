import os
import pickle

import pandas as pd
import scipy.cluster.hierarchy as shc
import seaborn as sns
from ete3 import Tree
from matplotlib import pyplot as plt
from tqdm import tqdm

# %% Multiprocessing the Robinson Foulds analysis (estimated time 7 days)


os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'MHC_Scripts/Phylogenetic_congruency/Prokka_Run/All_Trees')

done_tree = []
dict_storage = list([])
for file in tqdm(os.listdir()):
    if file.endswith('.contree'):
        for file2 in os.listdir():
            if file2.endswith('.contree') and (file2 not in done_tree):
                reference = Tree(file)
                accessory_tree = Tree(file2)
                rf = reference.robinson_foulds(accessory_tree, unrooted_trees=True)[0]
                rf_max = reference.robinson_foulds(accessory_tree, unrooted_trees=True)[1]
                shared_taxa = len(reference.robinson_foulds(accessory_tree, unrooted_trees=True)[2])
                if file == file2:
                    dict_storage.append([file.rsplit(".", 4)[0], file2.rsplit(".", 4)[0],
                                         0, shared_taxa])
                    done_tree.append(file)
                else:
                    if rf_max != 0:
                        percentage_dev = 1- (rf / rf_max)
                        dict_storage.append([file.rsplit(".", 4)[0], file2.rsplit(".", 4)[0],
                                             percentage_dev, shared_taxa])
                        done_tree.append(file)
                    else:
                        dict_storage.append([file.rsplit(".", 4)[0], file2.rsplit(".", 4)[0],
                                             1, shared_taxa])
                        done_tree.append(file)

with open('../Recollection_object_RobinsonFoulds.pickle', 'wb') as outfile:
    pickle.dump(dict_storage, outfile)

# %% Table of the raw results from previous script
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'MHC_Scripts/Phylogenetic_congruency/Prokka_Run/')
raw_data = pickle.load(open('Recollection_object_RobinsonFoulds.pickle', 'rb'))
keys_df = list(pickle.load(open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                                "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Dict_dataframe_taxa.pickle",
                                "rb")).keys())
recollection_df = pd.DataFrame(columns=keys_df, index=keys_df)

for gene in tqdm(raw_data):
    recollection_df[gene[0]][gene[1]] = gene[2]
    recollection_df[gene[1]][gene[0]] = gene[2]

#%% Upload the matrix with NaNs START FROM THIS POINT
recollection_df = pd.read_csv("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                              "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/RobinsonFoulds_corrected.csv", index_col=0)
recollection_df = recollection_df.fillna(0)
#%% Preparation of the data to match with previous dictionaries
table_fract = pd.read_csv("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                          "MHC_Scripts/SQL_Tables/SQL_ClusterGHTable.csv", header=None)
Dictionary_fraction = pd.Series(table_fract[3].values, index=table_fract[2]).to_dict()
Dictionary_fraction = {k.strip(" .\"\'"): v.strip(" .\"\'") for k, v in Dictionary_fraction.items()}

dictionary_index_gene = {}
for index in recollection_df.index:
    dictionary_index_gene.update({list(recollection_df).index(index): index})


# %% linkage method to use for calculating clusters: ward
lut = dict(zip(pd.Series(Dictionary_fraction.values()).unique(), ("#191970", "#C6E2FF", "#191970")))

genes = pd.Series(Dictionary_fraction).map(lut)

dendrogram_sns = sns.clustermap(recollection_df, metric="euclidean", method="ward", row_colors=genes, col_colors=genes)
plt.show()

list_reordered = list(
    x for x in Dictionary_fraction.keys() if 'Core' in Dictionary_fraction[x] and x in recollection_df.index)
list_reordered.extend(x for x in Dictionary_fraction.keys() if x not in list_reordered and x in recollection_df.index)

recollection_df = recollection_df.reindex(columns=list_reordered, index=list_reordered)

dend_reordered = shc.linkage(recollection_df, method='ward', optimal_ordering=True)
clustering_dend = sns.clustermap(recollection_df, metric="euclidean", method="ward", row_colors=genes,
                                 col_colors=genes,
                                 row_linkage=dend_reordered, col_linkage=dend_reordered,
                                 xticklabels=False, yticklabels=False,
                                 cmap="YlGnBu")
plt.savefig("/Users/u1866168/Documents/OneDrive - University of Warwick/"
            "MHC_Scripts/Phylogenetic_congruency/Prokka_Run/RobinsonFoulds_heatmap.pdf")
