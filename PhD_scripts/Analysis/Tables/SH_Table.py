import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from SQL_analysis import multi_robinsonfoulds_comparison, tree_distancematrix_normalized
from ete3 import Tree
from skbio.stats.distance import mantel

# %% Transformation of files into nexus

os.chdir('')
# Run through all files in the folder
for file in os.listdir():
    if file.endswith('.contree'):
        tree = Phylo.read(file, 'newick')  # Open the different newick trees
        tree.root_with_outgroup('Synechococcus_sp._RCC307')  # Outgroup the different files
        Phylo.write(tree, file.rsplit(".", 1)[0] + ".phyloxml", 'phyloxml')  # Close rooted tree file

for file in os.listdir():

    if file.endswith('.phyloxml'):
        Phylo.convert(file, 'phyloxml', file.rsplit(".", 1)[0] + ".nex", 'nexus')
# %% Loading of the core-genome tree parsididian distance matrices

core_pdm, ref_taxa = tree_distancematrix_normalized('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                                                    'MHC_Scripts/Phylogenetic_congruency/Concatenated_core_Prokka.newick')

# %% Provided a folder with gene trees to analyze, this script compare any gene tree with the core-genome tree

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Core/')
df_dict = {}
dict_mantel_accessory = {}
for i in os.listdir():
    if i.endswith('.contree'):
        pdm_gene, taxa_list = tree_distancematrix_normalized(i)
        taxa_list.sort()
        # extract distance between taxa from tree
        pdcsubtree_coregenome = core_pdm.loc[taxa_list, taxa_list]
        # Mantel Test
        mantel_test = mantel(pdcsubtree_coregenome, pdm_gene, method='pearson', permutations=999)
        dict_mantel_accessory.update({i.rsplit('.', 4)[0]: [mantel_test[0], mantel_test[1], mantel_test[2]]})
        df_dict.update({i.rsplit('.', 4)[0]: [mantel_test[0], mantel_test[1]]})
        # Update of the dictionary

mantel_dataframe = pd.DataFrame.from_dict(dict_mantel_accessory, orient='index', columns=['mantel_value', 'p_value'])
# dataframe = pd.DataFrame.from_dict(df_dict, orient='index_genome', columns=['Euclidean_Distance', 'WRF_Distance'])

# %% Robinson foulds test from accessory and core genome.

reference_tree = Tree('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                      'MHC_Scripts/Phylogenetic_congruency/Concatenated_core_Prokka.newick')
shell_softcore_list = os.chdir('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Shell-SoftCore')

multi_robinsonfoulds_comparison(shell_softcore_list, reference_tree, df_dict)

core_list = os.chdir('/Users/u1866168/Documents/OneDrive - University of '
                     'Warwick/MHC_Scripts/Phylogenetic_congruency/Prokka_Run/Core/')
multi_robinsonfoulds_comparison(core_list, reference_tree, df_dict)

dataframe = pd.DataFrame.from_dict(df_dict, orient='index', columns=['mantel_value', 'p_value', 'taxa',
                                                                     'Unweighted RF Distance',
                                                                     'Max Unweighted RF Distance', 'Percentage_dev',
                                                                     'Shared taxa'])

dataframe.drop(dataframe[dataframe.p_value > 0.05].index, inplace=True)

# %%
arrayX = np.array(dataframe['mantel_value'])
arrayY = np.array(dataframe['Percentage_dev'])
labels = dataframe.index

plt.scatter(arrayX, arrayY)

plt.annotate(labels, (arrayX, arrayY))

# %% Plot of the mantel coefficient and creation of the csv file

# Sort of the columns

mantel_dataframe_sorted = mantel_dataframe.sort_values(by=['mantel_value'])

x_label = mantel_dataframe_sorted.index_genome
y_label = mantel_dataframe_sorted['mantel_value']

plt.bar(x_label, y_label, align='center', alpha=0.5)
plt.xticks(x_label)
plt.ylabel('mantel_coeff')
plt.title('mantel_coeff_for_shell/core_genes')

plt.show()
mantel_dataframe_sorted.to_csv('Shell_core_mantel.csv')
