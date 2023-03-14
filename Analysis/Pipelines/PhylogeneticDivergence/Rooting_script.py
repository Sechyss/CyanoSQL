import os

from Bio import Phylo

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/Phylogenetic_tree/CorePhylo/List_trees_edited_No_freshwater')
# Run through all files in the folder
for file in os.listdir('/Users/u1866168/Documents/OneDrive - University of '
                       'Warwick/Experiments/Phylogenetic_tree/CorePhylo/List_trees_edited_No_freshwater'):
    if file.endswith('.treefile'):
        tree = Phylo.read(file, 'newick')  # Open the different newick trees
        tree.root_with_outgroup('Synechococcus_sp._RCC307')  # Outgroup the different files
        Phylo.write(tree, file.rsplit(".", 1)[0] + ".phyloxml", 'phyloxml')  # Close theady rooted

for file in os.listdir('/Users/u1866168/Documents/OneDrive - University of '
                       'Warwick/Experiments/Phylogenetic_tree/CorePhylo/List_trees_edited_No_freshwater'):

    if file.endswith('.phyloxml'):
        Phylo.convert(file, 'phyloxml', file.rsplit(".", 1)[0] + ".nex", 'nexus')
