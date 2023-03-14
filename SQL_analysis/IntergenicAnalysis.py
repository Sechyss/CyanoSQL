"""
The functions in this script are meant to be a helper when analysing the phylogenetic divergence of different fractions
of genes vs the core genome (also consider the consensus/average evolution).

For this to work it needs the resulting trees from IQtree with the sufix of '.transx.nt_ali.fasta.contree'.

"""
import dendropy
import pandas as pd
from ete3 import Tree


def tree_distancematrix_normalized(treefile):
    tree = dendropy.Tree.get(path=treefile, schema='newick')
    tree.is_rooted = False

    taxon_list = tree.taxon_namespace.labels()  # Extract the list of taxa in the comparing tree

    # Create distance matrix from any tree
    pdc = tree.phylogenetic_distance_matrix()
    tree_pdm_normalized = pd.DataFrame.from_dict(pdc.as_data_table()._data).div(pdc.sum_of_distances())  # Normalization
    tree_pdm_normalized_sorted = tree_pdm_normalized.reindex(sorted(tree_pdm_normalized.columns),
                                                             axis=1)  # Sort columns
    tree_pdm_normalized_sorted = tree_pdm_normalized_sorted.sort_index()  # Sort rows alphabetically
    # This create a parsidian distance matrix normalized by the total sum of distances.

    return tree_pdm_normalized_sorted, taxon_list


def multi_robinsonfoulds_comparison(directory, reference, dict_storage):
    for file in directory:
        if file.endswith('.contree'):
            accessory_tree = Tree(file)
            rf = reference.robinson_foulds(accessory_tree, unrooted_trees=True)[0]
            rf_max = reference.robinson_foulds(accessory_tree, unrooted_trees=True)[1]
            shared_taxa = len(reference.robinson_foulds(accessory_tree, unrooted_trees=True)[2])
            percentage_dev = 1 - (rf / rf_max)
            dict_storage[file.rsplit(".", 4)[0]].append(rf)
            dict_storage[file.rsplit(".", 4)[0]].append(rf_max)
            dict_storage[file.rsplit(".", 4)[0]].append(percentage_dev)
            dict_storage[file.rsplit(".", 4)[0]].append(shared_taxa)
