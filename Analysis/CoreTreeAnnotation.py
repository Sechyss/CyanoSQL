import os

import pandas as pd
from ete3 import Tree, TreeStyle

from CyanoScripts.MyPackage import color_node, color_leaves

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/Phylogenetic_congruency')
tree = Tree('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/'
            'Phylogenetic_congruency/Concatenated_core_Prokka.newick')
clade_table = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                          'MHC_Scripts/SQL_Tables/SQL_Clade_table_prokka.csv')
dict_clade = {}

for index, row in clade_table.iterrows():
    if row['Clade'] in dict_clade.keys():
        dict_clade[row['Clade']].append(row['Strain'])
    else:
        dict_clade.update({row['Clade']: [row['Strain']]})

# %%  Setting colours per clade

dict_colors = {'SeaGreen': dict_clade['Pro_LLII-III'],
               'LimeGreen': dict_clade['Pro_LLI'],
               'ForestGreen': dict_clade['Pro_LLIV'],
               'GreenYellow': dict_clade['Pro_HLII'],
               'Lime': dict_clade['Pro_HLI'],
               'Cyan': dict_clade['Syn_5.1_I'],
               'Orange': dict_clade['Syn_5.1_II'],
               'LightSalmon': dict_clade['Syn_5.1_III'],
               'MediumTurquoise': dict_clade['Syn_5.1_IV'],
               'Violet': dict_clade['Syn_5.1_V-VI'],
               'DarkViolet': dict_clade['Syn_5.1_VII'],
               'Olive': dict_clade['Syn_5.1_VIII'],
               'RosyBrown': dict_clade['Syn_5.1_CRD']}

dict_colors_leaves = {'DarkKhaki': dict_clade['Syn_5.1_IX'],
                      'Wheat': dict_clade['Syn_5.1_X'],
                      'Peru': dict_clade['Syn_5.1_UCA'],
                      'SaddleBrown': dict_clade['Syn_5.1_WPC-1'],
                      'DarkGoldenrod': dict_clade['Syn_5.1_WPC-2'],
                      'Maroon': dict_clade['Brackish'],
                      'Black': dict_clade['Undefined']}

# --- Color the nodes by clade
for key in dict_colors.keys():
    node = tree.get_common_ancestor(dict_colors[key])
    color_node(node, key)

for key in dict_colors_leaves.keys():
    for value in dict_colors_leaves[key]:
        leaf = tree & str(dict_colors_leaves[key][0])
        color_leaves(leaf, key)

ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180
ts.arc_span = 180

tree.show(tree_style=ts)
