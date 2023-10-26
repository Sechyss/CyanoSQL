import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

dfabundance = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/Psip1/'
                          'Psip1_phylo/OceanGeneAtlas_20211213-13_56_Psip1_hmm_homologs/abundance_matrix.tsv', sep='\t')

df_nutrient = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/Psip1/'
                          'Psip1_phylo/OceanGeneAtlas_20211213-13_56_Psip1_hmm_homologs/environmental_parameters.tsv',
                          sep='\t')

dict_nutrient = dict(zip(list(df_nutrient['sample_ID']), list(df_nutrient['PO4 (µmol/l)'])))
dict_abundance = {}

for row in dfabundance.columns:
    if 'TARA_' in row:
        dict_abundance.update({row:dfabundance[row].sum()})

x = []
y = []
for station in dict_nutrient.keys():
    if dict_nutrient[station] != 'nan':
        x.append(dict_nutrient[station])
        y.append(dict_abundance[station])
    else:
        continue
x = np.array(x)
y = np.array(y)
idx = np.isfinite(x) & np.isfinite(y)

plt.figure(figsize=(8,8))
plt.scatter(x[idx], y[idx])
plt.ylabel('Total abundance station')
plt.xlabel('Log PO4 (µmol/l)')