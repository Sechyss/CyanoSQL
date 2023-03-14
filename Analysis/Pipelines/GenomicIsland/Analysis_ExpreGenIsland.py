"""
Exploratory statistical analysis for the Genomic island genes using its AMT data
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

df_GI = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                    'MHC_Scripts/SQL_Tables/SQL_queries/GenomicIslandGenes_AMT.csv', index_col=0, low_memory=True)

stations = list(df_GI.columns)
colors = ['DarkGoldenRod', 'Coral', 'DarkGreen', 'FireBrick', 'LightSteelBlue', 'Navy', 'Olive', 'Purple',
          'Wheat', 'Gold', 'DarkSalmon', 'Chocolate', 'CadetBlue', 'BlueViolet', 'Cyan', 'GreenYellow']
# %% --- Where are the genes in genomic island mostly expressed?


values = []
pct = []

for i in stations:
    values.append(df_GI[i].sum())
    pct.append(str(round((df_GI[i].sum() * 100 / df_GI.values.sum()), 2)))

plt.figure(figsize=[10, 7])
plt.pie(values, labels=pct, colors=colors)
plt.suptitle('Percentage of TPM from genomic island genes per station', fontsize=12, fontweight="bold")
plt.legend(stations, loc='upper right', bbox_to_anchor=(1.325, 1), fontsize='small')
plt.show()

# %% --- PCA analysis

x = df_GI.values
x = StandardScaler().fit_transform(x)
np.mean(x), np.std(x)

feat_cols = ['station ' + str(i) for i in stations]
normalised_tpm = pd.DataFrame(x, columns=feat_cols)

pca = PCA(n_components=2)
principalComponents_AMT = pca.fit_transform(x)
principal_AMT_Df = pd.DataFrame(data=principalComponents_AMT
                                , columns=['PC1', 'PC2'])

plt.figure()
plt.figure(figsize=(10, 10))
plt.xticks(fontsize=12)
plt.yticks(fontsize=14)
plt.xlabel('Principal Component - 1', fontsize=20)
plt.ylabel('Principal Component - 2', fontsize=20)
plt.title("Principal Component Analysis of AMT-GI p/ Station", fontsize=20)


