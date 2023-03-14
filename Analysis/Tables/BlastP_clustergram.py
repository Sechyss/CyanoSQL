import math

import pandas as pd
import scipy.cluster.hierarchy as shc
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

blastp_output = pd.read_table('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                              'Thesis/Paper Draft/Table_distance_Blastp_APs_2.tsv')
Df_Ap = pd.read_excel('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                      'Thesis/Paper Draft/Dict_ID_AlkPh.xlsx', sheet_name='Sheet1')
ids = list(str(x) for x in Df_Ap['ID'].values)
Dict_AP = pd.Series(Df_Ap['AlkP'].values, index=ids).to_dict()
collector_df = pd.DataFrame(index=list(set(blastp_output['Taxa1'])), columns=list(set(blastp_output['Taxa1'])))

for index, row in blastp_output.iterrows():
    collector_df[blastp_output.loc[index]['Taxa1']][blastp_output.loc[index]['Taxa2']] = \
        math.log(1+ float(blastp_output.loc[index]['evalue']))

collector_df = collector_df.fillna(5)


lut = dict(zip(['Psip1', 'PhoX_Flavo', 'PhoX', 'PhoA', 'PhoD', 'PafA'],
               ['#FF6A6A', '#8968CD', '#8B3A62', '#6E8B3D', '#008B8B', '#FFB90F']))
phosphatases = pd.Series(Dict_AP).map(lut)
dend_reordered = shc.linkage(collector_df, method='ward', optimal_ordering=True)
clustering_dend = sns.clustermap(collector_df, metric="euclidean", method="ward", row_colors=phosphatases,
                                 col_colors=phosphatases,
                                 row_linkage=dend_reordered, col_linkage=dend_reordered,
                                 cmap="YlGnBu")
handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut, title='AlkPhos',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

plt.savefig('/Users/u1866168/Documents/OneDrive - University of Warwick/'
            'Thesis/Paper Draft/Clustergram_blastp-evalue.pdf')
plt.show()
