import pandas as pd
from scipy.stats import mannwhitneyu, kruskal

metaG = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/Psip1/'
                    'Psip1_phylo/OceanGeneAtlas_20211217-09_49_Hmm_trimcurated_homologs/cyanopsipphoxmg.csv')
metaT = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/Psip1/'
                    'Psip1_phylo/OceanGeneAtlas_Data_20220121-17_16_Hmm_trimcurated_homologs_MetaT/cyanopsipphoxmt.csv')
zonesG = list(set(metaG['Region']))
zonesT = list(set(metaT['Region']))
proteins = list(set(metaG['Gene']))


# Extraction of information for any table
def extraction_data_df(genes, regions, df):
    empty_dict = {}
    for gene in genes:
        empty_dict[gene] = {}
        for region in regions:
            abundance = list(df[(df['Region'] == region) & (df['Gene'] == gene)]['Abundance'])
            empty_dict[gene][region] = abundance

    return empty_dict


info = extraction_data_df(proteins, zonesG, metaG)

df_collector = pd.DataFrame(columns=['Abundance', 'Tag', 'Region', 'Protein'])
for i in info.keys():
    for y in info[i]:
        newdf = pd.DataFrame({'Tag': str(i) + '_' + str(y), 'Abundance': info[i][y], 'Region': str(y),
                              'Protein': str(i)})
        df_collector = pd.concat([df_collector, newdf], axis=0)


kruskal(*[group["Abundance"].values for name, group in df_collector.groupby("Tag")])

for i in zonesT:
    psip1 = info['psip1'][i]
    phoX = info['phoX'][i]

    try:
        U1, p = mannwhitneyu(psip1, phoX, method="exact")
        print(i, p)
    finally:
        continue
