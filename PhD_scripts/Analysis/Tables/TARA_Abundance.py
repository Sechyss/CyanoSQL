# %% Load the tables from TARA and collection of relevant data
import pandas as pd


abundancedf = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                          'Warwick/Experiments/Analysis_Genes/Psip1/Psip1_phylo/'
                          'OceanGeneAtlas_Data_20220121-17_16_Hmm_trimcurated_homologs_MetaT/abundance_matrix.csv',
                          sep=',', header=0, index_col=0, skiprows=[1])
df_nutrient = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/Psip1/'
                          'Psip1_phylo/OceanGeneAtlas_Data_20220121-17_16_Hmm_trimcurated_homologs_MetaT'
                          '/environmental_parameters.csv',
                          sep='\t')

taramedians = pd.read_excel('/Users/u1866168/Documents/OneDrive - University of '
                            'Warwick/Experiments/Analysis_Genes/Psip1/Psip1_phylo/Taramedians.xlsx',
                            sheet_name='Metatranscriptome', header=0, index_col=0)

# In case you want to select a particular set of microorganisms, change the abundancedf in the next part of the code
abundance_cyano = abundancedf[abundancedf['taxonomy'].str.contains('Cyanobacteria')]

# Creation of the dataframe with all the data for the abundance
dict_abundance = {}
for row in abundancedf.columns[1:]:
    if 'TARA_' in row:
        dict_abundance.update({row: [abundancedf[row][1:].sum(), str(row).rsplit('_')[2]]})

#   ---- Addition of the raw data of metagenome, location, median from Andrews unique genes and the normalised data
df_tocsv = pd.DataFrame.from_dict(dict_abundance, orient='index', columns=['Psip1', 'Location'])
df_tocsv = pd.concat([df_tocsv, taramedians], axis=1, join='inner')  # Addition of the normalized data
df_tocsv['Normalised_data(%)'] = (df_tocsv['Psip1']/df_tocsv['Median'])*100

#   ---- Finish off by adding the extra features to the table to the main matrix
df_tocsv['Region'] = pd.Series(dict(zip(list(df_nutrient['OGA_ID.1']), list(df_nutrient['Region']))))
df_tocsv['Region'] = list(str(x).rsplit(' ')[0].replace('[', '').replace(']', '') for x in df_tocsv['Region'])
df_tocsv['Province'] = pd.Series(dict(zip(list(df_nutrient['OGA_ID.1']), list(df_nutrient['Province']))))
df_tocsv['Province'] = list(str(x).rsplit(' ')[0].replace('[', '').replace(']', '') for x in df_tocsv['Province'])
df_tocsv['latitude'] = pd.Series(dict(zip(list(df_nutrient['OGA_ID.1']), list(df_nutrient['latitude']))))
df_tocsv['longitude'] = pd.Series(dict(zip(list(df_nutrient['OGA_ID.1']), list(df_nutrient['longitude']))))
df_tocsv['PO4 (µmol/l)'] = pd.Series(dict(zip(list(df_nutrient['OGA_ID.1']), list(df_nutrient['PO4 (µmol/l)']))))
df_tocsv['Iron_5m* (µmol/l)'] = pd.Series(dict(zip(list(df_nutrient['OGA_ID.1']),
                                                   list(df_nutrient['Iron_5m* (µmol/l)']))))

#  ---- Write the final output into a file with the final pivot table

df_tocsv.to_csv('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/Psip1/'
                'Psip1_phylo/OceanGeneAtlas_Data_20220121-17_16_Hmm_trimcurated_homologs_MetaT/'
                'Station_Gene_MetaTData_Cyano.csv')
