# Import necessary packages

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy import stats

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/Phylogenetic_congruency/')
# %% Upload of the different datasets from SH_table.py

shellcore = pd.read_csv('Intergenic_recombination_shell_softcore.csv', index_col=0)
strictcore = pd.read_csv('Intergenic_recombination_core.csv', index_col=0)

mantel_shellcore = np.array(shellcore['mantel_value'])
mantel_core = np.array(strictcore['mantel_value'])

rbfoulds_shellcore = np.array(shellcore['Percentage_dev'])
rbfoulds_core = np.array(strictcore['Percentage_dev'])

ticks = ['Core', 'SoftCoreShell']

test_mantel = [mantel_core, mantel_shellcore]
test_rbfoulds = [rbfoulds_core, rbfoulds_shellcore]

# %% Extraction Soft-Core vs Shell

mother_table = pd.read_excel('Intergenic_recombination_all.xlsx', sheet_name='Both_tests', index_col=0)

softcore_mantel = np.array(mother_table[mother_table['Shared taxa'] > 103]['mantel_value'])
shell_mantel = np.array(mother_table[mother_table['Shared taxa'] < 104]['mantel_value'])

softcore_rbfoulds = np.array(mother_table[mother_table['Shared taxa'] > 103]['Percentage_dev'])
shell_rbfoulds = np.array(mother_table[mother_table['Shared taxa'] < 104]['Percentage_dev'])

test_mantel_2 = [softcore_mantel, shell_mantel]
test_rbfoulds_2 = [softcore_rbfoulds, shell_rbfoulds]

genes_softcore = np.array(mother_table[mother_table['Shared taxa'] > 103].index)
genes_shell = np.array(mother_table[mother_table['Shared taxa'] < 104].index)
# %% Histogram

Xarray = np.array(mother_table['Shared taxa'])
Yarray = np.array(mother_table['Percentage_dev'])
Y2array = np.array(mother_table['mantel_value'])

# %%  Boxplot creation

# --- Labels for your data:
labels_list = ['Mantel distance                 ',
               'Robinson-Foulds distance             ']
xlocations = range(len(test_mantel))
width = 0.3
symbol = 'r+'
ymin = 0
ymax = 1.2

ax = plt.gca()
ax.set_ylim(ymin, ymax)
ax.set_xticklabels(labels_list, rotation=0)
ax.grid(True, linestyle='dotted')
ax.set_axisbelow(True)
ax.set_xticks(xlocations)
plt.xlabel('Intergenic tests')
plt.ylabel('Correlation')

colors = ['bisque', 'lightblue', 'steelblue', 'darkorange']
colorcore = dict(color=colors[3])
colorshell = dict(color=colors[2])
colorshellmedian = dict(color=colors[2])
colorcoremedian = dict(color=colors[3])

# plt.title('Intergenic tests')

# --- Offset the positions per group:
positions_group1 = [x - (width + 0.01) for x in xlocations]
positions_group2 = xlocations

plt.boxplot(test_mantel,
            sym=symbol,
            labels=[''] * len(labels_list),
            positions=positions_group1,
            widths=width,
            boxprops=colorcore,
            medianprops=colorcoremedian)

plt.boxplot(test_rbfoulds,
            labels=labels_list,
            sym=symbol,
            positions=positions_group2,
            widths=width,
            boxprops=colorshell,
            medianprops=colorshellmedian)

plt.legend(['SoftCore', 'Shell'],
           loc='best', fontsize='xx-small')

plt.savefig('Boxplot_intergenic_analysis_core_vs_shell-softcore.pdf')
plt.show()

# %% Test normality of the data

kolmogorov_1, pvalue_k_1 = stats.kstest(softcore_mantel, 'norm')
kolmogorov_2, pvalue_k_2 = stats.kstest(shell_mantel, 'norm')
kolmogorov_3, pvalue_k_3 = stats.kstest(softcore_rbfoulds, 'norm')
kolmogorov_4, pvalue_k_4 = stats.kstest(shell_rbfoulds, 'norm')

column_headers = ['Kolmogorov_Smirnov', 'pvalue_k', 'Levene', 'pvalue_l', '1-way_Anova', 'pvalue']

levene_m, pvalue_lm = stats.levene(softcore_mantel, shell_mantel)
levene_rb, pvalue_rb = stats.levene(softcore_rbfoulds, shell_rbfoulds)

anova_m, pvalue_am = stats.f_oneway(softcore_mantel, shell_mantel)
anova_rb, pvalue_arb = stats.f_oneway(softcore_rbfoulds, shell_rbfoulds)

dictionary = {'Levene': [levene_m, levene_rb],
              'pvalue_l': [pvalue_lm, pvalue_rb],
              '1-way_Anova': [anova_m, anova_rb],
              'pvalue': [pvalue_am, pvalue_arb]}

ANOVA_intergenic_summary = pd.DataFrame(dictionary, index=['Mantel_test', 'Robinson_Foulds'])

# %% Analysis of bins to extract the 0 to 0.2 fraction

allvalues = np.array(mother_table['mantel_value'])
histogramall = np.histogram(allvalues, bins=5)

print(histogramall)

genes_bins = list(mother_table[mother_table['mantel_value'] < 0.27325643].index)

with open('./Intergenic_statistics/List_genes_lowest_mantelv.txt', 'w') as file:
    for i in genes_bins:
        file.write("%s\n" % i.split('_', 1)[1].replace('_', ' '))

Table_clusters = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                             'Warwick/Genome_Database/Metadata/Shell_SoftCore_phylo_Prokka.csv', index_col="id")

redundant_clusters = []
for index, row in Table_clusters.iterrows():
    if row['ClusterGH'] in genes_bins and (row['ClusterGH'] not in redundant_clusters):
        with open('List_genes_lowest_mantelv.fasta', 'a') as outputfile:
            aaseq = Seq(row['AASEQ'].replace('#', ''))
            fastainput = SeqRecord(aaseq, id=str(row['ClusterGH']), description='')
            SeqIO.write(fastainput, outputfile, "fasta")
            redundant_clusters.append(row['ClusterGH'])

for index, row in Table_clusters.iterrows():
    if row['ClusterGH'] in genes_bins and ('8102' in row['Strain']):
        with open('ListLocustag_genes_lowest_mantelv.txt', 'a') as outputfile:
            outputfile.write("%s\n" % index)
