import pandas as pd

from CyanoScripts.MyPackage.FastaEditing import extract_common_sequences
from CyanoScripts.MyPackage.FileMethods import list_from_file

abundancedf = pd.read_csv('/Users/u1866168/Downloads/OceanGeneAtlas_Data_20211214-13_36_Psip1_hmm/abundance_matrix.csv',
                          sep='\t', header=None, index_col=0)
collector = {}
for index, row in abundancedf.iterrows():
    if 'OM' not in str(index):
        continue
    else:
        collector.update({index: str(row[1]).rsplit(';')[-1]})

test = list_from_file('/Users/u1866168/Documents/OneDrive - University of '
                      'Warwick/Experiments/Analysis_Genes/Psip1/Psip1_phylo/OceanGeneAtlas_20211213'
                      '-13_56_Psip1_hmm_homologs/Node_list.txt')

for i in test:
    try:
        print(i, collector[i], abundancedf.loc[i][1])
    except:
        print(i)

extract_common_sequences('/Users/u1866168/Documents/OneDrive - University of '
                         'Warwick/Experiments/Analysis_Genes/Psip1/Psip1_phylo/BlastP_Psip1_seqs_trimcurated.mafft',
                         '/Users/u1866168/Documents/OneDrive - University of '
                         'Warwick/Experiments/Analysis_Genes/Psip1/Psip1_phylo/Blastp_WH8102_seqs.faa',
                         '/Users/u1866168/Documents/OneDrive - University of '
                         'Warwick/Experiments/Analysis_Genes/Psip1/Psip1_phylo/BLASTP_TARA_Psip1_seqs.faa'
                         )
