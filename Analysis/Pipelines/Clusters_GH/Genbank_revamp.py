import glob
import os

import pandas as pd
from MyPackage.GenbankEditing import GenbankEditing

# %%


os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cya_clusterPipeline')

# This code will read through the genbank files with manually curated genes to assign a new locustag number to those
# genes
for files in glob.glob("/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Temporal/*"):
    if not files.startswith('.') and files.endswith('.gbk'):
        GenbankEditing(files).genbank_revamp(files.rsplit(".", 1)[0] + ".gb")

# %% Example of creation of gene table from SynWH7803

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Temporal')

genbank = SeqIO.parse('23615-Syn_WH7803.gb', 'genbank')

Dict = {'Locustag': ['Product', 'Feature']}

for genome in genbank:
    for feature in genome.features:
        if feature.type != 'source':
            Locustag = str(feature.qualifiers["locus_tag"][0])
            product = str(feature.qualifiers["product"][0])
            Dict.update({Locustag: [product, feature.type]})
Dataframe = pd.DataFrame.from_dict(Dict, orient='index_genome')
Dataframe.to_csv('SynWH7803_genes.csv')
