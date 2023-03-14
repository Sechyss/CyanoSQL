# %%

import csv
import os
import pickle

import pandas as pd
from Bio import SeqIO

# %%  Establishment of new files which will store the second part of the code data
GeneCsv = csv.writer(
    open("/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/SQL_Scripts/Genecontext_Prokka.csv", "w"))
# CSV file will store the data in a previous state to validation

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'Genome_Database/FinalSelectionGenomes/Prokka_selection/')

# Starting point for the iterator through the files in directory
for files in os.listdir():
    if not files.startswith('.') and files.endswith('gbk'):
        print(files)
        # Setting of variables that will contain the CDS data every run and opening of the files for appending that data
        gbank = SeqIO.parse(files, "genbank")  # Genbank parsing.
        Locustags = {}
        for genome in gbank:
            for feature in genome.features:  # Due to the diverse gene annotation, try statements are required.
                if feature.type == "CDS":  # Apply to those that are not manually created and might have CyaClusterID
                    ID = str(feature.qualifiers["locus_tag"][0]).replace(' ', '')
                    if feature.location.strand == -1:
                        Locustags.update({ID: "c"})
                    else:
                        Locustags.update({ID: "s"})
        x = list(Locustags.keys())
        for i in range(len(x)):
            if i == 0:
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + x[(i + 5):(i):-1]
                    Downstream = ["None"] * 5
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + ["None"] * 5
                    Downstream = x[(i + 1):(i + 6)]
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == 1:
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + x[(i + 5):(i):-1]
                    Downstream = [x[(i - 1)]] + ["None"] * 4
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + ["None"] * 4 + [x[(i - 1)]]
                    Downstream = x[(i + 1):(i + 6)]
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == 2:
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + x[(i + 5):i:-1]
                    Downstream = x[:i] + ["None"] * 3
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + ["None"] * 3 + x[:i]
                    Downstream = x[(i + 1):(i + 6)]
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == 3:
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + x[(i + 5):i:-1]
                    Downstream = x[:i] + ["None"] * 2
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + ["None"] * 2 + x[:i]
                    Downstream = x[(i + 1):(i + 6)]
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == 4:
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + x[(i + 5):(i):-1]
                    Downstream = x[:i] + ["None"]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + ["None"] + x[:i]
                    Downstream = x[(i + 1):(i + 6)]
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == 5:
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + x[(i + 5):(i):-1]
                    Downstream = x[:i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = x[(i + 1):(i + 6)]
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == (len(x) - 5):
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + ["None"] + x[len(x):i:-1]
                    Downstream = x[(i - 5):i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = x[(i + 1):] + ["None"]
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == (len(x) - 4):
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + ["None"] * 2 + x[len(x):i:-1]
                    Downstream = x[(i - 5):i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = x[(i + 1):] + ["None"] * 2
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == (len(x) - 3):
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + ["None"] * 3 + x[len(x):i:-1]
                    Downstream = x[(i - 5):i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = x[(i + 1):] + ["None"] * 3
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == (len(x) - 2):
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + ["None"] * 4 + x[len(x):i:-1]
                    Downstream = x[(i - 5):i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = x[(i + 1):] + ["None"] * 4
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == (len(x) - 1):
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + ["None"] * 5 + x[len(x):i:-1]
                    Downstream = x[(i - 5):i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = x[(i + 1):] + ["None"] * 5
                    GeneCsv.writerow(Upstream + Downstream)
            elif i == len(x):
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + ["None"]
                    Downstream = x[(i - 5):i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = ["None"] * 5
                    GeneCsv.writerow(Upstream + Downstream)
            else:
                if Locustags[x[i]] == "s":
                    Upstream = [x[i]] + x[(i + 5):i:-1]
                    Downstream = x[(i - 5):i]
                    GeneCsv.writerow(Upstream + Downstream)
                else:
                    Upstream = [x[i]] + x[(i - 5):i]
                    Downstream = x[(i + 1):(i + 6)]
                    GeneCsv.writerow(Upstream + Downstream)

# %% Transformation from locustag to clusters

Context = pd.read_csv("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                      "MHC_Scripts/SQL_Scripts/Genecontext_Prokka.csv", header=None, index_col=0)
Dictionary1 = open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                   "Experiments/Experiments_Python/Dict_Cluster.pickle", "rb")
Clusters_GH = pickle.load(Dictionary1)
Clusters_GH.update({'None': 'None'})
Dictionary_collector = {}
for index, row in Context.iterrows():
    key = index
    value = [Clusters_GH[row[1]], Clusters_GH[row[2]], Clusters_GH[row[3]], Clusters_GH[row[4]], Clusters_GH[row[5]],
             Clusters_GH[row[6]], Clusters_GH[row[7]], Clusters_GH[row[8]], Clusters_GH[row[9]], Clusters_GH[row[10]]]
    Dictionary_collector.update({key: value})

Header = ['Context+5', 'Context+4', 'Context+3', 'Context+2', 'Context+1',
          'Context-1', 'Context-2', 'Context-3', 'Context-4', 'Context-5']
Dict_clust = pd.DataFrame.from_dict(Dictionary_collector, orient='index', columns=Header)
Dict_clust.to_csv("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                  "MHC_Scripts/SQL_Scripts/Genecontext_Prokka_Clusters.csv")
