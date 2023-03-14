#!/usr/bin/env python

import csv
# Import all the packages necessary to do the analysis
import os
import pickle

import pandas as pd
from Bio import SeqIO

# %% Set of working directory and the files that are going to storage the data run by run
os.chdir(
    '/Users/u1866168/Documents/OneDrive - University of '
    'Warwick/Experiments/GET_HOMOLOGUES/GH_Prokka_run/OMCL_alg/Prokka_OMCL_clusters_fna/')

PickleCluster = open("/Users/u1866168/Documents/OneDrive - University of "
                     "Warwick/Experiments/Experiments_Python/Dict_Cluster.pickle", "wb")
pickle_dictCreation = {}
pickle.dump(pickle_dictCreation, PickleCluster)
PickleCluster.close()  # Pickle will store the data as Dictionary class. Input for the last part of the script.

for clusters in os.listdir('/Users/u1866168/Documents/OneDrive - University of '
                           'Warwick/Experiments/GET_HOMOLOGUES/GH_Prokka_run/OMCL_alg/Prokka_OMCL_clusters_fna/'):
    if not clusters.startswith('.') and clusters.endswith('.fna'):

        with open("/Users/u1866168/Documents/OneDrive - University of "
                  "Warwick/Experiments/Experiments_Python/Dict_Cluster.pickle", "rb") as output2:
            pickleDict = pickle.load(output2)
        cluster = str(clusters).replace('.fna', '')  # storage of the variable name and the code for the cluster
        print(cluster)
        cluster2 = cluster.replace("'", "")  # this step will remove ' from the cluster's name, later will crash if not.
        Array1 = [seq_record.id[3:] for seq_record in SeqIO.parse(clusters, "fasta")]  # Extraction of ID as
        # dictionary keys
        Array2 = [cluster2] * (len(Array1))  # Cluster associated to ID dictionary
        Final_dict = {}

        for i in range(len(Array1)):  # Construction of the dictionary and updating of the files
            Addition = {Array1[i]: Array2[i]}
            Final_dict.update(Addition)
        pickleDict.update(Final_dict)

        with open("/Users/u1866168/Documents/OneDrive - University of "
                  "Warwick/Experiments/Experiments_Python/Dict_Cluster.pickle", "wb") as wfp:
            pickle.dump(pickleDict, wfp)

# %% Establishment of new files which will store the second part of the code data
filecreation = open("/Users/u1866168/Documents/OneDrive - University of "
                    "Warwick/Experiments/Experiments_Python/Dict_LocusProtID.pickle", "wb")
pickle_intialdict = {}
pickle.dump(pickle_intialdict, filecreation)
filecreation.close()

GeneCsv = csv.writer(
    open("/Users/u1866168/Documents/OneDrive - University of "
         "Warwick/Experiments/Experiments_Python/FullListGenesCyanobacteria.csv", "w"))
# CSV file will store the data in a previous state to validation

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'Genome_Database/FinalSelectionGenomes/Prokka_selection/')

# Starting point for the iterator through the files in directory
for files in os.listdir(
        '/Users/u1866168/Documents/OneDrive - University of Warwick/'
        'Genome_Database/FinalSelectionGenomes/Prokka_selection/'):
    if not files.startswith('.'):
        # Setting of variables that will contain the CDS data every run and opening of the files for appending that data
        with open("/Users/u1866168/Documents/OneDrive - University of "
                  "Warwick/Experiments/Experiments_Python/Dict_LocusProtID.pickle", "rb") as output2:
            pickleDict = pickle.load(output2)
        gbank = SeqIO.parse(files, "genbank")  # Genbank parsing
        strain = []
        cluster = {}
        cyanorakclusterid = []
        CDS = []
        product = []
        aaseq = []
        ntseq = []
        genus = []
        if files.endswith('.gbk'):
            for genome in gbank:
                strain = files.split('-', 1)[1].split('.', 1)[0]
                if "Pro" in strain:
                    genus = str("Prochlorococcus")
                else:
                    genus = str("Synechococcus")
                for feature in genome.features:  # Due to the diverse gene annotation, try statements are required.
                    if feature.type == "CDS":
                        ID = str(feature.qualifiers["locus_tag"][0]).replace(' ', '')
                        try:
                            product = str(feature.qualifiers["product"][0])
                        except:
                            product = str("Product_no_annotated")
                        ntseq = str(feature.extract(genome.seq))
                        ntseq2 = feature.extract(genome.seq)
                        aaseq = str(ntseq2.translate(table=11, cds=False, stop_symbol="#", gap=None))

                        CDS = [product, ntseq, aaseq, genus, strain]
                        Addition = {ID: CDS}
                        AdditionList = [ID, product, ntseq, aaseq, genus, strain]
                        cluster.update(Addition)
                        GeneCsv.writerow(AdditionList)
            Final = str(cluster)
            pickleDict.update(cluster)
            with open("/Users/u1866168/Documents/OneDrive - University of "
                      "Warwick/Experiments/Experiments_Python/Dict_LocusProtID.pickle", "wb") as wfp:
                pickle.dump(pickleDict, wfp)

# %% Starting the third part of the code where it will match the locus_tag/prot_id and cyanorak and Get_homologue clusters

os.chdir("/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Experiments_Python/")

FinalOutput = csv.writer(open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                              "Experiments/Experiments_Python/GeneKeyTable.csv", "w"))
# Database with the final outcome for the genkey table.
ClusterGHtable = csv.writer(
    open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
         "Experiments/Experiments_Python/ClusterGHTable.csv", "w"))
# Creation of the cluster table for GH data


Dictionary1 = open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                   "Experiments/Experiments_Python/Dict_Cluster.pickle", "rb")
Clusters_GH = pickle.load(Dictionary1)  # Input from the first part of the script
Dictionary2 = open("/Users/u1866168/Documents/OneDrive - University of Warwick/"
                   "Experiments/Experiments_Python/Dict_LocusProtID.pickle", "rb")
Database = pickle.load(Dictionary2)  # Input from the second part of the script

Keys1 = set(list(Clusters_GH.keys()))
Keys2 = set(list(Database.keys()))

Intersec = list(set(Keys1) & set(Keys2))

# Data from Get_Homologues to fill the GH table.

Occupancy = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                        'Warwick/Experiments/Experiments_Python/pangenome_matrix_occupancy_corrected.csv',
                        index_col='Gene')

DictionaryOccupancy = dict(zip(Occupancy.index, Occupancy['Occupancy']))

clusterfraction = {}
for key in DictionaryOccupancy:
    fractiontemp = []
    if DictionaryOccupancy[key] == 110:
        fractiontemp = "Core"
    elif DictionaryOccupancy[key] <= 2:
        fractiontemp = "Cloud"
    elif DictionaryOccupancy[key] >= 104:
        fractiontemp = "Soft_Core"
    else:
        fractiontemp = "Shell"
    clusterfraction.update({key: [DictionaryOccupancy[key], fractiontemp]})

with open('/Users/u1866168/Documents/OneDrive - University of Warwick/'
          'Experiments/Experiments_Python/Dict_ClusterFraction', "wb") as PickleCluster:
    pickle.dump(clusterfraction, PickleCluster)

# Assembling of data from Get_Homologues and data obtained from genbank files
for i in Intersec:  # Creation of Genkey table and CKTable
    Completed_dataset1 = [i, Database[i][0], Clusters_GH[i], Database[i][1], Database[i][2], Database[i][3],
                          Database[i][4]]
    Completed_dataset4 = [i, Database[i][0], Clusters_GH[i], clusterfraction[Clusters_GH[i]][1]
        , clusterfraction[Clusters_GH[i]][0]]
    FinalOutput.writerow(Completed_dataset1)
    ClusterGHtable.writerow(Completed_dataset4)
