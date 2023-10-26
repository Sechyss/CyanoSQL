#!/usr/bin/env python

import collections
import csv
# Import all the packages necessary to do the analysis
import os
import pickle

import matplotlib.pyplot as plot
import pandas as pd
from Bio import SeqIO

# Set of working directory and the files that are going to storage the data run by run
os.chdir(
    '/Users/u1866168/Documents/OneDrive - University of '
    'Warwick/Experiments/GET_HOMOLOGUES/FinalDraftGetHomologues/OMCL_Clusters/Unknown_fraction')

PickleCluster = open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/Dict_ClusterUnk.pickle", "wb")
pickle_dictCreation = {}
pickle.dump(pickle_dictCreation, PickleCluster)
PickleCluster.close()  # Pickle will store the data as Dictionary class. Input for the last part of the script.

for clusters in os.listdir(
        '/Users/u1866168/Documents/OneDrive - University of '
        'Warwick/Experiments/GET_HOMOLOGUES/FinalDraftGetHomologues/OMCL_Clusters/Unknown_fraction'):
    if not clusters.startswith('.'):

        with open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/Dict_ClusterUnk.pickle",
                  "rb") as output2:
            pickleDict = pickle.load(output2)
        cluster = str(clusters)  # storage of the variable name and the code for the cluster
        cluster2 = cluster.replace("'", "")  # this step will remove ' from the cluster's name, later will crash if not.
        Array1 = [seq_record.id[3:] for seq_record in SeqIO.parse(clusters, "fasta")]  # Extraction of ID as
        # dictionary keys
        Array2 = [cluster2] * (len(Array1))  # Cluster associated to ID dictionary
        Final_dict = {}

        for i in range(len(Array1)):  # Construction of the dictionary and updating of the files
            Addition = {Array1[i]: Array2[i]}
            Final_dict.update(Addition)
        pickleDict.update(Final_dict)

        with open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/Dict_ClusterUnk.pickle", "wb") as wfp:
            pickle.dump(pickleDict, wfp)

os.chdir('/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results')

# Creation of the cluster table for Cyanorak
ClusterUnktable = csv.writer(
    open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/ClusterUnkTable.csv", "w"))
Dictionary1 = open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/Dict_ClusterUnk.pickle", "rb")
Clusters_Unk = pickle.load(Dictionary1)  # Input from the first part of the script
Dictionary2 = open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/Dict_LocusProtID.pickle", "rb")
Database = pickle.load(Dictionary2)  # Input from the second part of the script

Unknown_Dict = {}
Keys1 = set(list(Clusters_Unk.keys()))
Keys2 = set(list(Database.keys()))

Intersec = list(set(Keys1) & set(Keys2))

# Prior to this check the extension in the csv file
Occupancy = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                        'Warwick/Experiments/Experiments_Python/pangenome_matrix_occupancy_corrected.csv',
                        index_col='Gene')

DictionaryOccupancy = dict(zip(Occupancy.index_genome, Occupancy['Occupancy']))

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

for i in Intersec:
    Completed_dataset1 = {
        Database[i][0]: [Clusters_Unk[i], Database[i][1], Database[i][4], Database[i][5],
                         clusterfraction[Clusters_Unk[i]][1], clusterfraction[Clusters_Unk[i]][0]]}
    Completed_dataset2 = [Database[i][0], Clusters_Unk[i], Database[i][1], Database[i][4], Database[i][5],
                          clusterfraction[Clusters_Unk[i]][1], clusterfraction[Clusters_Unk[i]][0]]
    Unknown_Dict.update(Completed_dataset1)
    ClusterUnktable.writerow(Completed_dataset2)

with open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/Dict_ClusterUnk_SQL.pickle",
          "wb") as PickleCluster:
    pickle.dump(Unknown_Dict, PickleCluster)

Fraction = []
ClusterList = []
for i in Unknown_Dict.keys():
    if Unknown_Dict[i][0] not in ClusterList:
        Fraction.append(Unknown_Dict[i][4])
        ClusterList.append(Unknown_Dict[i][0])

FractionDiv = collections.Counter(Fraction)
FractionDiv["Soft_Core"] = (FractionDiv["Soft_Core"] + FractionDiv[
    "Core"])  # The sum of the soft core plus the core (that should be included based on its properties)
colors = ["yellow", "darkorange", "red", "black"]
Dictionary = plot.bar(FractionDiv.keys(), FractionDiv.values(), color=colors)
plot.xlabel("Fraction")
plot.ylabel("Number of clusters")
plot.show()
