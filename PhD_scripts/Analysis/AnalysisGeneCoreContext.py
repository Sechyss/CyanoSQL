import collections
import csv
import os
import pickle

import matplotlib.pyplot as plot

os.chdir('/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results')
# This will read the dictionary created during the extraction and translation into clusters of the different locustag
Dictionary1 = open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/Dict_ClusterAnalysis.pickle", "rb")
Cluster_Unk = pickle.load(Dictionary1)
# This part will open the CSV created by SQL where contains the Clade3 locustag
with open('/Users/u1866168/Documents/OneDrive - University of '
          'Warwick/Experiments/Experiments_Python/GeneContext_CoreUnk.csv', newline='') as genecontext:
    Context = list(csv.reader(genecontext))

ClusterUnkContext = csv.writer(
    open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/ClusterUnkContext.csv", "w"))

ClustersKeys = Cluster_Unk.keys()

for i in range(len(Context)):
    if i == 0:  # This is to avoid errors when reading the file, first row is the header
        pass
    else:
        ContextLocus = [Context[i][0]]
        for y in Context[i]:
            if y in ClustersKeys:
                ContextLocus.append(Cluster_Unk[y][0])  # Translation of the locustag into clusters and addition to list
                pass
        print(ContextLocus)
        ClusterUnkContext.writerow(ContextLocus)

with open('/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/ClusterUnkContext.csv',
          newline='') as genecontext:
    ContextLocus = list(csv.reader(genecontext))
# Reading of the csv already created
UnkContextPattern = csv.writer(
    open("/Users/u1866168/PycharmProjects/Cyanobacteria/Python_results/ClusterUnkContextPatterns.csv", "w"))
# Creation of the pattern based on surrounding clusters into a numeric code
ListCluster = []
Metadata_cluster = []
for i in range(len(ContextLocus)):
    Pattern = ContextLocus[i][1:]
    if Pattern not in ListCluster:
        ListCluster.append(Pattern)
        UnkContextPattern.writerow(Pattern)
        Metadata_cluster.append(ContextLocus[i][1])
# Count of the number of different patterns per cluster
Metadata_cluster = collections.Counter(Metadata_cluster)
Distribution = Metadata_cluster.values()

# Plot of the graph
Dictionary = plot.figure()
plot.hist(Distribution)
plot.xlabel('Number of patterns/cluster')
plot.ylabel('Number of cluster')
plot.title('Distribution of the number of patterns')
plot.show()
