#!/usr/bin/env python

import pickle
import sys

from Bio import SeqIO

output = open("Cluster.tab", "a")

with open("Dict_Cluster.pickle", "rb") as output2:
    pickleDict = pickle.load(output2)
cluster = str(sys.argv[1])
print(cluster)
Array1 = [seq_record.id for seq_record in SeqIO.parse(open(sys.argv[1], "rU"), "fasta")]

Array2 = [cluster] * (len(Array1))
Final_dict = {}

for i in range(len(Array1)):
    Addition = {Array1[i]: Array2[i]}
    Final_dict.update(Addition)
pickleDict.update(Final_dict)
Final = str(Final_dict)
output.write(Final)
output.close()

with open("Dict_Cluster.pickle", "wb") as wfp:
    pickle.dump(pickleDict, wfp)
