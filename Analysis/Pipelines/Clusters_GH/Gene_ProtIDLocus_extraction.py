#!/usr/bin/env python

import pickle
import random
import sys

from Bio import SeqIO

gbank = SeqIO.parse(open(sys.argv[1], "rU"), "genbank")  # Genbank parsing.
output_handle = open("LocusProtID.tab", "a")
with open("Dict_LocusProtID.pickle", "rb") as output2:
    pickleDict = pickle.load(output2)
strain = []
cluster = {}
cyanorakclusterid = []
CDS = []
product = []
aaseq = []
ntseq = []  # Creating a array with all the features associated to each rRNA.
for genome in gbank:
    for feature in genome.features:
        if feature.type == "source":
            strain = str(feature.qualifiers["organism"][0])

        if feature.type == "CDS":
            try:
                ID = str(feature.qualifiers["protein_id"][0])
                product = str(feature.qualifiers["product"][0])
                ntseq = feature.extract(genome.seq)
                aaseq = ntseq.translate(table="Bacterial", cds=True)
                CDS = [strain, product, ntseq, aaseq]
                Addition = {ID: CDS}
                cluster.update(Addition)
                print(ID)
            except:
                try:
                    if feature.qualifiers["note"][0] == "cyanorak_state: manually created":
                        ID = str(str(cluster.key(1))[0:4] + "_m" + str(random.randint(1, 10001)))
                        # ID = str(feature.type["LOCUS"][0]) + "_" + str(random(1, 100001))
                        product = str(feature.qualifiers["product"][0])
                        ntseq = feature.extract(genome.seq)
                        aaseq = ntseq.translate(table="Bacterial", cds=True)
                        CDS = [strain, product, ntseq, aaseq]
                        Addition = {ID: CDS}
                        cluster.update(Addition)
                        print(ID)

                    else:
                        ID = str(feature.qualifiers["locus_tag"][0])
                        product = str("Non-annotated")
                        ntseq = feature.extract(genome.seq)
                        aaseq = ntseq.translate(table="Bacterial", cds=True)
                        CDS = [strain, product, ntseq, aaseq]
                        Addition = {ID: CDS}
                        cluster.update(Addition)
                        print(ID)
                except:
                    ID = str(feature.qualifiers["locus_tag"][0])
                    product = str("Non-annotated")
                    ntseq = feature.extract(genome.seq)
                    aaseq = ntseq.translate(table="Bacterial", cds=True)
                    CDS = [strain, product, ntseq, aaseq]
                    Addition = {ID: CDS}
                    cluster.update(Addition)
                    print(ID)
Final = str(cluster)
pickleDict.update(cluster)

output_handle.write(Final)  # Write the programme to compile everything in one single file.
output_handle.close()

with open("Dict_LocusProtID.pickle", "wb") as wfp:
    pickle.dump(pickleDict, wfp)
