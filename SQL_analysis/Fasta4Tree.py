"""
Provided a table from the SQL database, this script will extract all sequence per GH cluster, giving the name of
the cluster as ID for the fasta file. It is necessary to provide the directory to store the fasta files.

"""

import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fasta_4_tree(directory_for_files, table_from_SQL):
    os.chdir(directory_for_files)
    list_to_collect = []  # This list will storage the different clusters from the CSV.
    for index, row in table_from_SQL.iterrows():  # Iteration through the table
        if row['ClusterGH'] not in list_to_collect:
            # When there is a new cluster it will create a new field withing the list
            NTSeq = Seq(row['NTSEQ'])  # Creates a Seq type object
            ID = str(row['Strain'])  # Storages the strain name
            # Creates the fastatype record in Nucleotide and Aminoacid level
            FastaNT = SeqRecord(NTSeq[:len(NTSeq) - 2], id=ID.replace(" ", "_"), description="")
            FastaAA = SeqRecord(NTSeq.translate(table=11, stop_symbol=''), id=ID.replace(" ", "_"),
                                description="")
            # Creation of the files that will contain the SeqRecord information
            f1 = open(row['ClusterGH'].rsplit(".", 1)[0] + ".fa", "w")
            f2 = open(row['ClusterGH'].rsplit(".", 1)[0] + ".fasta", "w")
            SeqIO.write(FastaAA, f1, "fasta")
            SeqIO.write(FastaNT, f2, "fasta")

            list_to_collect.append(
                row['ClusterGH'])  # Addition of the cluster to a list so it can follow a different path next time

            f1.close()
            f2.close()
        else:  # In case a file with SeqRecord was created previously this will overcome the potential error
            NTSeq = Seq(row['NTSEQ'])
            ID = str(row['Strain'])
            FastaNT = SeqRecord(NTSeq[:len(NTSeq) - 2], id=ID.replace(" ", "_"), description="")
            FastaAA = SeqRecord(NTSeq.translate(table=11, stop_symbol=''), id=ID.replace(" ", "_"),
                                description="")

            # The process is the same, but instead of overwriting it only adds the new information
            f1 = open(row['ClusterGH'].rsplit(".", 1)[0] + ".fa", "a")
            f2 = open(row['ClusterGH'].rsplit(".", 1)[0] + ".fasta", "a")
            SeqIO.write(FastaAA, f1, "fasta")
            SeqIO.write(FastaNT, f2, "fasta")
            f1.close()
            f2.close()


"""
It is important to mention that additional steps are required downstream:
1- Seqkit
2- Muscle
- TranslatorX

After this pipeline it is necessary to run MUSCLE to align the AA sequence prior to use translatorx:
    $ for file in ./*fa; do seqkit rmdup < $file > ${file%.*}.faa ; done
    $ for file in ./*fasta; do seqkit rmdup < $file > ${file%.*}.fna ; done
    $ for file in *; do res=$(ls -leaf | grep -c '>' $file); if (($res != 110)); then rm $file; fi; done
    $ for file in ./*.faa; do muscle -in $file -out ${file%.*}.muscle; done
    $ for file in ./*.fna; do perl translatorx.pl -i $file -a ${file%.*}.muscle -o ${file%.*}.transx -c 11; done
    $ seqkit concat *nt_ali.fasta > Concatenate_Core_Prokka.fasta

After this remove intermediate files like: 
    $ rm *fa; rm *fasta; rm *muscle

Using the same folder we have storage the TranslatorX result we can create a partition file
"""


def create_partition_file(directory_for_files, outfile):
    counter = 1
    for file in sorted(directory_for_files):
        if file.endswith('.nt_ali.fasta'):
            f2 = SeqIO.parse(file, 'fasta')
            gene = str(file.rsplit('.', 5)[0])  # Storage the name of the cluster
            length = 0
            for seq_record in f2:
                length = len(seq_record.seq) - 1
                break
            outfile.write("DNA, " + gene + " = " + str(counter) + "-" + str(
                counter + length) + '\n')  # New line of text with gene and length
            counter = counter + length + 1  # Keep record of the last position used
            outfile.close()
