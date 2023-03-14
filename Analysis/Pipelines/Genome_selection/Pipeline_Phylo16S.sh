#!/usr/bin/env bash
cd ~/Documents/Pipeline/Phylogenetic_tree16S

echo "Running Phylo16S"

for file in ~/Documents/Pipeline/Phylogenetic_tree16S/*.gbk  #Trying to apply script to every single genbank within the database

do    #Starting the code with Python language
 
python ~/Documents/Scripts/External_Scripts/rRNA_FASTA_Extraction.py $file

mv ~/Documents/Pipeline/Phylogenetic_tree16S/rRNAs.fa $file.fna

done

rename s/.gbk.fna/_rRNAs.fna/ *.gbk.fna

cat *_rRNAs.fna > Cya_rRNA.fna
sed -e "s/ /_/g" Cya_rRNA.fna> Cya_rRNA_sed.fna 

usearch_v10.0.240 -usearch_global Cya_rRNA_sed.fna -db ~/Documents/Pipeline/ReferenceFiles/16SrRNA_ProSS120.fna -id 0.85 -strand plus -matched Cya16S_rRNA.fna -query_cov 0.80

cat Cya16S_rRNA.fna 16SrRNA_Synelongatus.fna > Cya16S_rRNA_rooted.fna

muscle -in Cya16S_rRNA_rooted.fna -out Cya16S_rRNA_MUSCLE.fna

iqtree -s Cya16S_rRNA_MUSCLE.fna -st DNA -o Synechococcus_elongatus_16S_ribosomal_RNA_complete_sequence -nt 16 -m HKY -asr

newicktotxt Cya16S_rRNA_MUSCLE.fna.treefile

less Cya16S_rRNA_MUSCLE.txt  


echo "Enjoy your Wonderfull tree"


##!/usr/bin/env python

#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import SeqFeature
#import sys

#bank = SeqIO.parse(open(sys.argv[1],"rU"), "genbank") #Genbank parsing.
#output_handle = open("rRNAs.fa","w") 
#gbank2 = SeqIO.read(open(sys.argv[1],"rU"),"genbank")

#rRNAs = [] #Creating a list with all the features associated to each rRNA.
#for genome in gbank :   #Collection of features in rRNA list for each rRNA.
#    for feature in genome.features:
#        if (feature.type == "rRNA"):
#            ID =str(feature.qualifiers ["locus_tag"][0])
#            strain =str(gbank2.name)
#            seq = feature.extract(genome.seq)
#            record = SeqRecord (seq, id=str(strain+" "+ID), description="")
#            rRNAs.append(record)
            
#SeqIO.write (rRNAs, output_handle, "fasta") #Write the programme to compile everything in one single file.
#output_handle.close()


