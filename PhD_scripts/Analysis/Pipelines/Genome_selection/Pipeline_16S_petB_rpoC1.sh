#!/usr/bin/env bash
cd ~/Documents/Pipeline/Phylogenetic_treeCya_Final_draft || exit

echo "Running Phylogenetic_Pipeline"

for file in ~/Documents/Pipeline/Phylogenetic_treeCya_Final_draft/*.gbk  #Apply script to every single genbank within the database

do    #Starting the code with Python language

python ~/Documents/Scripts/External_Scripts/CDS_FASTA_Extraction.py $file  # Extraction of the CDS in nt and aa

mv ~/Documents/Pipeline/Phylogenetic_treeCya_Final_draft/CDS.fa $file.fna  # Change name of temporary file for CDS ntseq
mv ~/Documents/Pipeline/Phylogenetic_treeCya_Final_draft/AASeq.fa $file.fa  # Change name of temporary file for CDS aaseq

python ~/Documents/Scripts/External_Scripts/rRNA_FASTA_Extraction.py $file  # Extraction of rRNA from genbank

mv ~/Documents/Pipeline/Phylogenetic_treeCya_Final_draft/rRNAs.fa $file.fasta  # Change name of temporary file

done

# Rename files to remove extensions
rename s/.gbk.fna/_CDS.fna/ *.gbk.fna
rename s/.gbk.fa/_AASeq.fna/ *.gbk.fa
rename s/.gbk.fasta/_rRNAs.fna/ *.gbk.fasta

# Concatenation of every file and remove spaces from FASTA headers
cat *_CDS.fna > Cya_CDS.fna #Combine all the CDS from all strains
cat *_AASeq.fna > Cya_AASeq.fna
cat *_rRNAs.fna > Cya_rRNA.fna
sed -e "s/ /_/g" Cya_CDS.fna> Cya_CDS_sed.fna #remove spaces in headers
sed -e "s/ /_/g" Cya_AASeq.fna> Cya_AASeq_sed.fna
sed -e "s/ /_/g" Cya_rRNA.fna> Cya_rRNA_sed.fna

# Extraction of genes of interest using ProSS120 as reference sequence
usearch_v10.0.240 -usearch_global Cya_AASeq_sed.fna -db ~/Documents/Pipeline/ReferenceFiles/petB_AA_ProSS120.fna -id 0.7 -strand plus -matched Cya_petB_AA.fna -threads 16 -target_cov 0.7 #use of ProSS120 as reference sequence to remove all remaining CDS
usearch_v10.0.240 -usearch_global Cya_AASeq_sed.fna -db ~/Documents/Pipeline/ReferenceFiles/rpoC1_AA_ProSS120.fna  -id 0.7 -strand plus -matched Cya_rpoC1_AA.fna -threads 16 -target_cov 0.7 #use of ProSS120 as reference sequence to remove all remaining CDS
usearch_v10.0.240 -usearch_global Cya_rRNA_sed.fna -db ~/Documents/Pipeline/ReferenceFiles/16SrRNA_ProSS120.fna -id 0.8 -strand plus -matched Cya16S_rRNA.fna -threads 16 -target_cov 0.8

# Addition of S. elongatus as root
cat Cya_petB_AA.fna petB_AA_Synelongatus.fna > Cya_petB_AA_rooted.fna #Addition of freshwater cyanobacteria to use as root latter in tree
cat Cya_rpoC1_AA.fna rpoC1_AA_Synelongatus.fna > Cya_rpoC1_AA_rooted.fna #Addition of freshwater cyanobacteria to use as root latter in tree
cat Cya16S_rRNA.fna 16SrRNA_Synelongatus.fna > Cya16S_rRNA_rooted.fna

# All-against-all alignment. Preliminary step before setting the tree.
muscle -in Cya_petB_AA_rooted.fna -out Cya_petB_MUSCLE.fna
muscle -in Cya_rpoC1_AA_rooted.fna -out Cya_rpoC1_MUSCLE.fna
muscle -in Cya16S_rRNA_rooted.fna -out Cya16S_rRNA_MUSCLE.fna

# Set of parameters for the Phylogenetic tree.
iqtree -s Cya_petB_MUSCLE.fna -st AA -o Synechococcus_elongatus_petB_AASeq -nt AUTO -m TEST -asr
iqtree -s Cya_rpoC1_MUSCLE.fna -st AA -o Synechococcus_elongatus_rpoC1_AASeq -nt AUTO -m TEST -asr
iqtree -s Cya16S_rRNA_MUSCLE.fna -st DNA -o Synechococcus_elongatus_16S -nt AUTO -m HKY -asr



echo "Enjoy your Wonderfull tree"


