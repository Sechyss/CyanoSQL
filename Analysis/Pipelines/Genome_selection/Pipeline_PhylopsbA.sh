#!/usr/bin/env bash
# shellcheck disable=SC2164
cd ~/Documents/Pipeline/Phylogenetic_treeSecondaryGenes

echo "Running PhyloSecondaryGenes"

for file in ~/Documents/Pipeline/Phylogenetic_treepetB/*.gbk  #Trying to apply script to every single genbank within the database

do    #Starting the code with Python language

python ~/Documents/Scripts/External_Scripts/CDS_FASTA_Extraction.py $file

mv ~/Documents/Pipeline/Phylogenetic_treeSecondaryGenes/CDS.fa $file.fna
mv ~/Documents/Pipeline/Phylogenetic_treeSecondaryGenes/AASeq.fa $file.fa

done

rename s/.gbk.fna/_CDS.fna/ *.gbk.fna
rename s/.gbk.fa/_AASeq.fna/ *.gbk.fa

cat *_CDS.fna > Cya_CDS.fna #Combine all the CDS from all strains
cat *_AASeq.fna > Cya_AASeq.fna
sed -e "s/ /_/g" Cya_CDS.fna> Cya_CDS_sed.fna #remove spaces in headers
sed -e "s/ /_/g" Cya_AASeq.fna> Cya_AASeq_sed.fna

usearch_v10.0.240 -usearch_global Cya_CDS_sed.fna -db ~/Documents/Pipeline/ReferenceFiles/psbA_nt_ProSS120.fna -id 0.75 -strand plus -matched Cya_psbA_nt.fna -threads 16 -target_cov 0.7 #use of ProSS120 as reference sequence to remove all remaining CDS
usearch_v10.0.240 -usearch_global Cya_AASeq_sed.fna -db ~/Documents/Pipeline/ReferenceFiles/psbA_AA_ProSS120.fna  -id 0.75 -strand plus -matched Cya_rpoC1_nt.fna -threads 16 -target_cov 0.7 #use of ProSS120 as reference sequence to remove all remaining CDS









                    



