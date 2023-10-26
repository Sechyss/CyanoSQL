import os

import pandas as pd
from SQL_analysis import fasta_4_tree, create_partition_file

Table_clusters_Core = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                                  'Warwick/Genome_Database/Metadata/Core_phylo_Prokka.csv', index_col="id")
fasta_4_tree('/Users/u1866168/Documents/OneDrive - University of '
             'Warwick/Experiments/GET_HOMOLOGUES/GH_Prokka_run/OMCL_alg/TranslatorX/', Table_clusters_Core)

# After this pipeline it is necessary to run MUSCLE to align the AA sequence prior to use translatorx
# for file in ./*fa; do seqkit rmdup < $file > ${file%.*}.faa ; done
# for file in ./*fasta; do seqkit rmdup < $file > ${file%.*}.fna ; done
# $ for file in *; do res=$(ls -leaf | grep -c '>' $file); if (($res != 110)); then rm $file; fi; done
# $ for file in ./*.faa; do muscle -in $file -out ${file%.*}.muscle; done
# $ for file in ./*.fna; do perl translatorx.pl -i $file -a ${file%.*}.muscle -o ${file%.*}.transx -c 11; done

# %% Partition file creation
f1 = open('Partition_file.txt', 'a')  # This will create the partition file for the SH test and RDP4
filelist = os.listdir('/Users/u1866168/Documents/OneDrive - University of '
                      'Warwick/Experiments/GET_HOMOLOGUES/GH_Prokka_run/OMCL_alg/TranslatorX/')
create_partition_file(filelist, f1)
# seqkit concat *nt_ali.fasta > Concatenate_Core_Prokka.fasta

# %%

# for file in ./*fasta; do iqtree -s $file -o Synechococcus_sp_RCC307,Synechococcus_CB0205,Synechococcus_CB101
# -nt 16 -m GTR+I+G -bb 1000; done


# %% Analysis of accessory genome- fasta for tree creation

# It requires remove CB101 and CB0205 strains from the fasta files

Table_clusters_Shell = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                                   'Warwick/Genome_Database/Metadata/Shell_SoftCore_phylo_Prokka.csv', index_col="id")
fasta_4_tree('/Users/u1866168/Documents/OneDrive - University of '
             'Warwick/Experiments/GET_HOMOLOGUES/GH_Prokka_run/OMCL_alg/TranslatorX/', Table_clusters_Shell)
