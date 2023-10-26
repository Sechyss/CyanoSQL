#!/usr/bin/env bash

while read p; do cp /blastdb/Bacteriophage_monthly/GenomesDB/"$p"/*gbf Genomes_Andy_Podos; done <List_CyanoPodos_MASH.txt
while read p; do cp /blastdb/Bacteriophage_monthly/GenomesDB/"$p"/*gbf Genomes_Andy_Myos; done <List_Cyanomyos_MASH_1.txt
while read p; do cp /blastdb/Bacteriophage_monthly/GenomesDB/"$p"/*gbf Genomes_Andy_Siphos; done <List_CyanoSipho_MASH_4.txt