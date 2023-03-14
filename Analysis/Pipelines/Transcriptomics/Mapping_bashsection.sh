#!/usr/bin/env bash

cd ~/Documents/Pipeline/Transcriptomic_Analysis/Branko_Host_analysis/Branko_RNA_temporalseq/FASTQ_files || exit

for i in $(seq 1 60); do cat ${i}_*_R1_* > input1.fastq.gz; cat ${i}_*_R2_* > input2.fastq.gz; bbsplit.sh build=1 maxindel=100k minratio=0.7 ambiguous=toss in=input1.fastq.gz in2=input2.fastq.gz out=../Bamfiles/StrictParameters/${i}_mapped.bam; done

cd ~/Documents/Pipeline/Transcriptomic_Analysis/Branko_Host_analysis/Branko_RNA_temporalseq/Bamfiles/StrictParameters || exit

for i in $(seq 1 60); do samtools sort -o ${i}_onlymapped_sorted.bam -O BAM -l 5 --threads 16 ${i}_mapped_onlymapped.bam; done

for i in $(seq 1 60); do samtools index -@ 16 ${i}_onlymapped_sorted.bam; done

for i in $(seq 1 60); do bedtools bamtobed -i ${i}_onlymapped_sorted.bam > ../../Bedfiles/StrictParameters/${i}_tsm.bed; done