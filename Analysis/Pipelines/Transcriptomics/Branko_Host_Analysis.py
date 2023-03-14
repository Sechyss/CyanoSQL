# %% Import modules to work
import os

import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from MyPackage.GenbankEditing import GenbankEditing

# %% Creation of reference file to do the mapping
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cya_clusterPipeline')
WH7803_genome = SeqIO.parse('/Users/u1866168/Documents/OneDrive - University of '
                            'Warwick/Genome_Database/Synechococcus_Cyanorak/Syn_WH7803.gbk', 'genbank')

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Branko_Host_analysis')
f1 = open("Reference_Genome_Host_spmd2d.fna", "w")
f1.close()

for genome in WH7803_genome:
    f1 = open("Reference_Genome_Host_spmd2d.fna", "a")
    ntseq = Seq(genome.seq)
    FastaNT = SeqRecord(ntseq, id="CYKWH7803.1", description='')
    SeqIO.write(FastaNT, f1, "fasta")

SPM2_genome = SeqIO.parse('SPM2_genome.gb', 'genbank')
for genome in SPM2_genome:
    ntseq = genome.seq
    FastaNT = SeqRecord(ntseq, id="LN828717.1", description='')
    SeqIO.write(FastaNT, f1, "fasta")

# %% Samtools and transformation into bedfiles
# Move files into the server to analyse using the server
# bbsplit.sh ref=Reference_Genome_Host_spmd2d.fna build=1 path=./Branko_RNA_temporalseq

# for i in $(seq 1 60); do cat ${i}_*_R1_* > input1.fastq.gz; cat ${i}_*_R2_* > input2.fastq.gz; bbsplit.sh build=1
# in=input1.fastq.gz in2=input2.fastq.gz ambiguous=toss out=${i}_mapped.bam refstats=${i}_refstats.txt rpkm=${
# i}_rpkm.csv; done

# s samtools sort -o ${i}_mapped_sorted.bam -O BAM -leaf 5 --threads 16 ${i}_mapped.bam; done
# for i in $(seq 1 60); do samtools index_genome -@ 16 ${i}_mapped_sorted.bam; done

# for i in $(seq 1 60); do echo ${i}_mapped_sorted.bam;
# bedtools bamtobed -i ${i}_mapped_sorted.bam > ../Bedfiles/${i}_tsm.bed; done
# %% Creation of the GFF file for BEDTools
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Branko_Host_analysis')

WH7803_genome = SeqIO.parse('23615-Syn_WH7803.gb', 'genbank')
Spm2d_genome = SeqIO.parse('SPM2_genome.gb', 'genbank')
outfile = open('WH7803_Spm2d.gff', 'a')
GFF.write(WH7803_genome, outfile)
GFF.write(Spm2d_genome, outfile)
outfile.close()

# %% BEDtools script in bash shell
# gff2bed < Syn_WH7803.gff > Syn_WH7803.bed
# scp -r Syn_WH7803.* sls_alberto@137.205.68.160:~/Documents/Pipeline/Transcriptomic_Analysis/Branko_Host_analysis/
# for i in $(seq 1 60); do sed -i 's/LN828717/LN828717.1/g' ${i}_tsm.bed; done
# for i in $(seq 1 60); do bedtools intersect -c -bed -s -a ../../Syn_WH7803.bed -b ${i}_tsm.bed > ${i}_intersect.bed;
# ls ${i}_intersect.bed; done

# %% Analysis of Intersect bedfiles

# Include here the directory with the intersect

Dataframe_rw_reads = {}  # This Dictionary storages the data of all files
Dataframe_tpm = {}
Header = []
Gen_length = {}
for file in sorted(os.listdir()):
    if file.endswith('_intersect.bed') and not file.startswith('.'):
        print(file)
        Header.append(str(file).replace('_intersect.bed', ''))  # Extraction of file name for the sample

        bedfile = pd.read_csv(file, sep='\t', header=None)  # Read of CSV bedfiles

        Total_Reads = bedfile.iloc[2635][9] + bedfile.iloc[1][9]
        bedfile2 = bedfile[8].str.split(';', expand=True)  # Selection of column with locustag info

        Temporal_list = []
        for index, row in bedfile2.iterrows():  # Extraction of the locustag per position
            for i in row:
                if i is None:
                    continue
                if "locus_tag" in i and 'gene' not in bedfile.iloc[index][2]:
                    i = i.replace('locus_tag=', '')

                    if i not in Dataframe_rw_reads.keys():
                        Dataframe_rw_reads.update({i: [bedfile.iloc[index][9]]})
                        Gen_length.update({i: int(bedfile.iloc[index][4] - bedfile.iloc[index][3])})
                        Dataframe_tpm.update({i: [bedfile.iloc[index][9] / (
                                (Gen_length[i] / 1000) * (Total_Reads / 1000000))]})
                    else:
                        Dataframe_rw_reads[i].append(bedfile.iloc[index][9])
                        Dataframe_tpm[i].append(bedfile.iloc[index][9] / (
                                (Gen_length[i] / 1000) * (Total_Reads / 1000000)))

                    break

                if "locus_tag" not in i and \
                        ('tRNA' in bedfile.iloc[index][2] and 'LN828717' in bedfile.iloc[index][0]):
                    i = 'SPM2d-tRNA_' + str(bedfile.iloc[index][3]) + '_' + str(bedfile.iloc[index][4])

                    if i not in Dataframe_rw_reads.keys():
                        Dataframe_rw_reads.update({i: [bedfile.iloc[index][9]]})
                        Gen_length.update({i: int(bedfile.iloc[index][4] - bedfile.iloc[index][3])})
                        Dataframe_tpm.update({i: [bedfile.iloc[index][9] / (
                                (Gen_length[i] / 1000) * (Total_Reads / 1000000))]})
                    else:
                        Dataframe_rw_reads[i].append(bedfile.iloc[index][9])
                        Dataframe_tpm[i].append(bedfile.iloc[index][9] / (
                                (Gen_length[i] / 1000) * (Total_Reads / 1000000)))

                    break

# %%
HostAnalysis = pd.DataFrame.from_dict(Dataframe_rw_reads, orient='index')  # Creation of dataframe from Dict
TPM_normalized = pd.DataFrame.from_dict(Dataframe_tpm, orient='index')  # Creation of dataframe from Dict
# %% Creation of the CSV file from the dataframe
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Branko_Host_analysis/')
HostAnalysis.to_csv('All_features_Raw_reads.csv')  # Creation of the CSV file. It will require manual modification of
# 1st row
TPM_normalized.to_csv('All_features_TPM_normalized.csv')
# %% Creation of dictionary with gene length
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Branko_Host_analysis')
Gen_length = {}
GenbankEditing('23615-Syn_WH7803.gb').gene_length_dict(Gen_length)
GenbankEditing('SPM2_genome.gb').gene_length_dict(Gen_length)

# %% Preparation of Temporal models Plus P Non infected
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Branko_Host_analysis')
Dataframe = {}
Total_BAM = pd.read_excel('Total_Reads_BAM.xlsx', sheet_name='Total', Header=True)
Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P+_NoInf', index_col='Locustag')

for index, row in Data1.iterrows():
    #  Estimation of rpkm per gene.
    Serie1 = [int(row['0h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+P0-1'].sum() / 1000000)),
              int(row['3h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+P3-1'].sum() / 1000000)),
              int(row['6h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+P6-1'].sum() / 1000000)),
              int(row['12h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+P12-1'].sum() / 1000000)),
              int(row['15h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+P15-1'].sum() / 1000000))]
    Serie2 = [int(row['0h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+P0-2'].sum() / 1000000)),
              int(row['3h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+P3-2'].sum() / 1000000)),
              int(row['6h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+P6-2'].sum() / 1000000)),
              int(row['9h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+P9-2'].sum() / 1000000)),
              int(row['12h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+P12-2'].sum() / 1000000)),
              int(row['15h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+P15-2'].sum() / 1000000))]
    Serie3 = [int(row['0h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+P0-3'].sum() / 1000000)),
              int(row['3h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+P3-3'].sum() / 1000000)),
              int(row['6h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+P6-3'].sum() / 1000000)),
              int(row['9h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+P9-3'].sum() / 1000000)),
              int(row['12h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+P12-3'].sum() / 1000000)),
              int(row['15h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+P15-3'].sum() / 1000000))]
    # Extraction of coefficient, intercept and r score
    Dataframe.update({index: [Serie1[0], Serie2[0], Serie3[0], Serie1[1], Serie2[1], Serie3[1],
                              Serie1[2], Serie2[2], Serie3[2], Serie2[3], Serie3[3], Serie1[3],
                              Serie2[4], Serie3[4], Serie1[4], Serie2[5], Serie3[5]]})

# %% Preparation of Temporal model Minus P Non infected

Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P-_NoInf', index_col='Locustag')

for index, row in Data1.iterrows():
    #  Estimation of rpkm per gene and creation of Y array.
    Serie1 = [int(row['0h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-P0-1'].sum() / 1000000)),
              int(row['3h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-P3-1'].sum() / 1000000)),
              int(row['6h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-P6-1'].sum() / 1000000)),
              int(row['9h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-P9-1'].sum() / 1000000)),
              int(row['12h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-P12-1'].sum() / 1000000)),
              int(row['15h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-P15-1'].sum() / 1000000))]
    Serie2 = [int(row['0h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-P0-2'].sum() / 1000000)),
              int(row['3h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-P3-2'].sum() / 1000000)),
              int(row['6h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-P6-2'].sum() / 1000000)),
              int(row['9h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-P9-2'].sum() / 1000000)),
              int(row['12h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-P12-2'].sum() / 1000000)),
              int(row['15h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-P15-2'].sum() / 1000000))]
    Serie3 = [int(row['0h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-P0-3'].sum() / 1000000)),
              int(row['3h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-P3-3'].sum() / 1000000)),
              int(row['6h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-P6-3'].sum() / 1000000)),
              int(row['9h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-P9-3'].sum() / 1000000)),
              int(row['12h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-P12-3'].sum() / 1000000)),
              int(row['15h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-P15-3'].sum() / 1000000))]

    Dataframe[index].append(Serie1[0])  # Time 0 addition
    Dataframe[index].append(Serie2[0])
    Dataframe[index].append(Serie3[0])
    Dataframe[index].append(Serie1[1])  # Time 3 addition
    Dataframe[index].append(Serie2[1])
    Dataframe[index].append(Serie3[1])
    Dataframe[index].append(Serie1[2])  # Time 6 addition
    Dataframe[index].append(Serie2[2])
    Dataframe[index].append(Serie3[2])
    Dataframe[index].append(Serie1[3])  # Time 9 addition
    Dataframe[index].append(Serie2[3])
    Dataframe[index].append(Serie3[3])
    Dataframe[index].append(Serie1[4])  # Time 12 addition
    Dataframe[index].append(Serie2[4])
    Dataframe[index].append(Serie3[4])
    Dataframe[index].append(Serie1[5])  # Time 15 addition
    Dataframe[index].append(Serie2[5])
    Dataframe[index].append(Serie3[5])

# %% Preparation of the temporal model for Plus P Infected

Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P+_Inf', index_col='Locustag')

for index, row in Data1.iterrows():
    #  Estimation of rpkm per gene and creation of Y array.
    Serie1 = [int(row['0h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+P0-1'].sum() / 1000000)),
              int(row['3h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ3-1'].sum() / 1000000)),
              int(row['6h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ6-1'].sum() / 1000000)),
              int(row['9h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ9-1'].sum() / 1000000))]
    Serie2 = [int(row['0h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+P0-2'].sum() / 1000000)),
              int(row['3h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ3-2'].sum() / 1000000)),
              int(row['6h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ6-2'].sum() / 1000000)),
              int(row['9h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ9-2'].sum() / 1000000))]
    Serie3 = [int(row['0h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+P0-3'].sum() / 1000000)),
              int(row['3h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ3-3'].sum() / 1000000)),
              int(row['6h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ6-3'].sum() / 1000000)),
              int(row['9h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['+PΦ9-3'].sum() / 1000000))]

    Dataframe[index].append(Serie1[0])  # Time 0
    Dataframe[index].append(Serie2[0])
    Dataframe[index].append(Serie3[0])
    Dataframe[index].append(Serie1[1])  # Time 3
    Dataframe[index].append(Serie2[1])
    Dataframe[index].append(Serie3[1])
    Dataframe[index].append(Serie1[2])  # Time 6
    Dataframe[index].append(Serie2[2])
    Dataframe[index].append(Serie3[2])
    Dataframe[index].append(Serie1[3])  # Time 9
    Dataframe[index].append(Serie2[3])
    Dataframe[index].append(Serie3[3])

# %% Preparation of the temporal model for Minus P Infected

Data1 = pd.read_excel('Host_Analysis_RawReads.xls', sheet_name='P-_Inf', index_col='Locustag')

for index, row in Data1.iterrows():
    Serie1 = [int(row['0h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-P0-1'].sum() / 1000000)),
              int(row['6h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ6-1'].sum() / 1000000)),
              int(row['9h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ9-1'].sum() / 1000000)),
              int(row['12h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ12-1'].sum() / 1000000)),
              int(row['15h_1']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ15-1'].sum() / 1000000))]
    Serie2 = [int(row['0h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-P0-2'].sum() / 1000000)),
              int(row['3h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ3-2'].sum() / 1000000)),
              int(row['6h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ6-2'].sum() / 1000000)),
              int(row['9h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ9-2'].sum() / 1000000)),
              int(row['12h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ12-2'].sum() / 1000000)),
              int(row['15h_2']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ15-2'].sum() / 1000000))]
    Serie3 = [int(row['0h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-P0-3'].sum() / 1000000)),
              int(row['3h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ3-3'].sum() / 1000000)),
              int(row['6h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ6-3'].sum() / 1000000)),
              int(row['9h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ9-3'].sum() / 1000000)),
              int(row['12h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ12-3'].sum() / 1000000)),
              int(row['15h_3']) / ((Gen_length[index] / 1000) * (Total_BAM['-PΦ15-3'].sum() / 1000000))]

    Dataframe[index].append(Serie1[0])  # Time 0
    Dataframe[index].append(Serie2[0])
    Dataframe[index].append(Serie3[0])
    Dataframe[index].append(Serie2[1])  # Time 3
    Dataframe[index].append(Serie3[1])
    Dataframe[index].append(Serie1[1])  # Time 6
    Dataframe[index].append(Serie2[2])
    Dataframe[index].append(Serie3[2])
    Dataframe[index].append(Serie1[2])  # Time 9
    Dataframe[index].append(Serie2[3])
    Dataframe[index].append(Serie3[3])
    Dataframe[index].append(Serie1[3])  # Time 12
    Dataframe[index].append(Serie2[4])
    Dataframe[index].append(Serie3[4])
    Dataframe[index].append(Serie1[4])  # Time 15
    Dataframe[index].append(Serie2[5])
    Dataframe[index].append(Serie3[5])

TempDataframe = pd.DataFrame.from_dict(Dataframe, orient='index_genome')
TempDataframe.to_csv('Temporal_rpkm_Host.csv')
