# Pipeline followed to map genes from phages and then analyse the data to normalize the TPM

# Creation of the index_genome bbsplit:
# $ bbsplit.sh build=1 ref=all_genes_095cdhit

# Script for BBsplit:
# $ for file in ./*gz; do bbsplit.sh build=1 in=$file out=${file%.*}.bam ambiguous2=split rpkm=${file%.*}.csv; done

# !/usr/bin/env python

# Loading of the required packages
import os

import pandas as pd

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/Transcriptomic_Analysis/Andy_phage_tpm')
Transcript_dict = {}  # Empty structure for the final dataframe
List_columns = []  # Empty structure for the columns in the final dataframe
#  The CSV files require the transformation of the txt into columns
#  For loop to iterate through files in the folder
for file in sorted(os.listdir('/Users/u1866168/Documents/OneDrive - University of '
                              'Warwick/MHC_Scripts/Transcriptomic_Analysis/Andy_phage_tpm')):
    if file.endswith('.csv'):  # Loading of the files into dataframe to load the data
        Dataframe = pd.read_csv(file, low_memory=False, header=None)
        #  Setting of the parameters to calculate the TPM
        Total_Normalized_RPK = 0
        Temp_dict = {}
        print("Counting the total of normalized reads for " + file.rsplit(".", 1)[0])
        for index, row in Dataframe.iterrows():  # Estimation of the total number of (read/gene-length)
            if index > 4:
                TPM = float(row[4]) / float(row[1])  # Division (reads/genelength)
                Total_Normalized_RPK = Total_Normalized_RPK + TPM  # Sum of the total number per station
        print("     Normalizing the number of reads for " + file.rsplit(".", 1)[0])
        for index, row in Dataframe.iterrows():
            if index > 4:
                key = row[0]  # Setting of the gene
                rpkm = row[5]  # Value of RPKM
                # Estimation of TPM and addition of the values to a temporal dictionary
                Temp_dict.update({key: [(float(row[4]) / float(row[1])) / (Total_Normalized_RPK / 1000000), rpkm]})
        # Creation of the columns of the CSV files
        List_columns.append(file.rsplit(".", 1)[0] + ".TPM")
        List_columns.append(file.rsplit(".", 1)[0] + ".RPKM")
        print("     Adding data into the dictionary")
        # Here we add individual data of every station into the final dataframe
        for i in Temp_dict.keys():
            if i in Transcript_dict.keys():
                Transcript_dict[i].append(Temp_dict[i][0])
                Transcript_dict[i].append(Temp_dict[i][1])
            else:
                Transcript_dict.update({i: Temp_dict[i]})
print("Editing the final dataframe")
# Creation of the final dataframe from the dictionary
Dict = pd.DataFrame.from_dict(Transcript_dict, columns=List_columns, orient='index_genome')  #
#  Creation of the CSV file
Dict.to_csv("Andy_phage_AMT.csv")
