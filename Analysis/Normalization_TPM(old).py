"""This is an old version but this may work"""
# Don't run the full script at once, it requires to do some manual steps in between processes
import os
import pickle

import pandas as pd

# Creation of reads in TPM and RPKM for locustag. NOT CLUSTERS
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/Transcriptomic_Analysis/RPKM')
Transcript_dict = {}  # Empty structure for the final dataframe
List_columns = []  # Empty structure for the columns in the final dataframe
#  The CSV files require the transformation of the txt into columns
#  For loop to iterate through files in the folder
for file in sorted(os.listdir('/Users/u1866168/Documents/OneDrive - University of '
                              'Warwick/MHC_Scripts/Transcriptomic_Analysis/RPKM')):
    if file.endswith('.csv'):  # Loading of the files into dataframe to load the data
        Dataframe = pd.read_csv(file, low_memory=False, header=None)
        Total_Reads = int(Dataframe[1][0])
        #  Setting of the parameters to calculate the TPM
        Total_Normalized_RPK = 0
        Temp_dict = {}
        print("Counting the total of normalized reads for " + file.rsplit(".", 1)[0])
        for index, row in Dataframe.iterrows():  # Estimation of the total number of (read/gene-length)
            if index > 3:
                TPM = float(row[4]) / float(row[1])  # Division (reads/genelength)
                Total_Normalized_RPK = Total_Normalized_RPK + TPM  # Sum of the total number per station
        print("     Normalizing the number of reads for " + file.rsplit(".", 1)[0])
        for index, row in Dataframe.iterrows():
            if index > 3:
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
Dict = pd.DataFrame.from_dict(Transcript_dict, columns=List_columns, orient='index_genome')
Dict.to_csv("SQL_transcripts.csv")
#  The CSV file require after this a little bit of work, needs to set the columns and fields


# Assignment of locustag into clusters from GET_HOMOLOGUES
Dictionary_clusters = open("/Users/u1866168/Documents/OneDrive - University of "
                           "Warwick/MHC_Scripts/Cyanobacteria/Python_results/Dict_ClusterAnalysis.pickle", "rb")
Translator = pickle.load(Dictionary_clusters)  # Load dictionary with Locustag:Clusters

# Open the CSV previously created
Trancript_database = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                                 'Warwick/MHC_Scripts/Transcriptomic_Analysis/RPKM/SQL_transcripts.csv')
Temporal_dict = {}  # Empty dictionary that will sum total TPM into clusters
for index, row in Trancript_database.iterrows():  # This row will sum all locustag of a particular cluster/all stations
    Key = Translator[row['Locus_Tag']][0]  # Translation into cluster
    Value = row['Total_RPM']  # Value of all stations
    if Key in Temporal_dict.keys():
        Sum = float(Temporal_dict[Key]) + float(Value)
        Temporal_dict.update({Key: Sum})
    else:
        Temporal_dict.update({Key: Value})
Dict = pd.DataFrame.from_dict(Temporal_dict, orient='index')
Dict.to_csv("Clusters_Transcriptomics.csv")

temporal_dict = {}  # Empty dictionary to storage TPM per station into clusters
for index, row in Trancript_database.iterrows():  # Creation of dictionary and further dataframe for the CSV file
    key = Translator[row['Locus_Tag']][0]
    if key not in temporal_dict:
        addition = {key: [row['AMT22_18_18_CT.TPM'], row['AMT22_3_4_CT.TPM'], row['AMT22_53_53_CT.TPM'],
                          row['AMT22_75_75_CT.TPM'], row['AMT22_75_75_SS.TPM'], row['AMT22_CTD4.TPM'],
                          row['AMT22_CTD75.TPM'], row['AMT23_11_15_CT.TPM'], row['AMT23_11_15_DS.TPM'],
                          row['AMT23_37_46_CT.TPM'], row['AMT23_37_46_DS.TPM'], row['AMT23_39_49_DS.TPM'],
                          row['AMT23_3_5_CT.TPM'], row['AMT23_3_5_DS.TPM'], row['AMT23_54_65_CT.TPM'],
                          row['AMT23_54_65_DS.TPM']]}
        temporal_dict.update(addition)
    else:
        addition = {
            key: [temporal_dict[key][0] + row['AMT22_18_18_CT.TPM'], temporal_dict[key][1] + row['AMT22_3_4_CT.TPM'],
                  temporal_dict[key][2] + row['AMT22_53_53_CT.TPM'], temporal_dict[key][3] + row['AMT22_75_75_CT.TPM'],
                  temporal_dict[key][4] + row['AMT22_75_75_SS.TPM'], temporal_dict[key][5] + row['AMT22_CTD4.TPM'],
                  temporal_dict[key][6] + row['AMT22_CTD75.TPM'], temporal_dict[key][7] + row['AMT23_11_15_CT.TPM'],
                  temporal_dict[key][8] + row['AMT23_11_15_DS.TPM'], temporal_dict[key][9] + row['AMT23_37_46_CT.TPM'],
                  temporal_dict[key][10] + row['AMT23_37_46_DS.TPM'],
                  temporal_dict[key][11] + row['AMT23_39_49_DS.TPM'],
                  temporal_dict[key][12] + row['AMT23_3_5_CT.TPM'], temporal_dict[key][13] + row['AMT23_3_5_DS.TPM'],
                  temporal_dict[key][14] + row['AMT23_54_65_CT.TPM'],
                  temporal_dict[key][15] + row['AMT23_54_65_DS.TPM']]}
        temporal_dict.update(addition)
Dict_clust = pd.DataFrame.from_dict(temporal_dict, orient='index_genome')
Dict_clust.to_csv("Cluster_TPM_stations.csv")
