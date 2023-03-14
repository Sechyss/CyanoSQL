# %% Import of packages and functions
import csv
import glob
import os

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# %% Functions to updload

def genome_calling(a, b):  # Definition of the function to extract
    list_Cyanophage = {}  # Creation of the empty dictionary that will storage all the details.
    for index, row in Accession.iterrows():
        if ("Synechococcus" in row['Classification']) or ("Prochlorococcus" in row['Classification']) \
                or ("Cyanophage" in row['Classification']):  # First filter
            list_Cyanophage.update({index: [row['Description'], row['Classification']]})

    # Selection of a specific group -> a in this case

    for i in list_Cyanophage.keys():
        if a in list_Cyanophage[i][1]:
            b.update({i: list_Cyanophage[i]})


def mash_distance_estimation(a, b):  # This script to move to the folder containing the distances obtained after
    for file in os.listdir('../Distances_MASH/'):  # make list of files in the directory
        if file.endswith('tab'):
            acc_Dist = pd.read_csv(file, sep='\t', header=None)  # open the csv files
            for index, row in acc_Dist.iterrows():
                if row[3] < 1:  # if distance is not 1
                    acc_name = row[0]
                    if acc_name not in b:
                        a.append(acc_name)

    for i in a:  # this will store the information in a dictionary
        b.update({i: [Accession.at[i, 'Description'], Accession.at[i, 'Classification'],
                      Accession.at[i, 'Genome Length(bp)']]})


def locustag_extraction(a):
    for file in os.listdir():
        if file.endswith('.faa'):
            cluster = str(file).replace('.faa', '')  # storage of the variable name and the code for the cluster
            cluster2 = cluster.replace("'",
                                       "")  # this step will remove ' from the cluster's name, later will crash if not.
            Array1 = [seq_record.id[3:] for seq_record in SeqIO.parse(file, "fasta")]  # Extraction of ID as
            # dictionary keys
            Array2 = [cluster2] * (len(Array1))  # Cluster associated to ID dictionary

            for i in range(len(Array1)):  # Construction of the dictionary and updating of the dictionary
                Addition = {Array1[i]: Array2[i]}
                a.update(Addition)


def fraction_estimation(csv, a, b, Dict):
    for index, row in csv.iterrows():  # This will determine the fraction upon the frequency
        if row['Total'] == a:
            Dict.update({row['Gene']: 'Core'})
        elif row['Total'] == b:
            Dict.update({row['Gene']: 'SoftCore'})
        elif row['Total'] <= 2:
            Dict.update({row['Gene']: 'Cloud'})
        else:
            Dict.update({row['Gene']: 'Shell'})


def full_extraction_clusters(csv, Dict, FractDict):
    csv.writerow(['Locustag', 'ProteinID', 'Product', 'Cluster', 'Fraction', 'NTSEQ', 'AASEQ', 'Organism'])  # Headers
    organism = []  # This creates a variable outside the loop so it can be save for the try exception.
    for file in os.listdir():
        if file.endswith('.gb'):
            gbk = SeqIO.parse(file, 'genbank')
            for genome in gbk:
                for feature in genome.features:
                    if feature.type == 'source':
                        organism = str(feature.qualifiers["organism"][0])  # Extract the organism from source
                    if feature.type == 'CDS':
                        try:  # Extraction of main features, some genbank do not have locustag and it requires a try
                            ID = str(feature.qualifiers["protein_id"][0])
                            ID2 = str(feature.qualifiers["locus_tag"][0])
                            product = str(feature.qualifiers["product"][0])
                            ntseq = str(feature.extract(genome.seq))
                            ntseq2 = feature.extract(genome.seq)
                            aaseq = str(ntseq2.translate(table="Bacterial", cds=False, stop_symbol="#", gap=None))
                            cluster = []
                            fraction = []
                            if ID in Dict.keys():
                                cluster = Dict[ID]
                                fraction = FractDict[Dict[ID]]
                            elif ID2 in Dict.keys():
                                cluster = Dict[ID2]
                                fraction = FractDict[Dict[ID2]]
                            GeneCsv.writerow([ID2, ID, product, cluster, fraction, ntseq, aaseq, organism])

                        except:
                            ID = str(feature.qualifiers["protein_id"][0])
                            cluster = Dict[ID]
                            fraction = FractDict[Dict[ID]]
                            product = str(feature.qualifiers["product"][0])
                            ntseq = str(feature.extract(genome.seq))
                            ntseq2 = feature.extract(genome.seq)
                            aaseq = str(ntseq2.translate(table="Bacterial", cds=False, stop_symbol="#", gap=None))
                            GeneCsv.writerow(['None', ID, product, cluster, fraction, ntseq, aaseq, organism])


def sequence_extraction_annotation(directory):
    os.chdir(directory)
    for file in os.listdir():
        if file.endswith('.faa'):
            db = open(file.replace('.faa', '.fasta'), 'w')
            for seq_record in SeqIO.parse(file, 'fasta'):
                cluster = str(file.replace('.faa', ''))
                id = str(cluster)
                ntseq = seq_record.seq
                aa = SeqRecord(Seq(str(ntseq)), id=id, description='')
                SeqIO.write(aa, db, "fasta")
                db.close()
                break


# %% Extraction of Preliminary Cyanophage from Accession text file using genome_calling function
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cyanophage_Andy/Original_data/')

Accession = pd.read_csv('29Feb2020_phages.complete_genomes_accession.txt', index_col='Accession', sep='\t')

List_genomes = SeqIO.parse('29Feb2020.fasta', 'fasta')  # Fasta files with all the genomes
# %% Selection of only Podoviruses
CyanoPodos = {}
genome_calling('Podoviridae', CyanoPodos)  # Creates a dictionary with all the sequences from putative cyanopodos
# %% Creation of fasta files of the previous list to query the mash database
for genome in List_genomes:
    if genome.name in CyanoPodos.keys():
        print(genome.name)
        Podocya = open("CyanoPodos/" + genome.name.rsplit(".", 1)[0] + ".fna", 'w')
        SeqIO.write(genome, Podocya, "fasta")
        Podocya.close()

# Need to run mash with all these sequences vs the database
# %% Analysis of mash files to get close relatives

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Genome_Database/Cyanophage_Andy/CyanoPodos/Distances_MASH')

List_Cyanpodos = []
CyanoPodos = {}
mash_distance_estimation(List_Cyanpodos, CyanoPodos)
Dataframe_CyaPodo = pd.DataFrame.from_dict(CyanoPodos, orient='index_genome', columns=['Description', 'Classification',
                                                                                       'Gnm Length'])

Dataframe_CyaPodo.to_csv('First_draft_CyanoPodos.csv')

# Once we get the final list of genomes, we download them, and run GET_HOMOLOGUES
# %% Extraction of Locustag from the gene clusters

# To run the following pieces of code it requires the files from GET_HOMOLOGUES
os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/OMCL/Cluster_OMCL')

Final_dict_CyanoPodos = {}
locustag_extraction(Final_dict_CyanoPodos)

# %% Analysis of the fraction per gene-gene cluster

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/OMCL/')

# The following matrix requires to remove the '.faa' suffix and sum all cases per cluster
Dataframe_fraction = pd.read_csv('pangenome_matrix_t0.tr.csv')

# %% Fraction extraction for CyanoPodos
Fraction_dict_CyanoPodos = {}
fraction_estimation(Dataframe_fraction, 19, 18, Fraction_dict_CyanoPodos)

# %% Extraction of protein sequences from Podoviruses Cyanophage and creation of the csv table with all the data
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cyanophage_Andy/CyanoPodos_gbk')
Tmpfile = open("/Users/u1866168//Documents/OneDrive - University of "
               "Warwick/MHC_Scripts/Cyanobacteria/Python_results/FullListGenesCyanPodos.csv", "w")  # Creation of the file
GeneCsv = csv.writer(Tmpfile)
full_extraction_clusters(GeneCsv, Final_dict_CyanoPodos, Fraction_dict_CyanoPodos)
Tmpfile.close()
# %% Pipeline for Siphoviruses
CyanoSiphos = {}
genome_calling("Siphoviridae", CyanoSiphos)
List_genomes = SeqIO.parse('29Feb2020.fasta', 'fasta')
# %% Creation of fna files for mash query
for genome in List_genomes:
    if genome.name in CyanoSiphos.keys():
        print(genome.name)
        Siphocya = open("CyanoSiphos/" + genome.name.rsplit(".", 1)[0] + ".fna", 'w')
        SeqIO.write(genome, Siphocya, "fasta")
        Siphocya.close()
# for file in ./CyanoSiphos/*fna; do mash.2 dist ./29Feb2020.fasta.msh ${file} > ${file}.tab; done
# %% Analysis of mash files to get close relatives

os.chdir(
    '/Users/u1866168/Documents/OneDrive - University of '
    'Warwick/Genome_Database/Cyanophage_Andy/CyanoSiphos/Distances_MASH')
List_Cyanosiphos = []
Siphos = {}

mash_distance_estimation(List_Cyanosiphos, Siphos)
Dataframe_CyaSiph = pd.DataFrame.from_dict(Siphos, orient='index_genome', columns=['Description', 'Classification',
                                                                                   'Gnm Length'])

Dataframe_CyaSiph.to_csv('First_draft_CyanoSiphos.csv')

# Once we get the final list of genomes, we download them from NCBI and run GET_HOMOLOGUES
# %% Extraction of Locustag from the gene clusters

# To run the following pieces of code it requires the files from GET_HOMOLOGUES
os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Siphoviruses/OMCL/Cluster_OMCL')
Final_dict_CyanoSiphos = {}
locustag_extraction(Final_dict_CyanoSiphos)
# %% Analysis of the fraction per gene-gene cluster

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Siphoviruses/OMCL/')

# The following matrix requires to remove the '.faa' suffix and sum all cases per cluster
Dataframe_fraction_Sipho = pd.read_csv('pangenome_matrix_t0.tr.csv')

Fraction_dict_CyanoSiphos = {}
fraction_estimation(Dataframe_fraction_Sipho, 5, 4, Fraction_dict_CyanoSiphos)

# %% Extraction of protein sequences from Siphoviruses Cyanophage and creation of the csv table with all the data
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cyanophage_Andy/CyanoSiphos_gbk')

Tmpfile = open("/Users/u1866168//Documents/OneDrive - University of "
               "Warwick/MHC_Scripts/Cyanobacteria/Python_results/FullListGenesCyanSiphos.csv", "w")  # Creation of the file
GeneCsv = csv.writer(Tmpfile)
full_extraction_clusters(GeneCsv, Final_dict_CyanoSiphos, Fraction_dict_CyanoSiphos)
Tmpfile.close()

# %% Analysis for Myoviruses

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cyanophage_Andy')
CyanoMyos = {}
genome_calling('Myoviridae', CyanoMyos)

# %% Creation of fna files for mash query
for genome in List_genomes:
    if genome.name in CyanoMyos.keys():
        print(genome.name)
        Myocya = open("CyanoMyos/" + genome.name.rsplit(".", 1)[0] + ".fna", 'w')
        SeqIO.write(genome, Myocya, "fasta")
        Myocya.close()
# %% MASH distances
os.chdir(
    '/Users/u1866168/Documents/OneDrive - University of '
    'Warwick/Genome_Database/Cyanophage_Andy/CyanoMyos/Distances_MASH')
List_CyanMyos = []
CyanMyos = {}
mash_distance_estimation(List_CyanMyos, CyanMyos)
Dataframe_CyaMyos = pd.DataFrame.from_dict(CyanMyos, orient='index_genome', columns=['Description', 'Classification',
                                                                                     'Gnm Length'])

Dataframe_CyaMyos.to_csv('First_draft_CyanoMyos.csv')

# %% Extraction of all sequences from genbank files to annotate gene clusters from Myoviruses to create a db for blastp

os.chdir(
    '/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Cyanophage_Andy/CyanoMyos/CyanoMyos_gbk')

Podos = SeqIO.parse('Concatenated_Group1_Myos_genomes.gbk', 'genbank')

file = open('Concatenated_Group1_Myos_genomes.fasta', 'a')
for genome in Podos:
    for feature in genome.features:
        if feature.type == 'CDS':
            try:
                ID = str(feature.qualifiers["locus_tag"][0])
            except:
                ID = 'No_specified'
            try:
                protein_id = str(feature.qualifiers["protein_id"][0])
            except:
                protein_id = 'No_specified'
            try:
                product = str(feature.qualifiers["product"][0])
            except:
                product = 'No_specified'
            aaseq = str(feature.qualifiers["translation"][0])
            sequence = Seq(aaseq)
            FASTA = SeqRecord(sequence, id=(ID + '_' + protein_id + '_' + product).replace(' ', '_'), description='')
            SeqIO.write(FASTA, file, 'fasta')
file.close()

# %% Extraction of all sequences from genbank files to annotate gene clusters from Podpviruses to create a db for blastp

os.chdir(
    '/Users/u1866168/Documents/OneDrive - University of '
    'Warwick/Genome_Database/Cyanophage_Andy/CyanoPodos/CyanoPodos_gbk/')

Podos = SeqIO.parse('Concatenated_CyanoPodos.gbf', 'genbank')

file = open('Concatenated_CyanoPodos.fasta', 'a')
for genome in Podos:
    for feature in genome.features:
        if feature.type == 'CDS':
            try:
                ID = str(feature.qualifiers["locus_tag"][0])
            except:
                ID = 'No_specified'
            try:
                protein_id = str(feature.qualifiers["protein_id"][0])
            except:
                protein_id = 'No_specified'
            try:
                product = str(feature.qualifiers["product"][0])
            except:
                product = 'No_specified'
            aaseq = str(feature.qualifiers["translation"][0])
            sequence = Seq(aaseq)
            FASTA = SeqRecord(sequence, id=(ID + '_' + protein_id + '_' + product).replace(' ', '_'), description='')
            SeqIO.write(FASTA, file, 'fasta')
file.close()

# %% Simple extraction of sequences for annotation of each fraction of Myoviruses

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Myoviruses/Prokka_curated'
                               '/core_gene_clusters/')

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Myoviruses/Prokka_curated'
                               '/softcore_gene_clusters/')

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Myoviruses/Prokka_curated'
                               '/shell_gene_clusters/')

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Myoviruses/Prokka_curated'
                               '/cloud_gene_clusters/')

# %% Simple extraction of sequences for annotation of each fraction of Podoviruses

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/Prokka_curated'
                               '/core_gene_clusters/')

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/Prokka_curated'
                               '/softcore_gene_clusters/')

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/Prokka_curated'
                               '/shell_gene_clusters/')

sequence_extraction_annotation('/Users/u1866168/Documents/OneDrive - University of '
                               'Warwick/Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/Prokka_curated'
                               '/cloud_gene_clusters/')

# %% bash code for blastp to annotate the clusters

# for file in ../annotation/*fasta; do
# blastp -query ${file}  -db Cyanobacteria_db -out {file} -outfmt 6 -evalue 0.00001; done

# %% Extraction of annotation from blastp and HMMscan tables (Myoviruses)

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick'
         '/Experiments/GET_HOMOLOGUES/CyanoPhages/Myoviruses/Prokka_curated/annotation/')

data = pd.read_csv('blastp_softcore_myos.txt', sep="\t", header=None)

data2 = pd.read_csv('cloud_genes_clusters_pfam.csv', header=None)

data3 = pd.read_csv('shell_genes_clusters_pfam.csv', header=None)

dict_annotation = {}
for index, row in data.iterrows():
    if row[0] not in dict_annotation.keys():
        dict_annotation.update({str(row[0]).replace(' ', ''): str(row[1]).rsplit(".", 1)[1].replace('1_', '')})

for index, row in data2.iterrows():
    if row[2] not in dict_annotation.keys():
        dict_annotation.update({str(row[2]): [row[3]]})
    else:
        dict_annotation[row[2]].append(row[3])

for index, row in data3.iterrows():
    if row[2] not in dict_annotation.keys():
        dict_annotation.update({str(row[2]): [row[3]]})
    else:
        dict_annotation[row[2]].append(row[3])

for file in glob.glob('../*_gene_clusters/*.faa'):
    if file.rsplit(".", 1)[0].rsplit('/', 1)[1] not in dict_annotation.keys():
        dict_annotation.update({file.rsplit(".", 1)[0].rsplit('/', 1)[1]: "Unannotated_gene"})

dataframe = pd.DataFrame.from_dict(dict_annotation, orient='index_genome')
dataframe.to_csv('Annotation_collection.csv')

# %% Extraction of annotation from blastp and HMMscan tables (Podoviruses)

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick'
         '/Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/Prokka_curated/annotation/')

data = pd.read_csv('blastp_softcore_podos.txt', sep="\t", header=None)

data2 = pd.read_csv('cloud_genes_clusters_pfam.csv', header=None)

data3 = pd.read_csv('shell_genes_clusters_pfam.csv', header=None)

dict_annotation = {}
for index, row in data.iterrows():
    if row[0] not in dict_annotation.keys():
        dict_annotation.update({str(row[0]).replace(' ', ''): str(row[1]).rsplit(".", 1)[1].replace('1_', '')})

for index, row in data2.iterrows():
    if row[2] not in dict_annotation.keys():
        dict_annotation.update({str(row[2]): [row[4]]})
    else:
        dict_annotation[row[2]].append(row[4])

for index, row in data3.iterrows():
    if row[2] not in dict_annotation.keys():
        dict_annotation.update({str(row[2]): [row[4]]})
    else:
        dict_annotation[row[2]].append(row[4])

for file in glob.glob('../*_gene_clusters/*.faa'):
    if file.rsplit(".", 1)[0].rsplit('/', 1)[1] not in dict_annotation.keys():
        dict_annotation.update({file.rsplit(".", 1)[0].rsplit('/', 1)[1]: "Unannotated_gene"})

dataframe = pd.DataFrame.from_dict(dict_annotation, orient='index_genome')
dataframe.to_csv('Annotation_collection.csv')

# %% Analysis of the promoter region for the EMSA genes that gave shift

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'Genome_Database/Cyanophage_Andy/CyanoMyos/Genomes_Andy_Myos/')

# List of genes that gave shift
list_genes_EMSA = ['4_unannotated_protein.faa', '12451_unannotated_protein.faa', '12452_unannotated_protein.faa']

# Creation of variables to store the main information
list_genomes_cluster = []
list_emsa_gene = []
dict_emsa_gene = {}
dictionary_to_df = {}

# Collection of information from the different gene_clusters
for i in list_genes_EMSA:
    emsa_cluster = SeqIO.parse('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/'
                               'GET_HOMOLOGUES/CyanoPhages/Myoviruses/OMCL_clusters/OMCL/' + i, 'fasta')
    print(i)
    for seq_record in emsa_cluster:
        id = seq_record.id[3:11]  # This will extract the genomme name
        emsa_gene = seq_record.id[3:]  # This extract the locustag
        list_genomes_cluster.append(id)
        list_emsa_gene.append([emsa_gene])  # This list will be the reference when checking the locustag
        dict_emsa_gene.update({emsa_gene: i})  # Dictionary will store relationship locustag-genecluster
list_genomes_cluster = set(list_genomes_cluster)  # Remove of duplicates of the genomes

for y in list_genomes_cluster:
    gbk = SeqIO.parse('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                      'Genome_Database/Cyanophage_Andy/CyanoMyos/'
                      'Genomes_Andy_Myos/Myos_fastagenomes/' + y + '.gbf', 'genbank')
    for genome in gbk:
        for feature in genome.features:
            if feature.type == 'CDS':
                locustag = feature.qualifiers['locus_tag']
                if locustag in list_emsa_gene:
                    locustag = str(locustag).replace('\'', '').replace('[', '').replace(']', '')  # Remove list features
                    gene_location = feature.location.start  # Extraction start location of the gene of interest
                    sequence_300bp = gene_location - 300  # Extraction of the location of 100bp prior to gene
                    sequence_location = str(genome.seq[sequence_300bp:gene_location])  # 300bp 'Promoter' sequence
                    print(sequence_300bp, gene_location)
                    print(genome.seq[sequence_300bp:gene_location])
                    dictionary_to_df.update({locustag: [y, dict_emsa_gene[locustag]
                        , sequence_300bp, gene_location, sequence_location]})
dataframe_emsa = pd.DataFrame.from_dict(dictionary_to_df, orient='index', columns=['Genome', 'Cluster',
                                                                                   'Putative_promoter', 'Start_gene',
                                                                                   'Sequence'])
# Creation of the files to collect potential stepts of the data
# dataframe_emsa.to_csv('Putative_promoter_region_EMSA_shifted.csv')

multifasta = open('Putative_promoter_PhoBOX_fastas/Multifasta_putative_emsashifts.fasta', 'a')
spm2d004 = open('Putative_promoter_PhoBOX_fastas/Multifasta_putative_spm2d004.fasta', 'a')
spm2d133 = open('Putative_promoter_PhoBOX_fastas/Multifasta_putative_spm2d133.fasta', 'a')
spm2d134 = open('Putative_promoter_PhoBOX_fastas/Multifasta_putative_spm2d134.fasta', 'a')
for index, row in dataframe_emsa.iterrows():
    promoter_region = Seq(row['Sequence'])
    id_label = str(str(index) + '_' + row['Cluster'])
    seq_record = SeqRecord(promoter_region, id_label, description='')
    SeqIO.write(seq_record, multifasta, "fasta")
    if '4_unannotated_protein.faa' in row['Cluster']:
        SeqIO.write(seq_record, spm2d004, "fasta")
    elif '12451_unannotated_protein.faa' in row['Cluster']:
        SeqIO.write(seq_record, spm2d133, "fasta")
    else:
        SeqIO.write(seq_record, spm2d134, "fasta")
multifasta.close()
spm2d004.close()
spm2d133.close()
spm2d134.close()
