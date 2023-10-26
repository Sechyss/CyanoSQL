# %%
import os

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# %% Reference genomes
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Temporal')
WH7803_genome = SeqIO.parse('Syn_WH7803.gbk', 'genbank')
Spmd2d_genome = SeqIO.parse('SPM2_genome.gb', 'genbank')

# %%

for genome in WH7803_genome:
    f1 = open("Reference_Proteins_WH7803.faa", "w")
    for feature in genome.features:
        if feature.type == "CDS":  # Apply to those that are not manually created and might have CyaClusterID
            try:
                ID = str(feature.qualifiers["locus_tag"][0])
                ntseq2 = feature.extract(genome.seq)
                aaseq = str(ntseq2.translate(table="Bacterial", cds=False, stop_symbol="", gap=None))
                FastaAA = SeqRecord(Seq(aaseq),
                                    id=ID + '_' + str(feature.qualifiers['product'][0]).replace(' ', '_'),
                                    description='[Synechococcus WH7803]')
                SeqIO.write(FastaAA, f1, "fasta")
            except:
                pass

for genome in Spmd2d_genome:
    f1 = open("Reference_Proteins_SPM2.faa", "w")
    for feature in genome.features:
        if feature.type == "CDS":  # Apply to those that are not manually created and might have CyaClusterID
            ID = str(feature.qualifiers["locus_tag"][0])
            try:
                product = str(feature.qualifiers['product'][0]).replace(' ', '_')
            except:
                product = 'No_Product'
            ntseq2 = feature.extract(genome.seq)
            aaseq = str(ntseq2.translate(table=11, cds=False, stop_symbol="", gap=None))
            FastaAA = SeqRecord(Seq(aaseq), id=ID + '_' + product, description='[Phage SPM2]')
            SeqIO.write(FastaAA, f1, "fasta")

# %%

Accession_Phage = pd.read_csv('BLASTP_Phage.tab', sep='\t', header=None)
Blastp_WH7803 = Accession_Phage[1].to_list()

for i in range(len(Blastp_WH7803)):
    Blastp_WH7803[i] = Blastp_WH7803[i][0:14]
for genome in WH7803_genome:
    f1 = open("ToBlastp_Proteins_WH7803.faa", "w")
    for feature in genome.features:
        if feature.type == "CDS":  # Apply to those that are not manually created and might have CyaClusterID
            try:
                ID = str(feature.qualifiers["locus_tag"][0])
                if ID in Blastp_WH7803:
                    ntseq = feature.extract(genome.seq)
                    aaseq = str(ntseq.translate(table="Bacterial", cds=False, stop_symbol="", gap=None))
                    FastaAA = SeqRecord(Seq(aaseq),
                                        id=ID + '_' + str(feature.qualifiers['product'][0]).replace(' ', '_'),
                                        description='[Synechococcus WH7803]')
                    SeqIO.write(FastaAA, f1, "fasta")
                else:
                    pass
            except:
                pass

# %% Analysis of examples of gene clusters vs SPM2

f1 = open('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Temporal/Cyanobacteriadb.faa',
          'w')
f1.close()

os.chdir(
    '/Users/u1866168/Documents/OneDrive - University of '
    'Warwick/Experiments/GET_HOMOLOGUES/GH_Prokka_run/OMCL_alg/Prokka_OMCL_clusters_fna')
with open('/Users/u1866168/Documents/OneDrive - University of Warwick/'
          'Genome_Database/FinalSelectionGenomes/Prokka_selection/Database_Prokka_clusters.fna', 'a') as DB:
    for file in os.listdir():
        if file.endswith('.fna'):
            for seq_record in SeqIO.parse(file, 'fasta'):
                cluster = str(file.replace('.fna', ''))
                ID = cluster  # str(seq_record.id[3:] + '_from_' +
                print(ID)
                ntseq = seq_record.seq
                AA = SeqRecord(Seq(str(ntseq.translate(table=11, cds=False, stop_symbol="", gap=None))),
                               id=ID, description='')
                SeqIO.write(AA, DB, "fasta")
                break
# makeblastdb -dbtype prot -in CyanoDB.faa -input_type fasta -hash_index -out Cyanobacteria_db
# blastp -query Reference_Proteins_SPM2.faa  -db Cyanobacteria_db -out BLASTP_Phage_Cyano -outfmt 6 -evalue 0.00001
# %% Preparation of genes of AMGs between myos and Cyanos
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Genome_Database/Temporal')

Accession_Cyano = pd.read_csv('BLASTP_myoclusters_Cyano.txt', sep='\t', header=None)
Blastp_Cyano = Accession_Cyano[1].to_list()

f1 = open('ToBlast_Cyanos_orth.faa', 'w')
f1.close()

with open('ToBlast_Cyanos_orth.faa', 'a') as DB:
    for seq_record in SeqIO.parse('Cyanobacteriadb.faa', 'fasta'):
        if seq_record.id.replace('.', '') in Blastp_Cyano:
            ID = str(seq_record.id)
            aaseq = seq_record.seq
            AA = SeqRecord(Seq(str(aaseq)), id=ID, description='')
            SeqIO.write(AA, DB, "fasta")

# scp -r Reference_Proteins_*  sls_alberto@137.205.68.160:~/Documents/Pipeline/CyaPhagePipeline/Host_Phage_Orth

# %% Selection of genes from cyanos to double blastp agains podos

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'Experiments/GET_HOMOLOGUES/CyanoPhages/Podoviruses/Prokka_curated/AMGs_Podos/')

TablesBlastP = pd.read_csv('CyanoPodosvsHost_blastp.tsv', sep='\t', header=None)
Blastp_Cyano = TablesBlastP[1].to_list()
with open('ToBlast_CyanosvsPodos_orth.faa', 'a') as DB:
    for seq_record in SeqIO.parse('Cyanobacteriadb.faa', 'fasta'):
        if seq_record.id.replace('.', '') in Blastp_Cyano:
            ID = str(seq_record.id)
            aaseq = seq_record.seq
            AA = SeqRecord(Seq(str(aaseq)), id=ID, description='')
            SeqIO.write(AA, DB, "fasta")

# %% Extraction of seqs for Modelling of AMGs (Creation of fasta file for AMGs)
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick'
         '/Genome_Database/Phage_Host_Orth/Myos_Cyanobacteria/')

amg = pd.read_excel('Draft_Orth_Myos_Cyanobacteria.xlsx', sheet_name='Draft_AMGs')

list_amgs = list(amg['Cyanobacteria_cluster'])

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
         'Experiments/GET_HOMOLOGUES/FinalDraftGetHomologues/OMCL_Clusters/')
for file in os.listdir('/Users/u1866168/Documents/OneDrive - University of Warwick/'
                       'Experiments/GET_HOMOLOGUES/FinalDraftGetHomologues/OMCL_Clusters/'):

    db = open('/Users/u1866168/Documents/OneDrive - University of Warwick/'
              '/Genome_Database/Phage_Host_Orth/Myos_Cyanobacteria/Representatives_AMGs.fasta', 'a')

    if file.endswith('.faa') and str(file.replace('.faa', '')) in list_amgs:
        print(str(file.replace('.faa', '')))
        for seq_record in SeqIO.parse(file, 'fasta'):
            cluster = str(file.replace('.faa', ''))
            id = str(cluster)
            aaseq = seq_record.seq
            aa = SeqRecord(Seq(str(aaseq)), id=id, description='')
            SeqIO.write(aa, db, "fasta")
            db.close()
            break

# %% Extraction of seqs for modelling in Phyre2 (podovirus version)

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick'
         '/Genome_Database/Phage_Host_Orth/AMGs_Podos/')

Table_podo_amg_draft = pd.read_csv('FirstDraft_AMGs_PodosCyano.tsv', sep='\t', header=None)
list_amgs = Table_podo_amg_draft[0].to_list()

with open('Draft_AMGs_toModel.fasta', 'a') as outfile:
    for seq_record in SeqIO.parse('Cyanobacteriadb.faa', 'fasta'):
        if seq_record.id in list_amgs:
            ID = str(seq_record.id)
            aaseq = seq_record.seq
            AA = SeqRecord(Seq(str(aaseq)), id=ID, description='')
            SeqIO.write(AA, outfile, "fasta")
