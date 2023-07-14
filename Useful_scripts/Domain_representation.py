#!/usr/bin/env python3

import itertools
# Install all packages before running main script, look Documentation file on how to do it.
import os
import sys

import pandas as pd
from Bio import SeqIO
from ete3 import Tree, SeqMotifFace, TreeStyle

from Cyanopackage.FastaEditing import edit_NCBI_ids
from Cyanopackage.FileMethods import file_from_list


def itol_annotation(sequence_file, outfile, hexcode, domainname, dataframe, domainform):
    inputlist = []
    for sequence in sequence_file:
        if sequence.id in set(sorted(dataframe['FastaID'])):
            print(sequence.id)
            idx = set(dataframe.index[dataframe['FastaID'] == sequence.id])
            for index in idx:
                start = str(dataframe.loc[index]['From'])
                end = str(dataframe.loc[index]['To'])

                inputlist.append(str(sequence.id) + ' , ' + str(
                    len(sequence.seq)) + ' , ' + domainform
                                 + '|' + start + '|' + end + '|' + hexcode + '|' + domainname + ' , ')

    file_from_list(inputlist, outfile)


# %%

# Import of the list with antismash id checked using literature.
os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/'
         'Psip1/Psip1_phylo/Aug22_analysis_trimmedMAGs/')
HMMresults_df = pd.read_excel('HMMScan_Psip1_Extract.xlsx', sheet_name='Extract')

seq_file = SeqIO.parse('BLASTP_Psip1_seqs.faa', 'fasta')  # Open the genbank file
listoutfile = 'List for iTol annotation.txt'
color = '#007d7d'
domain = 'Psip1 like'
form = 'TR'

itol_annotation(seq_file, listoutfile, color, domain, HMMresults_df, form)

# %% Creation of raw tree with all the data and sequence lengths.

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/'
         'Psip1/Psip1_phylo/Aug22_RefSeq/')
edit_NCBI_ids('RSEQ_Psip1e-5raw.txt', 'RefSeq_Psip1_edited_evalue-5.faa')
colorphylum = pd.read_excel('Pythontreecolor.xlsx', sheet_name='Sheet1')
Psip1HMMresults_df = pd.read_excel('Psip1_like_HMMscan.xlsx', sheet_name='Extract')
seq_file = SeqIO.parse('RefSeq_Psip1_edited_evalue-5_RichSeqs.faa', 'fasta')  # Open the genbank file
CDDHits = pd.read_excel('CDDhit_RJP_excel.xlsx', sheet_name='CDDhit_RJP')

forms = ["[]"]
colors = ['#C83D95', '#107C10', 'brown', '#5C2D91', '#A80000', '#FFCC00', '#8A2BE2', '#876A7C', '#EE3B3B', '#228B22'
          '#004753', '#4B4C4E', '#8B8B00', '#4A588A', '#002050', ' #FFAC00', '#F08080', '#EE9572', '#000080', '#800080',
          '#8B4726', '#F28A90', '#FFD700'
          ]
domains = list(set(Psip1HMMresults_df['Domain']) | set(CDDHits['Short name']))
combi = list(itertools.product(forms, colors))
motifs = dict(zip(domains, combi))

tree = Tree('RefSeq-5.newick.txt', format=2)
# Calculate the midpoint node
R = tree.get_midpoint_outgroup()
# and set it as tree outgroup
tree.set_outgroup(R)
ts = TreeStyle()
nstyle = NodeStyle()

listleafs = set(tree.get_leaf_names())
# Draws nodes as small red spheres of diameter equal to 10 pixels
nstyle["shape"] = "sphere"
nstyle["size"] = 5
nstyle["fgcolor"] = "darkred"

alldomaindict = {}
for index, row in Psip1HMMresults_df.iterrows():
    if Psip1HMMresults_df.loc[index]['FastaID'] not in list(alldomaindict.keys()):
        alldomaindict.update(
            {Psip1HMMresults_df.loc[index]['FastaID']: {'MarineCyano_PsiP1': [[Psip1HMMresults_df.loc[index]['From'],
                                                                               Psip1HMMresults_df.loc[index]['To']]]}})
    else:
        alldomaindict[Psip1HMMresults_df.loc[index]['FastaID']]['MarineCyano_PsiP1'].append(
            [Psip1HMMresults_df.loc[index]['From'],
             Psip1HMMresults_df.loc[index]['To']])

for index, row in CDDHits.iterrows():
    if CDDHits.loc[index]['Query'] not in list(alldomaindict.keys()):
        alldomaindict.update(
            {CDDHits.loc[index]['Query']: {CDDHits.loc[index]['Short name']: [CDDHits.loc[index]['From'],
                                                                              CDDHits.loc[index]['To']]}})
    else:
        alldomaindict[CDDHits.loc[index]['Query']].update(
            {CDDHits.loc[index]['Short name']: [CDDHits.loc[index]['From'],
                                                CDDHits.loc[index]['To']]})

for sequence in seq_file:
    if sequence.id in set(listleafs):
        print(sequence.id)
        proteinseq: str = str(sequence.seq)
        mixedmotifs = []
        for putativedomains in alldomaindict[sequence.id].keys():
            for element in alldomaindict[sequence.id][putativedomains]:
                if putativedomains == 'MarineCyano_PsiP1':
                    for hit in alldomaindict[sequence.id][putativedomains]:
                        From = hit[0]
                        To = hit[1]
                        appendtomixed_motifs = [
                            # seq.start, seq.end, shape, width, height, fgcolor, bgcolor, name
                            From, To, motifs[putativedomains][0], None, 10, motifs[putativedomains][1],
                            motifs[putativedomains][1], 'times|2|white|' + putativedomains,
                        ]
                        mixedmotifs.append(appendtomixed_motifs)

                else:
                    From = alldomaindict[sequence.id][putativedomains][0]
                    To = alldomaindict[sequence.id][putativedomains][1]
                    appendtomixed_motifs = [
                        # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
                        From, To, motifs[putativedomains][0], None, 10, motifs[putativedomains][1],
                        motifs[putativedomains][1], 'times|2|white|' + putativedomains,
                    ]
                    mixedmotifs.append(appendtomixed_motifs)
        mixedmotifs.sort()
        mixedmotifs = list(k for k, _ in itertools.groupby(mixedmotifs))
        # noinspection PyTypeChecker
        seqFace = SeqMotifFace(proteinseq, motifs=mixedmotifs, seq_format="-")
        (tree & sequence.id).add_face(seqFace, 0, "aligned")

dict_nodecolor = {}
phylums = list(set(colorphylum['Phylum']))

for group in phylums:
    rslt_df = colorphylum[colorphylum['Phylum'] == group]
    key = list(rslt_df['Color'].unique())[0]
    taxa = list(rslt_df['Taxa'].unique())
    dict_nodecolor.update({key: taxa})

nstyle["hz_line_width"] = 2
nstyle["vt_line_width"] = 2

for n in tree.traverse():
    n.set_style(nstyle)
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True

tree.show()

# %% Creation of the tree with relative distances (Seq = x*half length)


os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Analysis_Genes/'
         'Psip1/Psip1_phylo/Aug22_RefSeq/')
edit_NCBI_ids('RSEQ_Psip1e-5raw.txt', 'RefSeq_Psip1_edited_evalue-5.faa')
colorphylum = pd.read_excel('Pythontreecolor.xlsx', sheet_name='Sheet1')
Psip1HMMresults_df = pd.read_excel('Psip1_like_HMMscan.xlsx', sheet_name='Extract')
seq_file = SeqIO.parse('RefSeq_Psip1_edited_evalue-5_RichSeqs.faa', 'fasta')  # Open the genbank file
CDDHits = pd.read_excel('CDDhit_RJP_excel.xlsx', sheet_name='CDDhit_RJP')

forms = ["[]"]
colors = ['#D33F49', '#D7C0D0', '#EFF0D1', '#5C2D91', '#77BA99', '#262730', '#FB8B24', '#D90368', '#820263', '#291720'
          '#04A777', '#424B54', '#E1CE7A', '#EBCFB2', '#C5BAAF', '#17BEBB', '#2E282A', '#CD5334', '#EDB88B', '#FAD8D6',
          '#ECA400', '#006992', '#00A6A6'
          ]
domains = list(set(Psip1HMMresults_df['Domain']) | set(CDDHits['Short name']))
combi = list(itertools.product(forms, colors))
motifs = dict(zip(domains, combi))

tree = Tree('RefSeq-5.newick.txt', format=2)
# Calculate the midpoint node
R = tree.get_midpoint_outgroup()
# and set it as tree outgroup
tree.set_outgroup(R)
ts = TreeStyle()
nstyle = NodeStyle()

listleafs = set(tree.get_leaf_names())
# Draws nodes as small red spheres of diameter equal to 10 pixels
nstyle["shape"] = "sphere"
nstyle["size"] = 5
nstyle["fgcolor"] = "darkred"

alldomaindict = {}
for index, row in Psip1HMMresults_df.iterrows():
    if Psip1HMMresults_df.loc[index]['FastaID'] not in list(alldomaindict.keys()):
        alldomaindict.update(
            {Psip1HMMresults_df.loc[index]['FastaID']: {'MarineCyano_PsiP1': [[Psip1HMMresults_df.loc[index]['From'],
                                                                               Psip1HMMresults_df.loc[index]['To']]]}})
    else:
        alldomaindict[Psip1HMMresults_df.loc[index]['FastaID']]['MarineCyano_PsiP1'].append(
            [Psip1HMMresults_df.loc[index]['From'],
             Psip1HMMresults_df.loc[index]['To']])

for index, row in CDDHits.iterrows():
    if CDDHits.loc[index]['Query'] not in list(alldomaindict.keys()):
        alldomaindict.update(
            {CDDHits.loc[index]['Query']: {CDDHits.loc[index]['Short name']: [CDDHits.loc[index]['From'],
                                                                              CDDHits.loc[index]['To']]}})
    else:
        alldomaindict[CDDHits.loc[index]['Query']].update(
            {CDDHits.loc[index]['Short name']: [CDDHits.loc[index]['From'],
                                                CDDHits.loc[index]['To']]})

for sequence in seq_file:
    if sequence.id in set(listleafs):
        print(sequence.id)
        lenprotseq = len(str(sequence.seq))
        proteinseq: str = str('x'*int(lenprotseq/5))
        mixedmotifs = []
        for putativedomains in alldomaindict[sequence.id].keys():
            for element in alldomaindict[sequence.id][putativedomains]:
                if putativedomains == 'MarineCyano_PsiP1':
                    for hit in alldomaindict[sequence.id][putativedomains]:
                        From = hit[0]
                        RelativeFrom = int(From / 5)
                        To = hit[1]
                        RelativeTo = int(To / 5)
                        appendtomixed_motifs = [
                            # seq.start, seq.end, shape, width, height, fgcolor, bgcolor, name
                            RelativeFrom, RelativeTo, motifs[putativedomains][0], None, 10, motifs[putativedomains][1],
                            motifs[putativedomains][1], None,
                        ]
                        mixedmotifs.append(appendtomixed_motifs)

                else:
                    From = alldomaindict[sequence.id][putativedomains][0]
                    RelativeFrom = int(From / 5)
                    To = alldomaindict[sequence.id][putativedomains][1]
                    RelativeTo = int(To / 5)
                    appendtomixed_motifs = [
                        # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
                        RelativeFrom, RelativeTo, motifs[putativedomains][0], None, 10, motifs[putativedomains][1],
                        motifs[putativedomains][1], None,
                    ]
                    mixedmotifs.append(appendtomixed_motifs)
        mixedmotifs.sort()
        mixedmotifs = list(k for k, _ in itertools.groupby(mixedmotifs))
        # noinspection PyTypeChecker
        seqFace = SeqMotifFace(proteinseq, motifs=mixedmotifs, seq_format="-")
        (tree & sequence.id).add_face(seqFace, 0, "aligned")

dict_nodecolor = {}
phylums = list(set(colorphylum['Phylum']))

for group in phylums:
    rslt_df = colorphylum[colorphylum['Phylum'] == group]
    key = list(rslt_df['Color'].unique())[0]
    taxa = list(rslt_df['Taxa'].unique())
    dict_nodecolor.update({key: taxa})

tree.show(tree_style=ts)
