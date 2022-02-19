import os

import pandas as pd
from collections import Counter
from typing import List, Any

print(os.listdir())
# save paths as strings
main_path = os.path.realpath('project 1')

human_genome_path = os.path.realpath('CovidSequenceHuman')

bat_genome_path = os.path.realpath('CovidSequenceBat')

# print paths
print(main_path)
print(human_genome_path)
print(bat_genome_path)

# Dictionary for Codon sequences that produce the same amino acid are in the same line

codon_dict = {'uuu': 'phe', 'uuc': 'phe',  # phe
              'uua': 'leu', 'uug': 'leu', 'cuu': 'leu', 'cuc': 'leu', 'cua': 'leu', 'cug': 'leu',  # leu
              'ucu': 'ser', 'ucc': 'ser', 'uca': 'ser', 'ucg': 'ser', 'agu': 'ser', 'agc': 'ser',  # ser
              'uau': 'tyr', 'uac': 'tyr',  # tyr
              'ugu': 'cys', 'ugc': 'cys',  # cys
              'ugg': 'trp',  # trp
              'ccu': 'pro', 'ccc': 'pro', 'cca': 'pro', 'ccg': 'pro',  # pro
              'cau': 'his', 'cac': 'his',  # his
              'caa': 'gin', 'cag': 'gin',  # gin
              'cgu': 'arg', 'cgc': 'arg', 'cga': 'arg', 'cgg': 'arg', 'aga': 'arg', 'agg': 'arg',  # arg
              'auu': 'lle', 'auc': 'lle', 'aua': 'lle',  # lle
              'aug': 'met',  # met
              'acu': 'thr', 'acc': 'thr', 'aca': 'thr', 'acg': 'thr',  # thr
              'aau': 'asn', 'aac': 'asn',  # asn
              'aaa': 'lys', 'aag': 'lys',  # lys
              'guu': 'val', 'guc': 'val', 'gua': 'val', 'gug': 'val',  # val
              'gcu': 'ala', 'gcc': 'ala', 'gca': 'ala', 'gcg': 'ala',  # ala
              'gau': 'asp', 'gac': 'asp',  # asp
              'gaa': 'glu', 'gag': 'glu',  # glu
              'ggu': 'gly', 'ggc': 'gly', 'gga': 'gly', 'ggg': 'gly',  # gly
              'uaa': 'stop', 'uag': 'stop', 'uga': 'stop'  # stop
              }
os.chdir(human_genome_path)
print(os.listdir())

infile = open("./sequenceHuman1.fasta")

print(infile)

human_sequence = 'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTG'

print(human_sequence)


def codon_counter(sequence):
    """

    :param sequence:
    :return:
    """
    codon_list = []
    for i in range(len(sequence) // 3):
        index: int = i * 3
        curr_sequence = sequence[index: index + 3].lower()

        if curr_sequence in codon_dict:
            codon_list.append(curr_sequence)
    codon_count = Counter(codon_list)
    return codon_count


def amino_acid_counter(sequence):
    """

    :param sequence:
    :return:
    """
    amino_acid_list = []
    for i in range(len(sequence) // 3):
        index: int = i * 3
        curr_sequence = sequence[index: index + 3].lower()

        if curr_sequence in codon_dict:
            amino_acid_list.append(codon_dict.get(curr_sequence.lower()))
    amino_acid_count = Counter(amino_acid_list)
    return amino_acid_count


# MAKE DataFrame
codon_counts = codon_counter(human_sequence)
amino_acid_counts = amino_acid_counter(human_sequence)
print(codon_counts)
print(amino_acid_counts)

# test method

# build Table
all_rows = [[]]


for elem in amino_acid_counts.elements():
    curr_amino_acid = elem
    curr_amino_acid_count = amino_acid_counts[curr_amino_acid]

    table_body = []
    for elem_codon in codon_counts.elements():
        if codon_dict.get(elem_codon) == elem:

            curr_codon_count = codon_counts[elem_codon]
            curr_codon_ratio = round((curr_codon_count / curr_amino_acid_count),3)

            table_body.append(elem)
            table_body.append(elem_codon)
            table_body.append(curr_codon_count)
            table_body.append(curr_codon_ratio)

            all_rows.append(table_body)
        continue    




# initialize DataFrame
df = pd.DataFrame(all_rows[1:], columns=("AminoAcid", "Codon", "Number", "Ratio"))

print(df)
