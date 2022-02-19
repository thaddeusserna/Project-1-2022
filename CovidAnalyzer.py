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


codon_counts = codon_counter(human_sequence)
amino_acid_counts = amino_acid_counter(human_sequence)
print(codon_counts)
print(amino_acid_counts)

# test method

table_header = ["AminoAcid", "Codon", "Number", "Ratio"]

table_body = []

table_body.append('arg')
table_body.append('cgu')
table_body.append(4)
table_body.append(2/8)



print(table_header)
print(table_body)

df = pd.DataFrame(table_body)
print(df)