# B DATA 200
# CovidAnalyzer 1.5
# 3-15-22

# imports
import os
from collections import Counter
from typing import List, Any
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from CAI import RSCU

# Dictionary for Codon sequences that produce the same amino acid are in the same line
from pandas import DataFrame, Series

codon_dict = {'uuu': 'phe', 'uuc': 'phe',  # phe
              'uua': 'leu', 'uug': 'leu', 'cuu': 'leu', 'cuc': 'leu', 'cua': 'leu', 'cug': 'leu',  # leu
              'ucu': 'ser', 'ucc': 'ser', 'uca': 'ser', 'ucg': 'ser', 'agu': 'ser', 'agc': 'ser',  # ser
              'uau': 'tyr', 'uac': 'tyr',  # tyr
              'ugu': 'cys', 'ugc': 'cys',  # cys
              'ugg': 'trp',  # trp
              'ccu': 'pro', 'ccc': 'pro', 'cca': 'pro', 'ccg': 'pro',  # pro
              'cau': 'his', 'cac': 'his',  # his
              'caa': 'gln', 'cag': 'gln',  # gln
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


# Credit for get_list geeksForGeeks
# https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
def get_list(dict):
    return dict.keys()


def codon_counter(sequence):
    """

    :param sequence:
    :return Counter of :
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
    amino_acid_list = []
    for i in range(len(sequence) // 3):
        index: int = i * 3
        curr_sequence = sequence[index: index + 3].lower()

        if curr_sequence in codon_dict:
            amino_acid_list.append(codon_dict.get(curr_sequence.lower()))
    amino_acid_count = Counter(amino_acid_list)
    return amino_acid_count


def fasta_to_string(path):
    infile = open(path)
    temp_seq = ""
    infile.readline()[1:]
    for line in infile:
        temp_seq += line

    # Make divisible by three
    if len(temp_seq) % 3 == 1:
        temp_seq = temp_seq[1:]
    elif len(temp_seq) % 3 == 2:
        temp_seq = temp_seq[2:]

    return temp_seq


def fasta_list_maker(base_path):
    return_list = []

    for entry in os.listdir(base_path):
        return_list.append(fasta_to_string(base_path + "/" + entry))

    return return_list


def amino_acid_list_maker(codon_list):
    amino_acid_list = []
    for entry in codon_list:
        amino_acid_list.append(codon_dict.get(entry.lower()))
    return amino_acid_list


def df_maker(base_path, name):
    data = fasta_list_maker(base_path)
    # Calculate RSCU
    temp_rscu = RSCU(data)

    # Make lists to form dataframe
    # Codon List
    codon_list = get_list(temp_rscu)
    # RSCU List
    RSCU_list = list(temp_rscu.values())
    # Amino Acid List
    amino_acid_list = amino_acid_list_maker(codon_list)
    # return_df = pd.DataFrame.from_dict(temp_rscu, orient='index', columns=[name])

    data = list(zip(amino_acid_list, codon_list, RSCU_list))

    return_df = pd.DataFrame(data, columns=['AminoAcid', 'Codon', 'RSCU'])
    return return_df


# make list from fasta Files

# Save path names
bat_path = "CovidSequenceBat"
human_path = "CovidSequenceHuman"

df = df_maker(bat_path, "Bat RSCU")
df_hum = df_maker(human_path, "Human RSCU")

# MAKING GRAPH
#sns.set_palette('Paired')
df['Virus'] = 'HKU8'
df_hum['Virus'] = 'Covid-19'

res: DataFrame | Series = pd.concat([df, df_hum])

#ax = sns.barplot(x='AminoAcid', y='RSCU', data=res, hue='AminoAcid')
#plt.xticks(rotation=90)
#
# Label bar
#for i in ax.containers:
    # ax.bar_label(i,)
 #   ax.bar_label(i, rotation=45)

#compare_ax = sns.barplot(x='Codon', y='RSCU', data=res, hue='Virus')


# Present Graph
#fig, axs = plt.subplots(ncols=2)
#sns.barplot(x='AminoAcid', y='RSCU', data=res, hue='AminoAcid', ax=axs[1])
#sns.barplot(x='Codon', y='RSCU', data=res, hue='Virus', ax=axs[0])


sns.scatterplot(data=res,x='RSCU',y='Codon',hue='AminoAcid', style = 'Virus')
plt.legend(bbox_to_anchor=(1.0, 0.9))
plt.show()
