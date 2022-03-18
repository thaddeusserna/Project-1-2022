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
    :return Counter of seqince [codon:count]
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
     :return Counter of sequence [amino acid:count]
     """
    amino_acid_list = []
    for i in range(len(sequence) // 3):
        index: int = i * 3
        curr_sequence = sequence[index: index + 3].lower()

        if curr_sequence in codon_dict:
            amino_acid_list.append(codon_dict.get(curr_sequence.lower()))
    amino_acid_count = Counter(amino_acid_list)
    return amino_acid_count


def fasta_to_string(path):
    """
       :param path: String path for folder
       :return:String
       """

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


def fasta_string_list_maker(base_path):
    """
          :param base_path: String path for folder
          :return:list of sequnces as strings
          """

    return_list = []

    for entry in os.listdir(base_path):
        return_list.append(fasta_to_string(base_path + "/" + entry))

    return return_list


def amino_acid_list_maker(codon_list):
    amino_acid_list = []
    for entry in codon_list:
        amino_acid_list.append(codon_dict.get(entry.lower()))
    return amino_acid_list


def df_maker(base_path):
    """
    :param base_path:
    :return:dataFrame
    """
    data = fasta_string_list_maker(base_path)
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

df = df_maker(bat_path)
df_hum = df_maker(human_path)

# MAKING GRAPH
sns.set_palette('Paired')
df['Virus'] = 'HKU8'
df_hum['Virus'] = 'Covid-19'

res: DataFrame | Series = pd.concat([df, df_hum])

# Present Graph
fig, axs = plt.subplots(ncols=2)

# make scatter plot
sns.scatterplot(data=res, x='RSCU', y='Codon', hue='AminoAcid', style='Virus', ax=axs[0])

# make bar plot
ax = sns.barplot(x='RSCU', y='Codon', data=res, hue='Virus', ax=axs[1])
# show graphs
plt.show()

# export data to csv
res.to_csv('data_results.csv')

# STAT ANALYSIS
# look for correlation between RSCU HUMAN AND BAT
from sklearn.metrics import r2_score

bat_rscu = df['RSCU']
hum_rscu = df_hum['RSCU']

r2 = r2_score(hum_rscu, bat_rscu)
print('r^2 value:', r2)
