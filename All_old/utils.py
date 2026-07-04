# Peptide Coverage Utilities ()

# Basic data processing Libraries
import pandas as pd
import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# For Data Visualization
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns


# For Data Visualization
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns

# For processing FASTA files
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Other Misc. libraries
import glob
import os
import math
import random
import re #regular expressions library
import warnings # To Turn off warnings in certain boxes or overall
import copy
from scipy.stats import ttest_ind
from aaindex import aaindex1


def fasta_import(app):
    records = []
    for fasta_file in app.fasta_files:
        records_i = list(SeqIO.parse(fasta_file, "fasta"))
        records = records + records_i

    # Convert to DataFrame (similar to as.data.frame(readAAStringSet()))
    fasta_df = pd.DataFrame({"names": [records.id for records in records],
        "prot_seq": [str(records.seq) for records in records]})

    # Split the "names" column to "review_type", "accession_number", and "protein_name"
    fasta_df[['review', 'accession', 'protein']] =  fasta_df['names'].str.split('|', expand = True)
    fasta_df = fasta_df.drop(['names'], axis = 1)
    fasta_df = fasta_df.drop_duplicates()

    # Add protein length ("prot_len") to fasta_df
    fasta_df['prot_len'] = fasta_df['prot_seq'].str.len()
    fasta_df

    return fasta_df


# Define a function to split and keep the first part (or add "area" back if needed)
def truncate_after_Area(col_name):
    # Split by 'area' (case-insensitive) and take the part before it
    parts = col_name.split('Area')
    if len(parts) > 1:
        # Keep the part before "Area"
        return parts[0] + 'Area'
    else:
        return col_name


def peptide_import_function(app):
    all_peptides_input_df = pd.DataFrame()
    if (app.input_type == "Manually Created"):
        print('"Manually Created" Import Method being used for Peptides.')
        rename_dict = {'Sequence': 'seq', '# PSMs': 'PSM', 'PSMs': 'PSM', 'Protein': 'protein',
               'Sequence': 'seq', 'Seq': 'seq', '# Peptides': 'peptides', 'Peptides': 'peptides',
               'Master Protein Accessions': 'accession', 'Modifications': 'mods', 'Mods': 'mods', 
               'Area': 'intensity', 'area': 'intensity', 'Intensity': 'intensity'}
        counter_i = 0  # zero a counter
        number_of_input_files = len(app.peptide_file_names)
        print(f"number of input files is: {number_of_input_files}")
        for sample_i, path_name_i in zip(app.peptide_file_names, app.peptide_files):
            print(path_name_i['path'])
            df_i = pd.read_csv(path_name_i['path'])
            #print(df_i.columns)
            df_i.rename(columns = rename_dict, inplace=True)
            df_i = df_i.reset_index()
            print(df_i.columns)
            df_i = df_i[['seq', 'accession', 'intensity', 'PSM']]
            df_i = df_i.groupby('seq').agg(PSM = ('PSM', 'sum'), accession = ('accession', 'first'),
                                   intensity =('intensity', 'sum')).reset_index()
            df_i['sample'] = sample_i
            all_peptides_input_df = pd.concat([all_peptides_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        print("Done importing files")
        return all_peptides_input_df
    elif (app.input_type == "MaxQuant"):
        ####### MaxQuant Results have no "mods Column
        print('"MaxQuant" Import method being used for peptides.')
        rename_dict = {'MS/MS Count': 'PSM', 'Proteins': 'accession',
               'Sequence': 'seq', 'Seq': 'seq', '# Peptides': 'peptides', 'Peptides': 'peptides',
               'Master Protein Accessions': 'accession', 'Modifications': 'mods', 'Mods': 'mods', 
               'Area': 'intensity', 'area': 'intensity', 'Intensity': 'intensity'}
        counter_i = 0  # zero a counter
        number_of_input_files = len(app.peptide_file_names)
        print(f"number of input files is: {number_of_input_files}")
        for sample_i, peptide_file_i in zip(app.peptide_file_names, app.peptide_files):
            path_name_i = peptide_file_i['path']
            df_i = pd.read_csv(path_name_i, sep = '\t', low_memory=False)
            #df_tab = pd.read_csv('data.txt', sep='\t')
            df_i.rename(columns = rename_dict, inplace = True)
            #df_i = df_i[['seq', 'mods', 'accession', 'intensity', 'PSM']]
            df_i = df_i[['seq', 'accession', 'intensity', 'PSM']]
            df_i = df_i.groupby('seq').agg(PSM = ('PSM', 'sum'), accession = ('accession', 'first'),
                                   intensity =('intensity', 'sum')).reset_index()
            df_i['sample'] = sample_i
            all_peptides_input_df = pd.concat([all_peptides_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        print("Done importing files")
        return all_peptides_input_df
    elif (app.input_type == "Proteome Discoverer"):
        print('"Proteome Discoverer" method being uses for peptides.')
        rename_dict = { 'Sequence': 'seq',
                        'Master Protein Accessions' : 'accession',
                       'Description': 'description',
                        'apQuant Area': 'intensity',
                       '# PSMs': 'PSM',
                       '# AAs': 'prot_len'}
        counter_i = 0  # zero a counter
        number_of_input_files = len(app.peptide_file_names)
        print(f"number of input files is: {number_of_input_files}")
        for sample_i, path_name_i in zip(app.peptide_file_names, app.peptide_files):
            print(path_name_i['path'])
            df_i = pd.read_excel(path_name_i['path'])
            # Truncate Column names
            df_i = df_i.rename(columns=lambda col: truncate_after_Area(col))
            df_i.rename(columns = rename_dict, inplace=True)
            df_i = df_i.reset_index()
            df_i = df_i[['seq', 'accession', 'intensity', 'PSM']]
            df_i = df_i.groupby('seq').agg(PSM = ('PSM', 'sum'), accession = ('accession', 'first'),
                                   intensity =('intensity', 'sum')).reset_index()
            df_i['sample'] = sample_i
            all_peptides_input_df = pd.concat([all_peptides_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        print("Done importing files")
        return all_peptides_input_df       
        
    elif(app.input_type == "MSFragger"):
        print('"MSFragger" Import Method being used for Peptides.')
        rename_dict = {'Peptide': 'seq', 
                       'Spectral Count': 'PSM', 
                       'Protein ID': 'accession',
                       'Intensity': 'intensity',
                       'Assigned Modifications': 'mods'}
        counter_i = 0  # zero a counter
        number_of_input_files = len(app.peptide_file_names)
        print(f"number of input files is: {number_of_input_files}")
        for sample_i, path_name_i in zip(app.peptide_file_names, app.peptide_files):
            print(path_name_i['path'])
            df_i = pd.read_csv(path_name_i['path'], sep='\t')
            df_i.rename(columns = rename_dict, inplace=True)
            df_i = df_i.reset_index()
            df_i = df_i[['seq', 'accession', 'intensity', 'PSM']]
            df_i = df_i.groupby('seq').agg(PSM = ('PSM', 'sum'), accession = ('accession', 'first'),
                                   intensity =('intensity', 'sum')).reset_index()
            df_i['sample'] = sample_i
            all_peptides_input_df = pd.concat([all_peptides_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        print("Done importing files")
        return all_peptides_input_df
    
    else:
        print("Method not recognized")
####### MaxQuant Results have no "mods Column in peptide results


# Calculating means and stds for groups
def calc_stats_by_group(df, property):
    groups = set(df['group'])
    print(groups)
    number_of_groups = len(groups)
    property_by_group = pd.DataFrame()
    for group in groups:
        # print(group)
        df_i = df[df['group'] == group]
        mean_i = df_i[property].mean()
        std_i = df_i[property].std()
        comment_i = "*std is a poor measure for this type of distribution"
        df_out_i = pd.DataFrame({'group': [group], 'mean': [mean_i], 'std': std_i, 'comment': comment_i})
        property_by_group = pd.concat([property_by_group, df_out_i], axis = 0)
    return(property_by_group)

def protein_import_function(app):
    all_proteins_input_df = pd.DataFrame()
    if (app.input_type == "Manually Created"):
        counter_i = 0  # zero a counter
        rename_dict = {'# PSMs': 'PSM', 'PSMs': 'PSM', 'MS/MS count': 'PSM',
               'Protein': 'accession', 'Accession': 'accession',
               'Description': 'description',
               '# AAs': 'prot_len',  '# Aas': 'prot_len', # pept_len is peptide length.
               '# Unique Peptides': 'numb_uniq_pept',
               'Area': 'intensity', 'area': 'intensity', 'Intensity': 'intensity'}
        print("Started Importing Protein Groups Files")
        number_of_input_files = len(app.protein_files)
        for sample_i, path_name_i in zip(app.protein_file_names, app.protein_files):
            print(f"Processing: {sample_i}")
            df_i = pd.read_csv(path_name_i)
            df_i.rename(columns=rename_dict, inplace=True)
            df_i = df_i[['accession', 'description', 'PSM', 'prot_len']]
            df_i['sample'] = sample_i
            all_proteins_input_df = pd.concat([all_proteins_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        return(all_proteins_input_df)
    elif (app.input_type == "MaxQuant"):
        counter_i = 0  # zero a counter
        rename_dict = {'Sequence length': 'prot_len', 'MS/MS count': 'PSM',
               'Protein IDs': 'accession', 'Protein names': 'description',  
               '# Unique Peptides': 'numb_uniq_pept'}
        print("Started Importing Protein Groups Files")
        number_of_input_files = len(app.protein_files)
        for sample_i, path_name_i in zip(app.protein_file_names, app.protein_files):
            print(f"Processing: {sample_i}")
            df_i = pd.read_csv(path_name_i, sep ='\t')
            df_i.rename(columns=rename_dict, inplace=True)
            print(df_i.columns)
            df_i = df_i[['accession', 'description', 'PSM', 'prot_len']]
            df_i['sample'] = sample_i
            all_proteins_input_df = pd.concat([all_proteins_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        return(all_proteins_input_df)
    
    elif (app.input_type == "Proteome Discoverer"):
        counter_i = 0  # zero a counter
        rename_dict = {'Description': 'description',
            'Accession': 'accession',
            '# PSMs': 'PSM', 
            '# AAs': 'prot_len'}
        print("Started Importing Protein Groups Files")
        number_of_input_files = len(app.protein_files)
        for sample_i, path_name_i in zip(app.protein_file_names, app.protein_files):
            print(f"Processing: {sample_i}")
            df_i = pd.read_excel(path_name_i)
            df_i.rename(columns=rename_dict, inplace=True)
            df_i = df_i[['accession', 'description', 'PSM', 'prot_len']]
            df_i['sample'] = sample_i
            all_proteins_input_df = pd.concat([all_proteins_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        return(all_proteins_input_df)
    elif (app.input_type == "MSFragger"):
        counter_i = 0  # zero a counter
        rename_dict = {'Protein Description': 'description',
                       'Protein ID': 'accession',
                       'Length': 'prot_len', 
                       'Total Spectral Count': 'PSM'}
        print("Started Importing Protein Groups Files")
        number_of_input_files = len(app.protein_files)
        for sample_i, path_name_i in zip(app.protein_file_names, app.protein_files):
            print(f"Processing: {sample_i}")
            df_i = pd.read_csv(path_name_i, sep ='\t')
            df_i.rename(columns=rename_dict, inplace=True)
            df_i = df_i[['accession', 'description', 'PSM', 'prot_len']]
            df_i['sample'] = sample_i
            all_proteins_input_df = pd.concat([all_proteins_input_df, df_i])
            counter_i = counter_i + 1
            if counter_i % 5 == 0:
                print(f"Processing file {counter_i} of {number_of_input_files}")
        return(all_proteins_input_df)
    else:
        print("Error in protein file processing!!!")









def protein_PSM_pivot_and_sort(protein_input_df):
    # Protein_input_df processing.  Pivot and sort decreasing order of abundance based on PSMs.
    # combining each protein into a single row.
    protein_PSM_df = protein_input_df.pivot(index=['accession', 'description', 'prot_len'], columns='sample', values='PSM')
    protein_PSM_df = protein_PSM_df.reset_index()
    protein_PSM_df = protein_PSM_df.fillna(0)
    protein_PSM_df['PSM'] =  protein_PSM_df.drop(columns = ['accession', 'description', 'prot_len']).sum(axis = 1)
    protein_PSM_df['PSM/AA'] = (protein_PSM_df['PSM'] + 0.01)/protein_PSM_df['prot_len']
    protein_PSM_df = protein_PSM_df.sort_values(by = 'PSM/AA', ascending = False)
    protein_PSM_df = protein_PSM_df.reset_index()
    protein_PSM_df = protein_PSM_df.drop('index', axis=1)  
    return(protein_PSM_df)


def create_group_dict(peptide_groups, peptide_file_names):
    """
    Adjusting/Correcting for non-entries and incorrect entries in peptide group boxes 
    and linking peptide_groups and peptide file names in "group_dict"
    """
    # Replace any blank value in list with 99
    peptide_groups = ['99' if item == '' else item for item in peptide_groups]
    unique_group_numbers = list(set(peptide_groups))
    number_of_groups = len(set(peptide_groups))
    # Create a dictionary with group_number as key and peptide files as values
    group_dict = {}
    # Iterate through the peptide_groups and peptide_file_names simultaneously using zip
    for group, file_name in zip(peptide_groups, peptide_file_names):
        # If the group number is not yet a key in the dictionary, add it with an empty list
        if group not in group_dict:
            group_dict[group] = []
        # Append the current name to the list associated with its group number
        group_dict[group].append(file_name)
    return(group_dict, number_of_groups, unique_group_numbers)


def calc_group_info_basic(peptide_file_names, protein_input_df, peptide_input_df):
    df_summary = pd.DataFrame()
    i = 0
    total_number_of_runs = len(peptide_file_names)
    print(f"Total number of runs to process: {total_number_of_runs}")
    for sample_i in peptide_file_names:
        protein_input_df_i = protein_input_df[protein_input_df['sample'] == sample_i]
        number_of_protein_groups = len(protein_input_df_i)
        peptide_input_df_i = peptide_input_df[peptide_input_df['sample'] == sample_i]
        number_of_peptides = len(peptide_input_df_i)
        number_of_PSM =  peptide_input_df_i['PSM'].sum()
        list_i =[sample_i, number_of_protein_groups, number_of_peptides, number_of_PSM] 
        columns=['sample', 'Protein Groups', 'Peptides', 'PSM']
        df_i = pd.DataFrame([list_i], columns=columns)
        df_summary = pd.concat([df_summary, df_i], ignore_index=True)
        i = i + 1
        print(f"Processing file {i} of {total_number_of_runs}")
    df_2 = df_summary.set_index('sample').T     
    return(df_2)


from scipy import stats

def calculate_mean(row, cols):
    return row[cols].mean()

def calculate_std(row, cols):
    return row[cols].std()

# Define a function to perform the t-test for each row
def calculate_p_value(row, group_1, group_2):
    t_statistic, p_value = stats.ttest_ind(row[group_1], row[group_2])
    return p_value


def getting_row_means_and_sds(df, group_dict, use_base_group =True):
    deep_copy_group_dict = copy.deepcopy(group_dict)
    unique_group_numbers = len(deep_copy_group_dict)
    key_list = list(deep_copy_group_dict.keys())
    df_out = pd.DataFrame()
    for key in key_list:
        columns_for_group_i = deep_copy_group_dict[key]
        mean_column_name = f"mean(Group {key})"
        df_out[mean_column_name] = df.apply(calculate_mean, axis=1, args=(columns_for_group_i,))

        # df[mean_column_name]  = df.apply(calculate_mean, axis=1, args=(columns_for_group_i, ))
        std_column_name = f"std(Group {key})"
        df_out[std_column_name] = df.apply(calculate_std, axis = 1, args = (columns_for_group_i, ))
    # Calculating p-values
    if len(deep_copy_group_dict) > 1:
        number_of_groups = len(deep_copy_group_dict)
        number_of_comparisons = number_of_groups -1
        print(f"Number of Comparisons is {number_of_comparisons}")
        for key in key_list:
            if (key == key_list[0]):
                print("skip first key")
            else:
                dict_for_modifying = copy.deepcopy(group_dict)
                columns_group_A = dict_for_modifying[key_list[0]]
                if use_base_group == False:
                    columns_group_A_modified = columns_group_A
                    del columns_group_A_modified[0]
                else: 
                    columns_group_A_modified = columns_group_A
                columns_group_B = deep_copy_group_dict[key]
                print(f'comparison group {key} versus base group {key_list[0]}.')
                col_name = f'P_value ({key} vs {key_list[0]})'
                df_out[col_name] = df.apply(calculate_p_value, axis=1, args=(columns_group_A_modified, columns_group_B))
                #df[col_name] = df.apply(calculate_p_value, axis=1, args=(columns_group_A, columns_group_B))
    df_out.index = df.index
    return(df_out)




def plots_for_protein_and_peptide_info(df, dict_1):
    # Create figure with higher DPI for sharpness
    fig, axs = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(14, 16),
        dpi=120,
        constrained_layout=True
    )

    colors = ['#2ca02c', '#1f77b4', '#9467bd', "#39caf2"]

    # -----------------------------
    # Protein Groups per sample
    # -----------------------------
    data_i = df.loc['Protein Groups']
    axs[0, 0].bar(
        data_i.index,
        data_i.values,
        color=colors[0],
        edgecolor='black',
        linewidth=1.0
    )
    axs[0, 0].set_title('Protein Groups per Sample', fontsize=12, weight='bold')
    axs[0, 0].set_xlabel('Sample / Run')
    axs[0, 0].set_ylabel('Protein Groups')
    axs[0, 0].tick_params(axis='x', rotation=90, labelsize=8)

    # -----------------------------
    # Peptides per sample
    # -----------------------------
    data_i = df.loc['Peptides']
    axs[1, 0].bar(
        data_i.index,
        data_i.values,
        color=colors[1],
        edgecolor='black',
        linewidth=1.0
    )
    axs[1, 0].set_title('Peptides per Sample', fontsize=12, weight='bold')
    axs[1, 0].set_xlabel('Sample / Run')
    axs[1, 0].set_ylabel('Peptides')
    axs[1, 0].tick_params(axis='x', rotation=90, labelsize=8)

    # -----------------------------
    # Protein Groups per group
    # -----------------------------
    df_pg = dict_1['Protein Groups']
    axs[0, 1].bar(
        df_pg.index,
        df_pg['mean'],
        yerr=df_pg['sd'],
        capsize=6,
        color=colors[2],
        edgecolor='black',
        linewidth=1.0,
        error_kw=dict(elinewidth=1.2, capthick=1.2)
    )
    axs[0, 1].set_title('Protein Groups per Group', fontsize=12, weight='bold')
    axs[0, 1].set_xlabel('Group')
    axs[0, 1].set_ylabel('Protein Groups')

    # -----------------------------
    # Peptides per group
    # -----------------------------
    df_pep = dict_1['Peptides']
    axs[1, 1].bar(
        df_pep.index,
        df_pep['mean'],
        yerr=df_pep['sd'],
        capsize=6,
        color=colors[3],
        edgecolor='black',
        linewidth=1.0,
        error_kw=dict(elinewidth=1.2, capthick=1.2)
    )
    axs[1, 1].set_title('Peptides per Group', fontsize=12, weight='bold')
    axs[1, 1].set_xlabel('Group')
    axs[1, 1].set_ylabel('Peptides')

    # -----------------------------
    # Overall title
    # -----------------------------
    fig.suptitle(
        'Summary Statistics per Sample and Group\n(Protein Groups and Peptides)',
        fontsize=16,
        weight='bold'
    )

    return fig


def plots_for_PSM_info(df, dict_1):
    # Create figure with higher DPI for sharpness
    fig, axs = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(10, 6),
        dpi=120,
        constrained_layout=True
    )

    colors = ["#1818f8", "red", '#9467bd', "#39caf2"]

    # -----------------------------
    # Protein Groups per sample
    # -----------------------------
    data_i = df.loc['PSM']
    axs[0].bar(
        data_i.index,
        data_i.values,
        color=colors[0],
        edgecolor='black',
        linewidth=1.0
    )
    axs[0].set_title('PSM per Sample', fontsize=12, weight='bold')
    axs[0].set_xlabel('Sample / Run')
    axs[0].set_ylabel('PSMs')
    axs[0].tick_params(axis='x', rotation=90, labelsize=8)


    # -----------------------------
    # PSMs per group
    # -----------------------------
    df_pg = dict_1['PSM']
    axs[1].bar(
        df_pg.index,
        df_pg['mean'],
        yerr=df_pg['sd'],
        capsize=6,
        color=colors[2],
        edgecolor='black',
        linewidth=1.0,
        error_kw=dict(elinewidth=1.2, capthick=1.2)
    )
    axs[1].set_title('PSMs per Group', fontsize=12, weight='bold')
    axs[1].set_xlabel('Group')
    axs[1].set_ylabel('PSMs')

    
    # -----------------------------
    # Overall title
    # -----------------------------
    fig.suptitle(
        'Summary Statistics per Sample and Group\n(PSMs)',
        fontsize=16,
        weight='bold'
    )

    return fig



def coverage_df_by_AA(POI_seq, df):
    df_calc_1 = pd.DataFrame()
    df_calc_1['AA_pos'] = list(range(len(POI_seq)))  # AA_pos is AA position
    df_calc_1['AA'] = list(POI_seq)
    for i in df.columns:
        df_calc_1[i] = 0.0 
        for seq_i, value_i in zip(df.index.get_level_values('seq'), df[i]):
            # Find the location of seq in the protein
            index = POI_seq.find(seq_i)
            if index == -1:  # This means that the peptide cannot be found in protein
                continue  # Skips to next peptide
            for j in range(0, len(seq_i)):
                df_calc_1.loc[index + j, i] = df_calc_1.loc[index + j, i] + value_i
                # df_calc_1.loc[index + j, i] += value_i
            # Get relative coverage ('rel_cov')
        max_num_value = df_calc_1[i].max()
        min_num_value = df_calc_1[i].min()
        if max_num_value == min_num_value:
            print("min and max are the same!")
            max_num_value = min_num_value + 1
        if min_num_value < 0:
            df_calc_1[i] = df_calc_1[i] - min_num_value
            df_calc_1[i] = df_calc_1[i] / (max_num_value - min_num_value)
        else:
            df_calc_1[i] = df_calc_1[i]/max_num_value  
    # This creates a column ['rand'] that is a randomization of the column for the first sample 
    df_calc_1['rand'] = df_calc_1[df_calc_1.columns[2]].sample(frac=1, random_state=42).reset_index(drop=True)
    return df_calc_1








##################################

def make_coverage_plots(coverage_df, app, POI_dict):
    number_of_runs = len(app.peptide_file_names) + 1  # 1 is added is for 'rand'
    print("Point B2", number_of_runs)
    runs = app.peptide_file_names
    if (number_of_runs > 20):
        print("Point D2")
        desired_runs = app.peptide_file_names[0:19]
        desired_runs.append('rand')
        messagebox.showinfo("Alert", f"The number of files exceeds amount allowed for visualization")
        print("Maximum Number of Comparisons is 20.  Only first 20 used.")
    else:
        desired_runs = runs
    print("point A2", desired_runs)

    rectangles_needed = len(desired_runs)
    # Create a figure and axes object
    fig, ax = plt.subplots(figsize=(8, (number_of_runs + 1)))
    # Calculate length of POI (protein of interest)
    POI_len = coverage_df.shape[0]
    print(f"Length Protein: {POI_len} Amino Acids") 
    print(desired_runs)
    for i in range(0, len(desired_runs)):
        # Calculate Similarity index versus first run which should be base run
        #similarity_index_i = compute_similarity_index(df, desired_runs[i], base_run, random_sample)
        similarity_index_i = compute_similarity_index(coverage_df, desired_runs[i], 
                                                      app.peptide_file_names[0], 'rand')
        box_height = (1.0 / (len(desired_runs) + 1))/2
        rect_name = Rectangle((0.0, (1.0-((i+1)*box_height*2))), 0.99, box_height, facecolor='white', 
                              edgecolor='black', linewidth = 2, fill = False)

        # Add the outline to the axes
        ax.add_patch(rect_name)
        ax.text(0.5, (1.0-((i+1.2)*box_height*2)), desired_runs[i], color='black', fontsize=9, ha='center', va='bottom')
        ax.text(1.05, (1.0-((i+0.8)*box_height*2)), f"SI = {similarity_index_i:.2f}", color='black', fontsize=9, ha='center', va='bottom')
        column_name = desired_runs[i]
        for index in coverage_df.index:
            alpha_1 = coverage_df.loc[index, column_name]   # transparency
            # print(index, alpha_1)
            rect_1B = Rectangle((((index /POI_len)*0.99), (1.0-((i+1)*box_height*2))), 0.95/(POI_len), box_height, 
                                facecolor = 'blue', alpha = alpha_1, edgecolor = 'none')
            ax.add_patch(rect_1B)
            ax.add_patch(rect_name)
    
    # Remove axes and borders
    ax.set_axis_off() # Hides the x and y axes
    ax.set_frame_on(False) # Hides the border of the plot

    max_i = len(desired_runs)  # max_i is number of runs.  This would be the index for the next run.
    rect_scale = Rectangle((0.0, (1.0-((max_i + 1)*box_height*2))), 0.99, box_height, facecolor='none', 
                        edgecolor='black', linewidth = 2)

    # Redoing Boxes for better image quality
    ax.add_patch(rect_scale)
    
    # Adding some additional information at top
    ax.text(0.0, 1.1, f"Accession Number is: {app.accession}", color = 'black', fontsize = 10, ha = 'left', 
            va = 'bottom')
    POI_description = POI_dict['description']
    POI_name = POI_description.split(" OS=")[0]
    ax.text(0.0, 1.05, f"Protein Description: {POI_name}", color = 'black', fontsize = 10, ha = 'left', 
            va = 'bottom', wrap=True)
    ax.text(0.0, 1.00, f"Calculation Method: {app.calculation_type}", color = 'black', fontsize = 10, ha = 'left', 
            va = 'bottom', wrap=True)

    # Adding scale at bottom
    #Label for scale:
    if app.calculation_type == "PSM":
        label_for_scale = "Number of PSMs (Normalized)"
    elif app.calculation_type == "intensity":
        label_for_scale = "Ion Intensity (Normalized)"
    elif app.calculation_type == "log_2_PSM":
        label_for_scale = "log-2 of PSM (Normalized)"
    else:
        label_for_scale = "log-2 of Ion Intensity (Normalized)"

    ax.text(0.5, (1.0-((max_i + 1.3)*box_height*2)), label_for_scale, color='black', fontsize=9, 
            ha='center', va='bottom')
    ax.text(0, (1.0-((max_i + 1.3)*box_height*2)), '0.00', color = 'black', fontsize = 8, ha = 'center', va = 'bottom')
    ax.text(0.99, (1.0-((max_i + 1.3)*box_height*2)), '1.00', color = 'black', fontsize = 8, ha = 'center', va = 'bottom')
    for i_3 in range(0, 99):
        alpha_scale = i_3/100   # transparency
        rect_scale_B = Rectangle((0.0 + (0.01 * i_3), (1.0-((max_i + 1)*box_height*2))), 0.01, box_height, facecolor='blue', 
                                 alpha = alpha_scale, edgecolor='none')
        ax.add_patch(rect_scale_B)
    return(fig)


def compute_similarity_index(df, run_i, base_run, random_sample):
    POI_len = df.shape[0]
    Nonscaled_similarity = 1 - (abs(df[run_i] - df[base_run])).sum()/POI_len      
    random_similarity = 1 - (abs(df[random_sample] - df[base_run])).sum()/POI_len
    scaled_similarity = (Nonscaled_similarity - random_similarity)/(1.0000000001-random_similarity) 
    return(scaled_similarity)


def processing_alternatives(protein_PSM_df, app): #processing_method, accession, rank):
    if (app.processing_method == "accession_input"):  # Option 1 is for User provided Accession Number
        list_of_accession_numbers = list(protein_PSM_df['accession'])
        if app.accession not in list_of_accession_numbers:
            messagebox.showerror("File error", f"{app.accession} was not found in protein results")
        print("Processing method is protein accession number")        
        POI = app.accession
        return(POI)

    elif (app.processing_method == "random"):  # Option 2 is Random Protein
        # To select at random protein, use the following code: 
        # Select a random protein from list of top 500 most abundant proteins based on spectral count
        max_rank = 500
        if max_rank > len(protein_PSM_df):
            max_rank = len(protein_PSM_df)
        random_protein_rank = random.randint(0, max_rank) - 1
        POI = protein_PSM_df.iloc[random_protein_rank]['accession']
        return POI
    
    elif (app.processing_method == "rank_input"):  # Option 3 is Rank based
        # For selection based on abundance, simply use index.  Proteins were sorted based on the number of
        # PSMs/(length of protein in number of AAs) in descending order.  Therefore, the most abundant protein 
        # is in the first row (index value = 0).
        if not bool(app.rank):
            app.rank = 10 # Sets default rank to 10
        print(f"rank is {int(app.rank)}")
        protein_index = int(app.rank)
        POI = protein_PSM_df.iloc[(protein_index - 1)] ['accession'] 
        print(f"Processing method is rank, POI is: {POI}")
        return POI
   
    elif (app.processing_method == "top_N"):
        print("Processing method is top_N so there is no protein of interest (POI)")
    else:
        print("Something went wrong with processing method")
        sys.exit(1)




def calculation_alternatives(POI_df, app):
    if app.calculation_type == "PSM":
        pivoted_POI_df = POI_df.pivot(index=['seq'], columns='sample', values='PSM').fillna(0)
    
    elif app.calculation_type == "intensity":
        pivoted_POI_df = POI_df.pivot(index=['seq'], columns='sample', values='intensity').fillna(0)
    
    elif app.calculation_type == "log_2_PSM":
        pivoted_POI_df = POI_df.pivot(index=['seq'], columns='sample', values='PSM').fillna(0)
        # Add a small constant to avoid log of zero and then calculate log-2.
        cols_to_update = pivoted_POI_df.columns[pivoted_POI_df.columns != 'seq']
        pivoted_POI_df[cols_to_update] = np.log2(pivoted_POI_df[cols_to_update] + 0.1)
    
    elif app.calculation_type == "log_2_intensity":
        pivoted_POI_df = POI_df.pivot(index=['seq'], columns='sample', values='intensity').fillna(0)
        # Add a small constant to avoid log of zero and then calculate log-2.
        cols_to_update = pivoted_POI_df.columns[pivoted_POI_df.columns != 'seq']
        pivoted_POI_df[cols_to_update] = np.log2(pivoted_POI_df[cols_to_update] + 10)

    else:
        raise ValueError(f"Unknown calculation type: {app.calculation_type}")
    
    # An error occurs if I have a sample with no peptides for a given protein.  This solves it.
    if(set(POI_df['sample']) != set(app.peptide_file_names)):
        print ("Lists are not the same")
        missing_elements = list(set(app.peptide_file_names) - set(POI_df['sample']))
        print(f"length of missing elements is: {len(missing_elements)}")

        print(f"Elements in list_b but not in list_a: {missing_elements}")
        for element in missing_elements:
            pivoted_POI_df[element] = 0
    return(pivoted_POI_df) 


def create_similarity_matrix(df):
    POI_len = len(df)
    df_calc = df
    runs_to_compare = df.columns.drop(['AA_pos', 'AA', 'rand'])
    
    # Define the dimensions of the DataFrames
    number_of_rows = len(runs_to_compare)
    number_of_columns = len(runs_to_compare)   

    # Create a NumPy array filled with ones for the similarity index. 
    np_similarity = np.ones((number_of_rows, number_of_columns))

    # Convert the NumPy arrays to a Pandas DataFrame
    df_similarity = pd.DataFrame(np_similarity)

    # Label the rows and columns with sample information.
    df_similarity.columns = runs_to_compare
    df_similarity.index = pd.Index(runs_to_compare, name='sample')
    remaining_runs = list(runs_to_compare.copy())
    for run_x in runs_to_compare:
        remaining_runs.pop(0)  # removes the first item from a list.
        for run_y in remaining_runs:
            #print(run_x, run_y)
            x_y_similarity = 1 - (abs(df_calc[run_y] - df_calc[run_x])).sum()/POI_len
            random_similarity = 1 - (abs(df_calc['rand'] - df_calc[run_x])).sum()/POI_len
            scaled_x_y_similarity = (x_y_similarity - random_similarity)/(1.0000000001-random_similarity)
            #print(scaled_x_y_similarity)
            df_similarity.loc[run_x, run_y] = scaled_x_y_similarity
            df_similarity.loc[run_y, run_x] = scaled_x_y_similarity 
    return(df_similarity)


def save_figure_unique(fig, base_filename, directory=".", format="png"):
    """
    Saves a Matplotlib figure to a file, adding a counter to the filename
    if a file with the same name already exists.

    Args:
        fig (matplotlib.figure.Figure): The Matplotlib figure to save.
        base_filename (str): The base name for the file (without extension).
        directory (str): The directory to save the file in. Defaults to the current directory.
        format (str): The file format (e.g., "png", "jpg", "pdf"). Defaults to "png".
    Created by Google AI
    """
    counter = 0
    while True:
        if counter == 0:
            filename = f"{base_filename}.{format}"
        else:
            filename = f"{base_filename}_{counter}.{format}"

        filepath = os.path.join(directory, filename)

        if not os.path.exists(filepath):
            fig.savefig(filepath, format=format)
            print(f"Figure saved to: {filepath}")
            break
        counter += 1



def save_dataframe_with_counter(df, base_filename, file_format='csv', index=False):
    """
    Saves a Pandas DataFrame to a file, adding a counter to the filename 
    if a file with the same name already exists.

    Args:
        df (pd.DataFrame): The DataFrame to save.
        base_filename (str): The desired base filename (without extension).
        file_format (str): The desired file format ('csv', 'xlsx', 'parquet', etc.).
        index (bool): Whether to write the DataFrame index to the file.
    Created by Google AI
    """
    counter = 0
    while True:
        if counter == 0:
            filename = f"{base_filename}.{file_format}"
        else:
            filename = f"{base_filename}_{counter}.{file_format}"

        if not os.path.exists(filename):
            if file_format == 'csv':
                df.to_csv(filename, index=True)
            elif file_format == 'xlsx':
                df.to_excel(filename, index=True)
            elif file_format == 'parquet':
                df.to_parquet(filename, index=True)
            # Add more formats as needed
            print(f"DataFrame saved to: {filename}")
            break
        else:
            counter += 1



def show_dataframe_in_new_window_treeview(dataframe):
    # This function creates a new Tkinter window and displays the given DataFrame in a 
    # Treeview widget.
    # This function was created by Google AI and then modified.
    new_window = tk.Toplevel()
    new_window.title("Overall Similarity Matrix")

    tree = ttk.Treeview(new_window)
    tree.pack(fill=tk.BOTH, expand=True)

    # Define columns
    columns = list(dataframe.columns)
    tree["columns"] = columns
    #tree["show"] = "headings" # Hide the default first column

    # Set headings
    for col in columns:
        tree.heading(col, text=col)
        tree.column(col, width=100, anchor="w") # Adjust width as needed

    # Insert data
    for index, row in dataframe.iterrows():
        tree.insert("", "end", values=list(row))


def plot_similarity_table(df, app, title="Overall Similarity Index", figsize=(8, 4),
    font_size=10, scale=(1.2, 1.2), header_color="#c7d4ff", column_widths=None):
    # list of widths for each column
    # Made with the help of Chat GPT
    """
    Plots a DataFrame as a formatted table with 3 significant digits,
    colored headers, and adjustable column widths.

    Parameters:
        df : pandas DataFrame
        app.N_value : int for number of proteins used for comparision.
        app.calculation_type : str for method used ("PSM", "log_2_PSM", "intensity", "log_2_intensity"
        title : str
        figsize : tuple
        font_size : int
        scale : tuple
        header_color : str (Hex color)
        column_widths : list or None (width per column)
    """

    # Reset index for display
    df = df.reset_index()

    # Format all numeric entries to 3 significant digits
    df = df.map(lambda x: f"{x:.3g}" if isinstance(x, (int, float, np.number)) else x)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('off')
    ax.axis('tight')

    # Create table
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        loc='center',
        cellLoc='center'
    )

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    table.scale(*scale)

    # === Header Styling ===
    for (row, col), cell in table.get_celld().items():
        if row == 0:  # header row
            cell.set_facecolor(header_color)
            cell.set_text_props(weight='bold')  

    # === Column Width Adjustment ===
    if column_widths is not None:
        for col_index, width in enumerate(column_widths):
            table.auto_set_column_width(col=col_index)
            for row_index in range(len(df) + 1):  # +1 for header
                table[(row_index, col_index)].set_width(width)

    # Add title and descriptors
    ax.text(0.5, 1.10, title, horizontalalignment='center',
            verticalalignment='bottom', transform=ax.transAxes,
            fontsize=12, color='black', fontweight='bold')

    ax.text(0.0, 1.0, f'Top N proteins ({app.N_value} proteins)',
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes, fontsize=9)

    ax.text(0.0, 0.95, f'Calculation Type = {app.calculation_type}',
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes, fontsize=9)
    ax.text(0.0, 0.90, "This table is saved in results folder as a CSV file.",
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes, fontsize=9)

    fig.tight_layout()





def summarize_by_group(
    df: pd.DataFrame,
    group_dict: dict,
    ref_group=None,
    assume_equal_variance=True,
    exclude_columns=None
):
    """
    Compute mean, SD, and p-values vs reference group.
    Output: one table per property with groups as rows.

    If a t-test cannot be performed, p-value is set to 1.0.
    If SD cannot be computed (fewer than 2 values), SD is set to 0.0.
    """

    if exclude_columns is None:
        exclude_columns = set()
    else:
        exclude_columns = set(exclude_columns)

    if ref_group is None:
        ref_group = list(group_dict.keys())[0]

    results = {}

    for prop in df.index:
        rows = []

        # Reference values
        ref_cols = [c for c in group_dict[ref_group] if c not in exclude_columns]
        ref_vals = df.loc[prop, ref_cols].dropna()

        for group, cols in group_dict.items():
            use_cols = [c for c in cols if c not in exclude_columns]
            values = df.loc[prop, use_cols].dropna()

            mean = values.mean()

            # SD: set to 0.0 if not enough values
            if len(values) >= 2:
                sd = values.std(ddof=1)
            else:
                sd = 0.0

            # Default p-value
            pval = 1.0

            if group != ref_group:
                try:
                    if (
                        len(ref_vals) >= 2
                        and len(values) >= 2
                        and ref_vals.var(ddof=1) > 0
                        and values.var(ddof=1) > 0
                    ):
                        _, pval = ttest_ind(
                            ref_vals,
                            values,
                            equal_var=assume_equal_variance,
                        )
                except Exception:
                    pval = 1.0

            rows.append(
                {
                    "group": group,
                    "mean": mean,
                    "sd": sd,
                    f"p_value_vs_{ref_group}": pval,
                }
            )

        results[prop] = pd.DataFrame(rows).set_index("group")

    return results




def plot_runs_and_stats(df, df_2, property):
    # Create a figure and two subplots side-by-side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # --- FIRST PLOT ---
    sample_names = df['sample']
    values_1 = df[property]

    ax1.bar(sample_names, values_1, color='blue')
    ax1.tick_params(axis='x', rotation=90)
    ax1.set_title(f"{property} by sample")
    ax1.set_ylabel(property)


    # --- SECOND PLOT (error bars) ---
    df_2_indexed = df_2.set_index('group')

    df_2_indexed.plot.bar(
        y='mean',
        yerr='sd',
        ax=ax2,
        capsize=4,
        rot = 0,  # Don't rotate
        legend = False,
        title = f'{property} Chart with Error Bars*',
        ylabel = property
    )
    ax2.text(
        0.5, -0.15,  # x, y in axis coordinates
        "* For calculation of mean and sd, base run results have been removed",
        ha='center',  # horizontal alignment
        va='top',     # vertical alignment
        fontsize=10,
        transform=ax2.transAxes  # coordinates relative to the axis
    )
    plt.tight_layout()
    return(fig)



def calc_aver_property_per_seq(seq, dict):
    # Function to Calculate the average property of a peptide sequence
    # using AA values from AAindex.
    
    # Args:
    #    seq (str): Peptide sequence using one-letter amino acid codes.
    #    dictionary for the various properties.
    # Returns: float: Average value for that property
    values = []
    seq = seq.upper()
    for aa in seq:
        if aa in dict:
            values.append(dict[aa.upper()])
        else:
            raise ValueError(f"Unknown amino acid: {aa}")
    return sum(values) / len(values)


 
# Creates the AA_index_df for any property
def AA_index_df(app, POIs_PSM_df, group_dict, property_dict):
    XX_calc_df = POIs_PSM_df
    # Apply a lambda function to the 'numbers' column to square each value
    XX_calc_df['property_X'] = \
        XX_calc_df['seq'].apply(lambda x: calc_aver_property_per_seq(x, property_dict))
    # Add group column
    # Create a group column to dataframe
    sample_to_group = {sample: group
        for group, sample_list in group_dict.items()
        for sample in sample_list}
    XX_calc_df['group'] = XX_calc_df['sample'].map(sample_to_group)
            
    return (XX_calc_df)

# converts from peptide per row to PSM per row

def peptide_to_PSM_per_row(df):
    return (df.loc[df.index.repeat(df['PSM'])]   # repeat rows
        .drop(columns='PSM')               # drop PSM column
        .reset_index(drop=True)
    )


def AA_index_calc_function(app, AA_index_reference, AA_index_XXX, POIs_PSM_df, 
                          group_dict):
    value_for_test =getattr(app, AA_index_XXX)
    print(value_for_test, "Point L")
    if getattr(app, AA_index_XXX) == 1:
        X_property_description = aaindex1[AA_index_reference].description.split(' (', 1)[0]
        print(X_property_description)
        X_property_dict = aaindex1[AA_index_reference].values
        print(X_property_dict)
        X_property_cite = aaindex1[AA_index_reference].description.split(' (', 1)[1].split(')')[0]
        print(X_property_cite)

        X_property_df = AA_index_df(app, POIs_PSM_df, group_dict, X_property_dict)
        print(X_property_df)

        # Calculate stats for PSM_x_property df
        X_property_stats_df = calc_stats_by_group(X_property_df, "property_X")
        print(f"X_property_df for {AA_index_XXX}:")
        print()
        print(X_property_stats_df)

        property_x_plot_by_run = plot_property_by_run(X_property_df, "property_X", group_dict, 
                                AA_index = True, property_description = X_property_description)

        property_x_plot_by_group = plot_property_by_group(X_property_df, "property_X", group_dict, 
                                AA_index = True, property_description = X_property_description)

        # Save the all the data and the figures
        AA_index_property_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            '_' + AA_index_reference 
        save_dataframe_with_counter(X_property_df, AA_index_property_df_path_name, file_format='csv')
        print(AA_index_property_df_path_name)

        AA_index_property_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            '_' + AA_index_reference + "_stats"
        save_dataframe_with_counter(X_property_stats_df, AA_index_property_stats_df_path_name, file_format='csv')
        print(AA_index_property_stats_df_path_name)

        # Save the AA_index_property_by_run_figure
        save_figure_unique(property_x_plot_by_run, app.results_file_name + '_' + AA_index_reference + '_per_run', 
                       directory=app.results_subfolder_path, format="png")


        # Save the AA_index_property_by_group_figure
        save_figure_unique(property_x_plot_by_group, app.results_file_name + '_' + AA_index_reference + '_per_group', 
                       directory=app.results_subfolder_path, format="png")
    
    else: 
        X_property_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        X_property_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        property_x_plot_by_run =  make_blank_plot()
        property_x_plot_by_group =  make_blank_plot()

    return(X_property_df, X_property_stats_df, 
           property_x_plot_by_run, property_x_plot_by_group)




    



# Make blank plot for unrequested plots
def make_blank_plot():
    fig, ax = plt.subplots(figsize=(4, 2))
    ax.text(0.5, 0.5, 'This plot was not requested!', horizontalalignment='center', 
            verticalalignment='center', fontsize=20, color='black', transform=ax.transAxes)
    ax.axis('off')
    return fig





def backup_AA_index_calc_function(app, AA_index_reference, AA_index_XXX, POIs_peptide_df, 
                          group_dict):
    X_property_description = aaindex1[AA_index_reference].description.split(' (', 1)[0]
    print(X_property_description)
    X_property_dict = aaindex1[AA_index_reference].values
    print(X_property_dict)
    X_property_cite = aaindex1[AA_index_reference].description.split(' (', 1)[1].split(')')[0]
    print(X_property_cite)

    X_property_df = AA_index_df(app, POIs_peptide_df, group_dict, X_property_dict)
    print(X_property_df)

    PSM_X_property_df = peptide_to_PSM_per_row(X_property_df)
    print(f"Shape of PSM_X_property_df = {len(PSM_X_property_df)}")

    # Calculate stats for PSM_x_property df
    PSM_X_property_stats_df = calc_stats_by_group(PSM_X_property_df, "property_X")
    print(f"PSM_X_property_df for {AA_index_XXX}:")
    print()
    print(PSM_X_property_stats_df)

    property_x_plot_by_run = plot_property_by_run(PSM_X_property_df, "property_X", group_dict, 
                                AA_index = True, property_description = X_property_description)

    property_x_plot_by_group = plot_property_by_group(PSM_X_property_df, "property_X", group_dict, 
                                AA_index = True, property_description = X_property_description)

    # Save the all the data and the figures
    AA_index_property_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            '_' + AA_index_XXX 
    save_dataframe_with_counter(X_property_df, AA_index_property_df_path_name, file_format='csv')
    print(AA_index_property_df_path_name)

    AA_index_property_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            '_' + AA_index_XXX + "_stats"
    save_dataframe_with_counter(PSM_X_property_stats_df, AA_index_property_stats_df_path_name, file_format='csv')
    print(AA_index_property_stats_df_path_name)

    # Save the AA_index_property_by_run_figure
    save_figure_unique(property_x_plot_by_run, app.results_file_name + '_' + AA_index_XXX + '_per_run', 
                       directory=app.results_subfolder_path, format="png")


    # Save the AA_index_property_by_group_figure
    save_figure_unique(property_x_plot_by_group, app.results_file_name + '_' + AA_index_XXX + '_per_group', 
                       directory=app.results_subfolder_path, format="png")
    return(X_property_df, PSM_X_property_stats_df, 
           property_x_plot_by_run, property_x_plot_by_group)

def create_peptide_length_df(POIs_PSM_df, group_dict): 
    # invert the group_dict
    sample_to_group = {sample: group for group, samples in group_dict.items() for sample in samples}
    POIs_PSM_df["group"] = POIs_PSM_df["sample"].map(sample_to_group)

    pept_len_df = POIs_PSM_df.reset_index()
    
    pept_len_df['pept_len'] = pept_len_df['seq'].str.len()
    PSM_pept_len_df = (pept_len_df.groupby(["group", "sample", "pept_len"]).size().
        reset_index(name="count"))
    return(PSM_pept_len_df)

import matplotlib.pyplot as plt
import seaborn as sns

def plot_peptide_length_by_run(peptide_length_df, group_dict):
    fig, ax = plt.subplots(figsize=(8, 5))

    sns.barplot(data=peptide_length_df, x="pept_len", y="count", hue="sample", ax=ax )

    ax.set_xlabel("Peptide Length")
    ax.set_ylabel("Number of Peptides")
    ax.set_title("Peptide Length Distribution by Sample")

    fig.tight_layout()
    return fig




def plot_peptide_length_by_group(peptide_length_df, group_dict):
    fig, ax = plt.subplots(figsize=(8, 5))

    # Get unique groups and peptide lengths
    groups = peptide_length_df["group"].unique()
    pept_lens = peptide_length_df["pept_len"].unique()

    # Draw barplot WITHOUT seaborn error bars
    sns.barplot(
        data=peptide_length_df,
        x="pept_len",
        y="mean",
        hue="group",
        errorbar=None,
        ax=ax
    )

    # --- Add SD error bars manually ---
    n_groups = len(groups)
    bar_width = 0.8 / n_groups  # seaborn default total width ≈ 0.8

    for i, group in enumerate(groups):
        group_df = peptide_length_df[peptide_length_df["group"] == group]

        x_positions = (
            np.arange(len(pept_lens))
            - 0.4
            + bar_width / 2
            + i * bar_width
        )

        ax.errorbar(
            x=x_positions,
            y=group_df["mean"],
            yerr=group_df["std"],
            fmt="none",
            ecolor="black",
            capsize=2,
            linewidth=0.6
        )

    ax.set_xlabel("Peptide Length")
    ax.set_ylabel("Number of Peptides")
    ax.set_title("Peptide Length Distribution by Group")

    fig.tight_layout()
    return fig




def plot_property_by_run(df, property, group_dict, AA_index = False, property_description = None):
    # Create a figure and two subplots side-by-side (1 row, 2 columns)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6)) # figsize adjusts the overall size
    
    # Plot on the first subplot (ax1)
    # Create the box plot for various runs/samples
    samples = get_all_values(group_dict)

    df_kde = df[df['sample'].isin(samples)]
    
    # samples = list(set(df['sample']))   
    number_of_samples = len(samples)
    list_data = []
    for sample in samples:
        df_i = df[df['sample'] == sample]
        data_i = list(df_i[property])
        list_data.append(data_i)

    # Plot on the first subplot (ax1)
    # Create the bar chart for various runs/samples
    ax1.boxplot(list_data, tick_labels = samples)
    
    # Add labels, title, horizontal line, and adjust y-scale.
    ax1.set_xlabel('Sample/Run')
    if AA_index == False:
        ax1.set_ylabel(property.title())
        if property == "pI":
            ax1.set_title("Isoelectric Point (pI) per Sample")    
        else:
            ax1.set_title(f'{property.title()} per Sample')
    else:
        ax1.set_ylabel(property_description)
        ax1.set_title(f'{property_description} per Sample')
    
    ax1.tick_params(axis='x', labelrotation=90)
    
    # Create the density plot as a second subplot on ax2
    sns.kdeplot(data=df_kde, x=property, hue='sample', fill=True, common_norm=False, alpha=0.2, ax=ax2)

    ax2.set_title('Density Plot of Value by Sample/Run')
    if AA_index == False:
        ax2.set_xlabel(property.title())
    else:
        ax2.set_xlabel(property_description)
    ax2.set_ylabel('Density')
    #ax2.legend()

    # Adjust layout to prevent overlapping titles/labels
    fig.tight_layout()
    return(fig)


def get_all_values(d):
    # This function is directly from Google AI
    all_values = []
    for value in d.values():
        if isinstance(value, dict):
            # Recursively call the function for nested dictionaries
            all_values.extend(get_all_values(value))
        elif isinstance(value, list):
            # Iterate through list items and check for further nesting
            for item in value:
                if isinstance(item, dict):
                    all_values.extend(get_all_values(item))
                elif isinstance(item, list):
                    # For lists within lists, you would need another layer of handling 
                    # or a more robust general flattening utility (like one using 
                    # itertools.chain.from_iterable recursively for all iterables)
                    all_values.append(item) # Appends the inner list as is
                else:
                    all_values.append(item)
        else:
            all_values.append(value)
    return all_values



def plot_property_by_group(df, property, group_dict, AA_index = False, property_description = None):
    # Create a figure and two subplots side-by-side (1 row, 2 columns)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6)) # figsize adjusts the overall size
    groups = group_dict.keys()
    number_of_groups = len(groups)
    list_data = []
    list_group = []
    for group in groups:
        df_i = df[df['group'] == group]
        data_i = list(df_i[property])
        list_data.append(data_i)
        list_group.append(group)

    # Plot on the first subplot (ax1)
    # Create the bar chart for various runs/samples
    ax1.boxplot(list_data, tick_labels = list_group)
    #y_tick_locations = [-1, 0, 1]
    ##y_tick_labels = ['Min Value', 'Zero', 'Max Value']

    # Set the y-ticks and their labels
    #plt.yticks(y_tick_locations, y_tick_labels)


    # Add labels, title, horizontal line, and adjust y-scale.
    ax1.set_xlabel('Group')
    if AA_index == False:
        ax1.set_ylabel(property)
        ax1.set_title(f'{property} per Group')
    else:
        ax1.set_ylabel(property_description)
        ax1.set_title(f'{property_description} per Group')
    ax1.tick_params(axis='x', labelrotation=90)


    # Create the density plot as a second subplot on ax2
    group_order = list(group_dict.keys())
    sns.kdeplot(data = df, x = property, hue='group', hue_order=group_order, fill=True, 
                common_norm=False, alpha=0.2, ax = ax2)
    ax2.set_title('Density Plot Per Group')
    if AA_index == False:
        ax2.set_xlabel(property.title())
    else:
        ax2.set_xlabel(property_description)
    ax2.set_ylabel('Density')
    #ax2.legend(loc='upper right')

    # Adjust layout to prevent overlapping titles/labels
    fig.tight_layout()

    return(fig)

