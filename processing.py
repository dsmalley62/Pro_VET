import os
from pprint import pp
from tkinter import messagebox
import pandas as pd
import numpy as np
from utils import *
from processing_top_N import top_N_processing_pipeline
from processing_single_protein import single_protein_processing_pipeline

import copy


# -------------------------
# FASTA FILE NAME EXTRACTION
# -------------------------
def getting_fasta_file_names(app):
    fasta_file_names = [
        long_name.rsplit("/", 1)[-1]
        for long_name in app.fasta_files
    ]
    app.fasta_file_names = fasta_file_names
    return app


# --------------------------------
# PEPTIDE FILE NAME EXTRACTION
# --------------------------------
def get_peptide_file_names(app):
    peptide_file_names = []
    peptide_groups = []

    suffix_map = {
        "Manually Created": "_peptides.csv",
        "MaxQuant": "_peptides.txt",
        "Proteome Discoverer": "_PeptideGroups.xlsx",
        "MSFragger": "_peptide.tsv"
    }

    if app.input_type not in suffix_map:
        print("Error: Unknown input type for peptide files.")
        return app

    suffix = suffix_map[app.input_type]

    for row in app.peptide_files:
        file_path = row["path"]
        file_name = file_path.rsplit("/", 1)[-1].replace(suffix, "")
        peptide_file_names.append(file_name)
        peptide_groups.append(row["group"])

    app.peptide_file_names = peptide_file_names
    app.peptide_groups = peptide_groups
    return app


# --------------------------------
# PROTEIN FILE NAME EXTRACTION
# --------------------------------
def get_protein_file_names(app):
    suffix_map = {
        "Manually Created": "_proteins.csv",
        "MaxQuant": "_proteinGroups.txt",
        "Proteome Discoverer": "_Proteins.xlsx",
        "MSFragger": "_protein.tsv"
    }

    suffix = suffix_map.get(app.input_type, "_proteins.csv")

    protein_file_names = [
        pf.rsplit("/", 1)[-1].replace(suffix, "")
        for pf in app.protein_files
    ]

    app.protein_file_names = protein_file_names
    return app


# ----------------------------
# BASE RUN NAME EXTRACTION
# ----------------------------
def get_base_run_file_names(app):
    if not app.base_run_file_path:
        # default to first peptide file
        app.base_run_name = app.peptide_file_names[0]
        return app

    suffix_map = {
        "Manually Created": "_peptides.csv",
        "MaxQuant": "_peptides.txt",
        "Proteome Discoverer": "_PeptideGroups.xlsx",
        "MSFragger": "_peptides.tsv"
    }

    suffix = suffix_map.get(app.input_type, "")
    app.base_run_name = app.base_run_file_path.rsplit("/", 1)[-1].replace(suffix, "")

    return app


# ----------------------------
# SORT BASE RUN FIRST
# ----------------------------

def old_sort_base_run_first(app):

    if not app.base_run_name:
        return app

    base_name = app.base_run_name

    def custom_sort_key(item):
        return (0 if item == base_name else 1, item)

    combined = list(zip(app.peptide_file_names, app.peptide_files, app.peptide_groups))

    sorted_combined = sorted(combined, key=lambda x: custom_sort_key(x[0]))

    # unzip
    app.peptide_file_names, app.peptide_files, app.peptide_groups = map(list, zip(*sorted_combined))

    return app


def sort_base_run_first(app):
    if app.base_run_name == '':
        return app

    combined = zip(app.peptide_file_names, app.peptide_files, app.peptide_groups)

    base_name = app.base_run_name

    sorted_combined = sorted(combined, key=lambda x: (
            x[0] != base_name,        # False (0) for base run → comes first
            x[2] != app.peptide_groups[app.peptide_file_names.index(base_name)],  # base group first
            x[2] )                     # then alphabetical by group
    )
    (app.peptide_file_names, app.peptide_files, app.peptide_groups) = map(list, zip(*sorted_combined))
    return app


######################################################################

def import_data(app):
    #### Import and Initial Processing of Fasta File(s)
    fasta_df = fasta_import(app)
    print("Finished importing fasta_df", fasta_df.shape)
    print(fasta_df.head())
    print("point CCC")
    return fasta_df


    #### Import and Initial Processing of Peptide File(s)
    #peptide_input_df = peptide_import_function(app)
    #print("Finished importing Peptide Files")

    #### Import and Initial Processing of Protein File(s)
    #protein_input_df = protein_import_function(app)
    # print("Finished importing Protein Files")

    #return fasta_df #, peptide_input_df, protein_input_df







# ----------------------------
# RESULTS FOLDER CREATION
# ----------------------------
def create_results_folder(app):
    print("point B")
    results_path = "./results"
    os.makedirs(results_path, exist_ok=True)

    #full_path_with_project = os.path.join(results_path, app.results_file_name)
    full_path_with_project = results_path + "/" + app.results_file_name
    print(f"Full results file path: {full_path_with_project}")

    os.makedirs(full_path_with_project, exist_ok=True)
    print(f"Folder '{full_path_with_project}' created or already exists.")
    return full_path_with_project

# ----------------------------
# For Testing Purposes
# ----------------------------
def testing(app):
    # Insert preset testing parameters if in testing mode
    if app.testing_multiple == 1:
        app.analysis_type='multiple_proteins' 
    else:
        app.analysis_type='single_protein'
    app.input_type='Proteome Discoverer' 
    app.processing_method='rank_input' 
    app.calculation_type= 'log_2_intensity'  # 'PSM' or 'log_2_PSM' or 'intensity' or 'log_2_intensity' 
    app.N_value='55' 
    app.rank='77' 
    app.accession='P4242' 
    app.results_file_name='test_001_folder' 
    app.base_run_file_path = 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_3_A_1_PeptideGroups.xlsx'
    app.fasta_files=['C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/Human.fasta']
    """
    app.protein_files  =['C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_1_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_2_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_3_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_2_A_1_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_2_A_2_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_2_A_3_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_3_A_1_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_3_A_2_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_3_A_3_Proteins.xlsx']
    app.peptide_files= [{'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_1_PeptideGroups.xlsx', 'group': '3'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_2_PeptideGroups.xlsx', 'group': '3'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_3_PeptideGroups.xlsx', 'group': '3'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_2_A_1_PeptideGroups.xlsx', 'group': '2'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_2_A_2_PeptideGroups.xlsx', 'group': '2'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_2_A_3_PeptideGroups.xlsx', 'group': '2'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_3_A_1_PeptideGroups.xlsx', 'group': '5'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_3_A_2_PeptideGroups.xlsx', 'group': '5'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_3_A_3_PeptideGroups.xlsx', 'group': '5'}]  
    """
    # Short version for loading less data during testing
    app.protein_files  =['C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_1_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_2_Proteins.xlsx', 
                         'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_3_Proteins.xlsx']
    app.peptide_files= [{'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_1_PeptideGroups.xlsx', 'group': '3'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_2_PeptideGroups.xlsx', 'group': '3'}, 
                        {'path': 'C:/Users/dsmal/DMS/Projects/M03_Pro_Vet/data/NIST_data/1mg_1_A_3_PeptideGroups.xlsx', 'group': '3'}]  
    
    app.MW_var=1 
    app.peptide_length_var=1
    app.similarity_var=1 
    app.AA_composition_var=1 
    app.AA_type_var = 1 
    app.missed_cleavages_var=1 
    app.hydrophobicity_var=1 
    app.instability_var=1 
    app.pI_var=1
    app.helix_propensity_var=1 
    # app.pI_calc='Lehninger' 
    app.pKa_dict='EMBOSS'  # options: 'Bjellqvist', 'Lehninger', 'Stryer', 'EMBOSS'
    app.helix_type='pi' # options: 'alpha', 'pi'
    app.AA_index_001_var = 1
    app.AA_index_002_var = 1
    app.AA_index_003_var = 1
    app.AA_index_004_var = 1
    app.AA_index_005_var = 1
    app.AA_index_006_var = 1
    app.AA_index_007_var = 1
    app.AA_index_008_var = 1
    app.AA_index_009_var = 1
    app.AA_index_010_var = 1
    app.AA_index_011_var = 1
    app.AA_index_012_var = 1
    app.AA_index_013_var = 1
    app.AA_index_014_var = 1
    app.AA_index_015_var = 1
    app.AA_index_016_var = 1
    app.AA_index_017_var = 1
    app.AA_index_018_var = 1
    app.AA_index_019_var = 1
    app.AA_index_020_var = 1
    app.AA_index_999_var=1 
    AA_index_number_var='FASG760104'


    app.results_path = "./results"
    os.makedirs(app.results_path, exist_ok=True)
    full_path_with_project = app.results_path + "/" + app.results_file_name
    print(f"Full results file path: {full_path_with_project}")

    os.makedirs(full_path_with_project, exist_ok=True)
    print("Point A")
    print(f"Folder '{full_path_with_project}' created or already exists.")
    return app





# ----------------------------
# INITIAL PROCESSING PIPELINE
# ----------------------------
def initial_processing(app):
    # Load preset testing parameters if testing mode is on
    if app.testing_single == 1 or app.testing_multiple == 1:
        print("Testing mode is ON")
        print("loading Pre-set Testing Parameters")
        app = testing(app)

    print()
    print(app)
    print()
    print("Creating Results Path")
    app.results_subfolder_path = create_results_folder(app)
    print(f"app.results_subfolder path is: {app.results_subfolder_path}")

    print("Initial processing of input variables")

    app.N_value = int(app.N_value)
    app.rank = int(app.rank)

    # Extract names
    app = getting_fasta_file_names(app)
    print("fasta file names:", app.fasta_file_names)

    app = get_peptide_file_names(app)
    print("peptide file names:", app.peptide_file_names)
    print("groups:", app.peptide_groups)

    app = get_protein_file_names(app)
    print("protein file names:", app.protein_file_names)

    app = get_base_run_file_names(app)
    print("base run name:", app.base_run_name)

    # Sorting
    app = sort_base_run_first(app)

    print("\nAfter Sorting Base Run First:")
    print(app.peptide_file_names)
    
    print(app.peptide_files)
    print(app.peptide_groups)

    #fasta_df, peptide_input_df, protein_input_df = import_data(app)
    try:
        fasta_df = fasta_import(app)
        print ("shape of fasta_df: ", fasta_df.shape)
        print(fasta_df.head())
    except Exception as e:
        messagebox.showerror("Import Error", f"Error importing fasta files: {e}")
        print("Error importing fasta files:", e)
    
    try:
        peptide_input_df = peptide_import_function(app)
        print ("shape of peptide_input_df: ", peptide_input_df.shape)
        print(peptide_input_df.head())
    except Exception as e:
        messagebox.showerror("Import Error", f"Error importing peptide files: {e}")
        print("Error importing peptide files:", e)

    try:
        protein_input_df = protein_import_function(app)
        print ("shape of protein_input_df: ", protein_input_df.shape) 
        print(protein_input_df.head())
    except Exception as e:
        messagebox.showerror("Import Error", f"Error importing protein files: {e}")
        print("Error importing protein files:", e)

    # Protein_input_df processing for PSM values.  Pivot and sort decreasing order of abundance based on PSMs.
    protein_PSM_df = protein_PSM_pivot_and_sort(protein_input_df)
    print("Finished Pivoting protein_input_df to create protein_PSM_df")
    print(protein_PSM_df)

    # Create group dictionary from peptide groups and peptide file names
    
    group_dict, number_of_groups, unique_group_numbers = create_group_dict(app.peptide_groups, app.peptide_file_names)
    print("Created group dictionary from peptide groups and peptide file names")
    print(group_dict)
    print(number_of_groups)
    print(unique_group_numbers)

    # Create a deep copy of the group_dict
    deep_copy_group_dict = copy.deepcopy(group_dict)
    #deep_copied_dict = copy.deepcopy(original_dict)
    print(deep_copy_group_dict) 

    # Calculating Summary Information for Each Group
    print("Calculating summary information for each group...")
    summary_info_df = calc_group_info_basic(app.peptide_file_names, protein_input_df, peptide_input_df)
    print(summary_info_df)   
    save_dataframe_with_counter(summary_info_df, 
                                app.results_subfolder_path + "/" + app.results_file_name + "_summary_info", 
                                file_format='csv')
    print("Saved summary_info_df to file.")

    # Check to make sure there is at least two samples per group required for stats
    groups_with_1_entry = [k for k, v in group_dict.items() if len(v) < 2]

    if groups_with_1_entry:
        print("Groups with fewer than 2 values:", groups_with_1_entry)
        messagebox.showwarning("Insufficient Data Warning", 
                               f"The following groups have fewer than 2 values, which may affect statistical calculations: {groups_with_1_entry}")
        app.all_stats = False
    else:
        print("All groups have at least 2 values.")
        app.all_stats = True

    # purposefully error
    summary_stats_dict = summarize_by_group(summary_info_df, group_dict)        
    summary_stats_df = (pd.concat(summary_stats_dict, names=["property", "group"]).reset_index())
    print("Point D:", summary_stats_dict)

    summary_stats_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_summary_stats"
    save_dataframe_with_counter(summary_stats_df, summary_stats_path_name, 
                                file_format='csv')
    print(summary_stats_df)


    protein_and_peptide_figure = plots_for_protein_and_peptide_info(df = summary_info_df , dict_1 = summary_stats_dict)
    save_figure_unique(protein_and_peptide_figure, app.results_file_name + "_protein_and_peptide_figure", 
                       directory=app.results_subfolder_path, format="png")
    
    PSM_figure = plots_for_PSM_info(df = summary_info_df , dict_1 = summary_stats_dict)
    save_figure_unique(PSM_figure, app.results_file_name + "_PSM_figure", 
                       directory=app.results_subfolder_path, format="png")
    
    
    print("protein_and_peptide_figure and PSM_figure has been saved.")
    
    print(f"app.processing_method is: {app.processing_method}")

    if app.analysis_type == "multiple_proteins":
        print("Selected the multiple protein processing pipeline")
        top_N_processing_pipeline(app, fasta_df, peptide_input_df, 
                protein_input_df, group_dict, protein_PSM_df,
                summary_info_df = summary_info_df, summary_stats_df = summary_stats_df, 
                protein_and_peptide_figure=protein_and_peptide_figure, PSM_figure=PSM_figure)

    else:
        print("The Single Protein Processing Pipeline has been selected")
        single_protein_processing_pipeline(app, fasta_df, peptide_input_df, 
                protein_input_df, group_dict, protein_PSM_df,
                summary_info_df = summary_info_df, summary_stats_df = summary_stats_df, 
                protein_and_peptide_figure = protein_and_peptide_figure, PSM_figure = PSM_figure)        

    print("Processing complete.")

    return app

  ######################################################################





























"""
import os

# Get Fasta file names alone (# remove path from fasta file names)
def getting_fasta_file_names(app):
    fasta_file_names = []
    for long_fasta_name in app.fasta_files:
        fasta_file_name_i = long_fasta_name.rsplit("/", 1)[-1]
        fasta_file_names.append(fasta_file_name_i)
    app.fasta_file_names = fasta_file_names  # remove path from fasta file names
    return(app)


# Get peptide/protein file names from peptide/protein file names and paths
def old_get_peptide_file_names(app):
    peptide_file_names = []
    peptide_groups = []
    if app.input_type == "Manually Created":
        for i in range(0,len(app.peptide_files)):
            row_i = app.peptide_files[i]
            file_path_i = row_i['path']
            file_name_i = file_path_i.rsplit("/", 1)[-1].replace("_peptides.csv", "")
            peptide_file_names.append(file_name_i)
            group_i = row_i['group']
            peptide_groups.append(group_i)
        app.peptide_file_names = peptide_file_names
        app.peptide_groups = peptide_groups

        return(app)
    elif app.input_type == "MaxQuant":
        for i in range(0,len(app.peptide_files)):
            row_i = app.peptide_files[i]
            file_path_i = row_i['path']
            file_name_i = file_path_i.rsplit("/", 1)[-1].replace("_peptides.txt", "")
            peptide_file_names.append(file_name_i)
            group_i = row_i['group']
            peptide_groups.append(group_i)
            app.peptide_file_names = peptide_file_names
            app.peptide_groups = peptide_groups
            return(app)
        
    elif app.input_type == "Proteome Discoverer":
        for i in range(0,len(app.peptide_files)):
            row_i = app.peptide_files[i]
            file_path_i = row_i['path']
            file_name_i = file_path_i.rsplit("/", 1)[-1].replace("_PeptideGroups.xlsx", "")
            peptide_file_names.append(file_name_i)
            group_i = row_i['group']
            peptide_groups.append(group_i)
        app.peptide_file_names = peptide_file_names
        app.peptide_groups = peptide_groups
        return(app)
    else:
        print("Error: Unknown input type for peptide files.")
        return(app)
    


def old_get_protein_file_names(app):
    protein_file_names = []
    if app.input_type == "Manually Created":
        for long_protein_name in app.protein_files:
            protein_file_name_i = long_protein_name.rsplit("/", 1)[-1].replace("_proteins.csv", "")
            protein_file_names.append(protein_file_name_i)
        app.protein_file_names = protein_file_names
        return(app)
    
    elif app.input_type == "MaxQuant":
        for long_protein_name in app.protein_files:
            protein_file_name_i = long_protein_name.rsplit("/", 1)[-1].replace("_proteinGroups.txt", "")
            protein_file_names.append(protein_file_name_i)
        app.protein_file_names = protein_file_names
        return(app)
    
    elif app.input_type == "Proteome Discoverer":
        for long_protein_name in app.protein_files:
            protein_file_name_i = long_protein_name.rsplit("/", 1)[-1].replace("_Proteins.xlsx", "")
            protein_file_names.append(protein_file_name_i)
        app.protein_file_names = protein_file_names
        return(app)

def old_get_base_run_file_names(app):
    if app.base_run_file_path == "":
        app.base_run_file_name = app.peptide_file_names[0]
        #app.base_run_group = app.peptide_groups[0]
        return(app)
    
    if app.input_type == "Manually Created":
        app.base_run_name = app.base_run_file_path.rsplit("/", 1)[-1].replace("_peptides.csv", "")
        return(app)
    
    elif app.input_type == "MaxQuant":
        app.base_run_name = app.base_run_file_path.rsplit("/", 1)[-1].replace("_peptides.txt", "")
        return(app)

    elif app.input_type == "Proteome Discoverer":
        app.base_run_name = app.base_run_file_path.rsplit("/", 1)[-1].replace("_PeptideGroups.xlsx", "")
        return(app)


# Make base run is the first element in peptide_file_names, peptide_groups, and peptide_files
# create a small function to help sort
def custom_sort_key(item):
    if item == app.base_run_name:
        return (-1, item) # Assigns a primary sort key of -1 for 'aaa'
    else:
        return (0, item)  # Assigns a primary sort key of 0 for all others


def XXX_sort_base_run_first(app):
    if app.base_run_name == '':
        return app

    combined = zip(app.peptide_file_names, app.peptide_files, app.peptide_groups)

    base_name = app.base_run_name

    sorted_combined = sorted(combined, key=lambda x: (
            x[0] != base_name,        # False (0) for base run → comes first
            x[2] != app.peptide_groups[app.peptide_file_names.index(base_name)],  # base group first
            x[2] )                     # then alphabetical by group
    )
    (app.peptide_file_names, app.peptide_files, app.peptide_groups) = map(list, zip(*sorted_combined))
    return app



def old_initial_processing(app):
    create_results_folder(app)
    print("Point G")
   
    print("Initial processing of input variables")

    # convert top_N and rank to integers
    app.N_value = int(app.N_value)
    app.rank = int(app.rank)
    print("Point H")

    # Getting file names from file_paths
    app = getting_fasta_file_names(app)
    print("fasta file names are: ", app.fasta_file_names)
    app = get_peptide_file_names(app)
    print("peptide file names are: ", app.peptide_file_names)
    print("groups are: ", app.peptide_groups)
    app = get_protein_file_names(app)
    print("protein file names are: ", app.protein_file_names) 
    app = get_base_run_file_names(app)
    print("base run file names are: ", app.base_run_name)

    
    # Sort 
    app = sort_base_run_first(app)
    print()
    print("After Sorting Base Run First:")
    print(app.peptide_file_names)
    print(app.peptide_files)
    print(app.peptide_groups)

    print("Processing complete.")   


def create_results_folder(app):
    # Create a results folder and subfolder if they don't exist
    results_path = "./results"  # Define path to folder
    os.makedirs(results_path, exist_ok=True)
    results_folder_name = app.results_file_name
    full_path = results_path + "/" + results_folder_name
    print(f"Full results file path: {full_path}")
    os.makedirs(full_path, exist_ok=True)
    print(f"Folder '{full_path}' created or already exists.")
          
"""