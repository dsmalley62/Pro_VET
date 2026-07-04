import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import sys

from AA_index_dict import *

from results_viewer import ResultsViewer

#from utils import calculation_alternatives, create_similarity_matrix, getting_row_means_and_sds, \
#    make_coverage_plots, save_dataframe_with_counter, coverage_df_by_AA, save_figure_unique, \
#        create_similarity_matrix, plot_similarity_table, summarize_by_group, plot_runs_and_stats

from utils import *


# Make blank plot for unrequested plots
def old_make_blank_plot():
    fig, ax = plt.subplots(figsize=(4, 2))
    ax.text(0.5, 0.5, 'This plot was not requested!', horizontalalignment='center', 
            verticalalignment='center', fontsize=20, color='black', transform=ax.transAxes)
    ax.axis('off')
    return fig


# Create dataframe with only peptides from top-N proteins
def create_POIs_peptide_df(app, protein_PSM_df, peptide_input_df):
    POIs = []
    for i in range(int(app.N_value)):
        POI_i = protein_PSM_df.iloc[i] ['accession']
        POIs.append(POI_i)
        POIs_peptide_df = peptide_input_df[peptide_input_df['accession'].isin(POIs)].copy()
    return POIs_peptide_df

def create_AA_differences_df(POIs_PSM_df, app):
    print("POIs_PSM_df size: ", POIs_PSM_df.shape)
    # Get set of unique AAs
    unique_AAs = sorted(list(set("".join(POIs_PSM_df['seq'].astype(str)))))
    print(unique_AAs)
    overall_counts_per_AA = pd.DataFrame()
    for file_name in app.peptide_file_names:
        df_i = POIs_PSM_df[POIs_PSM_df['sample'] == file_name]
        print(file_name, df_i.shape)
        AA_df = pd.DataFrame()
        AA_PSM_df = pd.DataFrame()
        for AA in unique_AAs:
            AA_df[AA] = df_i['seq'].str.count(AA)
        counts_per_AA = pd.DataFrame(AA_df.sum())
        counts_per_AA.columns = [file_name]
        overall_counts_per_AA = pd.merge(overall_counts_per_AA, counts_per_AA, left_index=True, right_index=True, how='outer')
    return(overall_counts_per_AA)    



def calculate_overall_similarity(app, fasta_df, peptides_input_df, protein_PSM_df):
    # Drops rows with no accession.  (This is typically occurs for the reverse (incorrect) sequences)
    # (These should have already been removed but just to be safe.)
    peptides_input_df = peptides_input_df.dropna(subset=['accession'])
    
    all_value_means = []
    print("Please Wait.  This may take some time.")
    array_side_dim = len(app.peptide_files)
    similarity_array = np.zeros((int(app.N_value), array_side_dim, array_side_dim))
    for i in range(0, int(app.N_value)):
        POI_i = protein_PSM_df.iloc[(i)] ['accession']
        POI_fasta = fasta_df[fasta_df['accession'] == POI_i]
        POI_df = peptides_input_df[peptides_input_df['accession'].str.contains(POI_i)]
        #print("POint L:", POI_fasta['prot_seq'].item())
        #print(POI_fasta)
        POI_seq = POI_fasta['prot_seq'].item()
        calculations_df = calculation_alternatives(POI_df, app)
        coverage_df = coverage_df_by_AA(POI_seq, calculations_df)
        sim_matrix = create_similarity_matrix(coverage_df)
        similarity_array[i, :, :] = sim_matrix
        if i % 10 == 0:
            print(f"Processing: {i + 1} of {app.N_value}")
    average_similarity_array = np.sum(similarity_array, axis=0)/int(app.N_value)
    average_similarity_df = pd.DataFrame(average_similarity_array, columns = calculations_df.columns,
                                        index = calculations_df.columns)
    return(average_similarity_df)



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



def normalize_AA_differences(AA_df):
    column_sum = AA_df.sum(axis = 0)
    aver_col_sum = column_sum.mean()
    df_normalized_by_column = AA_df * aver_col_sum/column_sum
    df_normalized_by_column
    row_means = df_normalized_by_column.mean(axis=1)
    normalized_AA_df = df_normalized_by_column.div(row_means, axis=0)
    return(normalized_AA_df)    


def create_AA_type_df(AA_df):
    aliphatic_rows_to_sum = ['A', 'G', 'I', 'L', 'P', 'V']
    aromatic_rows_to_sum = ['F', 'W', 'Y']
    sulfur_containing_rows_to_sum = ['C', 'M']
    hydroxyl_rows_to_sum = ['S', 'T', 'Y']  # Y is both aromatic and hydroxyl
    basic_rows_to_sum = ['R', 'H', 'K']
    acidic_rows_to_sum = ['D', 'E']
    amide_rows_to_sum = ['N', 'Q']  

    #columns_to_sum = AA_df.columns
    # Calculate the column sums for the specified rows
    AA_type_df = pd.DataFrame()
    AA_type_df['aliphatic'] = AA_df.loc[aliphatic_rows_to_sum, : ].sum(axis=0)
    AA_type_df['aromatic'] = AA_df.loc[aromatic_rows_to_sum, : ].sum(axis=0)
    AA_type_df['sulfur'] = AA_df.loc[sulfur_containing_rows_to_sum, : ].sum(axis=0)
    AA_type_df['hydroxyl'] = AA_df.loc[hydroxyl_rows_to_sum, : ].sum(axis=0)
    AA_type_df['basic'] = AA_df.loc[basic_rows_to_sum, : ].sum(axis=0)
    AA_type_df['acidic'] = AA_df.loc[acidic_rows_to_sum, : ].sum(axis=0)
    AA_type_df['amide'] = AA_df.loc[amide_rows_to_sum, : ].sum(axis=0)
    # Transpose DataFrame
    AA_type_df = AA_type_df.T   
    return(AA_type_df)


def detecting_missed_cleavages(seq):
    seq_minus_N_termini = seq[:-1]
    count_K = seq_minus_N_termini.count('K')
    count_R = seq_minus_N_termini.count('R')
    count_KP = seq_minus_N_termini.count('KP')
    count_RP = seq_minus_N_termini.count('RP')
    missed_cleavages = (count_K + count_R) - (count_KP + count_RP)
    return(missed_cleavages)



def create_missed_cleavages_df(app, POIs_PSM_df):     
    print("Calculating Missed Cleavages Differences")
    # Examining Sequences and making dataframe indicating number of missed cleavages and PSM for each
    missing_cleavages_per_sample = []
    frequency_df = pd.DataFrame()
    for file_name in app.peptide_file_names:
        df_i = POIs_PSM_df[POIs_PSM_df['sample'] == file_name]
        output_i = []
        missing_i =  list(df_i['seq'].apply(detecting_missed_cleavages))
        for i in range(len(missing_i)):
            output_i.append(missing_i[i])
        frequency_i_df = pd.DataFrame(pd.Series(output_i).value_counts()).reset_index()
        #frequency_i_df = pd.DataFrame(output_i)
        frequency_i_df['sample'] = file_name
        frequency_df = pd.concat([frequency_df, frequency_i_df])
    frequency_df = frequency_df.rename(columns = {'index': 'missed cleavages'})
    return(frequency_df)




def create_missed_cleavages_plot(df_1, df_2, group_dict):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6)) 

    # -------- LEFT PLOT --------
    sns.barplot(
        x='missed cleavages',
        y='count',
        hue='sample',
        data=df_1,
        ax=ax1
    )

    ax1.set_xlabel('Missed Cleavages')
    ax1.set_ylabel("Counts (PSM × number)")
    ax1.set_title('Missed Cleavages per Run/Sample')

    # -------- RIGHT PLOT --------
    df_2 = df_2.reset_index()

    groups = list(group_dict.keys())
    n_groups = len(groups)

    missed = sorted(df_2["missed cleavages"].unique())
    x = np.arange(len(missed))

    width = 0.8 / n_groups  # auto-scale width

    for i, g in enumerate(groups):
        sub = (
            df_2[df_2["group"] == g]
            .sort_values("missed cleavages")
        )

        ax2.bar(
            x + i * width,
            sub["mean"],
            width=width,
            yerr=sub["sd"],
            capsize=5,
            label=f"Group {g}",
            edgecolor="black"
        )

    # Center ticks
    ax2.set_xticks(x + width * (n_groups - 1) / 2)
    ax2.set_xticklabels(missed)

    ax2.set_xlabel("Missed Cleavages")
    ax2.set_ylabel("Mean")
    ax2.set_title("Mean ± SD by Missed Cleavages and Group")
    ax2.legend()

    plt.tight_layout()
    return fig
 

def barplot_for_AA_figures(df: pd.DataFrame, norm_df: pd.DataFrame, save_path,
    graph_type = "AA_composition", dpi: int = 300,):
    """
    Create AA barplot figure.

    Parameters
    ----------
    df : pandas DataFrame
        Dataframe with AA counts.
    save_path : str or None
        File path to save the figure (e.g. 'AA_plot.png').
        If None, the figure is not saved.
    property : str or None
        options are 'AA composition', 'AA_type'
    dpi : int
        Resolution of saved figure.
    """
    # Adjust Dataframe
    df = df.reset_index()
    norm_df = norm_df.reset_index()
    df = df.rename(columns={'index': 'AA'})
    norm_df = norm_df.rename(columns={'index': 'AA'})
    df_long = pd.melt(df,
                  id_vars=['AA'],
                  var_name='sample',
                  value_name='value')
    norm_df_long = pd.melt(norm_df,
                  id_vars=['AA'],
                  var_name='sample',
                  value_name='value')
    palette = sns.color_palette("colorblind")

    # Create 2 subplots, stacked
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 10), sharex=True)
    # fig.suptitle('AA composition')

    # Plot the raw values of AA composition on top axes (axes[0])
    sns.barplot(x='AA', y='value', hue='sample', data=df_long, ax=axes[0], palette=palette)
    axes[0].tick_params(axis='x', labelbottom=True)
    if graph_type == "AA_type":
        axes[0].set_title('AA Types')
    else:
        axes[0].set_title('AA Composition of Detected Peptides')
    
    # axes[0].set_xlabel('Amino Acid')
    axes[0].set_ylabel('Relative Number of Each AA')
    axes[0].legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Plot the normalized values of AA composition on bottom axes (axes[1])
    sns.barplot(x='AA', y='value', hue='sample', data=norm_df_long, ax=axes[1], palette=palette)
    if graph_type == "AA_type":
        axes[1].set_title('AA Types (normalized)')
    else:
        axes[1].set_title('Normalized AA Composition')
    #axes[1].set_xlabel('Amino Acid')
    axes[1].set_ylabel('Relative AA Composition (normalized)')
    axes[1].legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    
    
    # ---- Save if requested ----
    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"Figure saved to {save_path}")

    return fig


def AA_stats_bar_chart(AA_stats_df, normalized_AA_stats_df, 
                       graph_type = "AA_composition", dpi: int = 300,):

    AA_stats_df = AA_stats_df.reset_index(drop=True)
    normalized_AA_stats_df = normalized_AA_stats_df.reset_index(drop=True)

    aas = AA_stats_df['AA'].unique()
    groups = sorted(AA_stats_df['group'].unique())

    n_aas = len(aas)
    n_groups = len(groups)

    x = np.arange(n_aas)

    # ---- Dynamic bar sizing ----
    max_cluster_width = 0.8
    bar_width = max_cluster_width / n_groups
    bar_width = max(bar_width, 0.05)

    # ---- Figure sizing ----
    fig_width = max(10, n_groups * 1.3)

    fig, axes = plt.subplots(
        2, 1,
        figsize=(fig_width, 9),
        sharex=True
    )

    #fig.suptitle('AA composition (per group)', fontsize=14, y=0.97)

    # ---- Plot helper ----
    def plot_panel(ax, df, title):
        for i, group in enumerate(groups):
            group_df = df[df['group'] == group]

            means = group_df['mean'].values
            sds = group_df['sd'].values

            offsets = (
                x
                - (bar_width * n_groups) / 2
                + i * bar_width
                + bar_width / 2
            )

            ax.bar(offsets, means, bar_width, yerr=sds, capsize=2, label=group, edgecolor='black',
                   linewidth=0.6, error_kw={'alpha': 0.7})
            

        ax.set_title(title, fontsize=11)
        ax.set_ylabel('Relative Number', fontsize=10)
        ax.tick_params(axis='y', labelsize=9)

        # ---- Legend for THIS axis ----
        ax.legend(
            title='Sample',
            loc='upper left',
            bbox_to_anchor=(1.01, 1.0),
            frameon=True,
            fontsize=9,
            title_fontsize=10,
            borderaxespad=0
        )

    # ---- Panels ----
    plot_panel(axes[0], AA_stats_df, 'AA Detection')
    plot_panel(axes[1], normalized_AA_stats_df, 'AA Detection (normalized)')

    # ---- X-axis formatting ----
    axes[0].tick_params(axis='x', labelbottom=True)
    
    # axes[1].set_xlabel('Amino Acid', fontsize=10)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(aas, rotation=45, ha='right', fontsize=9)

    # ---- Layout adjustments ----
    plt.subplots_adjust(
        left=0.08,
        right=0.82,   # space for legends
        top=0.92,
        bottom=0.15,
        hspace=0.25
    )

    return fig



def calc_MW_df(app, POIs_PSM_df, group_dict):
    # Residual mass from ionsource.com
    residual_MW	=pd.DataFrame({'G' : [57.021464], 'D' : [115.02694], 'A' : [71.037114], 'Q' : [128.05858], 
                           'S' : [87.032029], 'K': [128.09496], 'P' : [97.052764], 'E' : [129.04259], 
                           'V': [99.068414], 'M' : [131.04048], 'T' : [101.04768], 'H' : [137.05891], 
                           'C' : [160.0307], 'F' : [147.06841], 'L' : [113.08406], 'R' : [156.10111], 
                           'I' : [113.08406], 'Y' : [163.06333], 'N' : [114.04293],  'W' :  [186.07931]})
    MW_dict	= {'G' : 57.021464, 'D' : 115.02694, 'A' : 71.037114, 'Q' : 128.05858, 
                           'S' : 87.032029, 'K': 128.09496, 'P' : 97.052764, 'E' : 129.04259, 
                           'V': 99.068414, 'M' : 131.04048, 'T' : 101.04768, 'H' : 137.05891, 
                           'C' : 160.0307, 'F' : 147.06841, 'L' : 113.08406, 'R' : 156.10111, 
                           'I' : 113.08406, 'Y' : 163.06333, 'N' : 114.04293,  'W' :  186.07931} 

    # note C mass if for the carboxyamdiomethylated form, not the unmodified Cys.
    # The thiol group (–SH) is replaced by carbamidomethyl group (–CH2CONH2).
    # Delta mass is 57.02146
    # Cys unmodified residual mass is 103.00919 
    
    unique_AAs = list(set("".join(POIs_PSM_df['seq'].astype(str))))
    overall_MW_df = pd.DataFrame()
    for file_name in app.peptide_file_names:
        df_i = POIs_PSM_df[POIs_PSM_df['sample'] == file_name].reset_index()
        MW_df = pd.DataFrame()
        MW_2_list = []
        AA_df = pd.DataFrame()
        for AA in unique_AAs:
            # print(df_i['seq'])
            AA_df[AA] = df_i['seq'].str.count(AA) 
        for col_name, residual_weight in MW_dict.items():            
            if col_name in AA_df.columns:  # Check if the column exists in the DataFrame
                MW_df[col_name] = AA_df[col_name] * residual_weight
        MW_per_i = MW_df.sum(axis = 1) + 18.01056  # for N-terminal and C-terminal ends.
        df_MW_i = pd.DataFrame(MW_per_i)
        df_MW_i['sample'] = file_name 
        overall_MW_df = pd.concat([overall_MW_df, df_MW_i])
        
    overall_MW_df['MW'] = overall_MW_df[0]
    overall_MW_df = overall_MW_df[["MW", "sample"]] 
    
    # Add group column to dataframe
    # Build sample → group mapping
    sample_to_group = {sample: group
        for group, sample_list in group_dict.items()
        for sample in sample_list}

    overall_MW_df['group'] = overall_MW_df['sample'].map(sample_to_group)
    
    return(overall_MW_df.reset_index(drop = True))


def calc_gravy_index(app, POIs_PSM_df, group_dict):
    """
     Grand average of hydropathicity (GRAVY) of the linear polypeptide sequence is calculated as the sum of 
     hydropathy values of all amino acids, divided by the number of residues in the sequence. Increasing positive 
     score indicats greater hydrophobicity, but no account is taken of the way the protein folds in three 
     dimensions or the percentage of residues buried in the hydrophobic core of the protein. The calculation 
     is based on the Kyte-Doolittle scale.
     From: Kyte J, Doolittle RF (1982) A simple method for displaying the hydropathic character of a protein. J Mol Bio 157 (1), 105-32. 
    """
    hydrophobicity_dict = {'A' : 1.8, 'R' : -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q' : -3.5,'E': -3.5, 
                           'G' : -0.4, 'H' : -3.2, 'I': 4.5, 'L' : 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 
                           'P': -1.6, 'S' : -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    overall_gravy_df = pd.DataFrame( )
    unique_AAs = list(set("".join(POIs_PSM_df['seq'].astype(str))))
    for file_name in app.peptide_file_names:
        df_i = POIs_PSM_df[POIs_PSM_df['sample'] == file_name]
        AA_df = pd.DataFrame()
        gravy_2_list = []
        gravy_df = pd.DataFrame()
        for AA in unique_AAs:
            AA_df[AA] = df_i['seq'].str.count(AA)
            # print(AA_df)
        for col_name, multiplier_value in hydrophobicity_dict.items():
            #print(col_name, multiplier_value)
            if col_name in AA_df.columns:  # Check if the column exists in the DataFrame
                gravy_df[col_name] = AA_df[col_name] * multiplier_value
                gravy_i = list(gravy_df.sum(axis = 1)/AA_df.sum(axis = 1))
                #print("gravy_i", gravy_i)
        #gravy_2_list.append(gravy_i)
        df_gravy_i = pd.DataFrame(gravy_i)
        df_gravy_i['sample'] = file_name
        #print(df_gravy_i)

        overall_gravy_df = pd.concat([overall_gravy_df, df_gravy_i])
    overall_gravy_df = overall_gravy_df.rename(columns = {0:'Gravy Score'})       
    
    # Create a group column to dataframe
    sample_to_group = {sample: group
        for group, sample_list in group_dict.items()
        for sample in sample_list}

    overall_gravy_df['group'] = overall_gravy_df['sample'].map(sample_to_group)
    
    return(overall_gravy_df.reset_index(drop = True))



def calc_pI_per_seq(app, seq):
    # Possible dictionaries for calculating pIs.  All assume Cys is reduced and alklyated 
    # For Cys reduced and modified to have no net charge
    if app.pKa_dict == "EMBOSS":
        pKa_dict = {"Nterm": 8.0, "Cterm": 3.1, "D": 3.9, "E": 4.1, "H": 6.0, "K": 10.5, 
                    "R": 12.5, "Y": 10.1} #EMBOSS_dict
    elif app.pKa_dict == "Bjellqvist":
        pKa_dict = {"Nterm": 7.5, "Cterm": 3.55, "D": 4.05, "E": 4.45, "H": 5.98, "K": 10.0, 
                    "R": 12.0, "Y": 10.0} #Bjellqvist dictionary
    elif app.pKa_dict == "Lehninger":
        pKa_dict = {"Nterm": 9.6, "Cterm": 2.3, "D": 3.9, "E": 4.3, "H": 6.0, "K": 10.5, 
                    "R": 12.5, "Y": 10.1}
    elif app.pKa_dict == "Stryer": 
        pKa_dict = {"Nterm": 8.0, "Cterm": 3.1, "D": 3.9, "E": 4.1, "H": 6.5, "K": 10.8, 
                    "R": 12.5, "Y": 10.1} # Styer Dictionary
    else:
        print("Dictionary wasn't Selected Properly!!!")
        #sys.exit(1)

    if app.pI_method == "method 1":
         # Count residues
        counts = {AA: seq.count(AA) for AA in set(seq)}

        # Charge calculation at a given pH
        def net_charge(pH):
            charge = 0.0
            # N-terminus
            charge += 1 / (1 + 10**(pH - pKa_dict["Nterm"]))
            # C-terminus
            charge -= 1 / (1 + 10**(pKa_dict["Cterm"] - pH))
            # Side chains
            for aa, n in counts.items():
                if aa in ["K", "R", "H"]:  # basic
                    charge += n * (1 / (1 + 10**(pH - pKa_dict[aa])))
                elif aa in ["D", "E", "Y"]:  # acidic
                    charge -= n * (1 / (1 + 10**(pKa_dict[aa] - pH)))
            return charge
        # Binary search for pH where net charge = 0
        low, high = 0.0, 14.0
        for _ in range(50):  
            mid = (low + high) / 2
            if net_charge(mid) > 0:
                low = mid
            else:
                high = mid
        return round((low + high) / 2, 4)
    elif app.pI_method == "method 2":
        count_N_terminal = 1
        count_C_terminal = 1
        count_D = seq.count("D")
        count_E = seq.count("E")
        count_H = seq.count("H")
        count_K = seq.count("K")
        count_R = seq.count("R")
        count_Y = seq.count("Y")
        list_of_lists = [pKa_dict['Nterm']] * count_N_terminal, \
            [pKa_dict['Cterm']] * count_C_terminal, \
            [pKa_dict['D']] * count_D, \
            [pKa_dict['E']] * count_E, \
            [pKa_dict['H']] * count_H, \
            [pKa_dict['K']] * count_K, \
            [pKa_dict['R']] * count_R, \
            [pKa_dict['Y']] * count_Y
        single_list = [item for sublist in list_of_lists for item in sublist]
        sorted_list = sorted(single_list)
        len_sorted_list = len(sorted_list)
        pH_in_acid = count_N_terminal + count_K + count_R + count_Y + count_H
        # pKa_1 is where the net charge transitions from a positive value to zero.
        pKa_1 = sorted_list[pH_in_acid - 1 ]
        # pKa_2 is where the net charge transitions from zero to negative.
        pKa_2 = sorted_list[pH_in_acid ]
        pI = (pKa_1 + pKa_2)/2
        return(pI)
    else:
        print("method not selected")
        sys.exit(1)
 

def create_pI_df(app, POIs_PSM_df, group_dict):
    pI_calc_df = POIs_PSM_df
    # Apply a lambda function to the 'numbers' column to calculate pI for each value
    pI_calc_df['pI'] = pI_calc_df['seq'].apply(lambda x: calc_pI_per_seq(app, x))
    print(pI_calc_df)
    # Add group column
    # Create a group column to dataframe
    sample_to_group = {sample: group
        for group, sample_list in group_dict.items()
        for sample in sample_list}
    pI_calc_df['group'] = pI_calc_df['sample'].map(sample_to_group)
            
    return (pI_calc_df)


def calc_instability_index_per_seq(app, seq: str) -> float:
    # This function Generated by chatGPT and modified
    """
    Calculate the instability index of a peptide using
    Guruprasad et al. (1990) dipeptide weights.

    Parameters
    ----------
    app  # This is not used but keeping it for consistency
    seq : str
        Amino acid sequence (1-letter codes)

    Returns
    -------
    float
        Instability index
    """
    # Full 400 dipeptide instability weights  #### Check the numbers!!!!
    # Guruprasad et al. (1990) Protein Engineering 4(2):155–161

    instability_weights = {"AA":  1.0,  "AC": -1.5, "AD":  0.5, "AE": -0.5, "AF":  1.0, "AG":  0.0, "AH": -1.0, 
        "AI":  1.0, "AK": -1.0, "AL":  1.0, "AM":  1.0, "AN":  0.5, "AP":  0.0, "AQ": -0.5, "AR": -1.0, 
        "AS":  0.5, "AT":  0.5, "AV":  1.0, "AW":  1.0, "AY":  0.0,
        "CA": -1.0, "CC": -1.0, "CD": -1.0, "CE": -1.0, "CF":  1.0, "CG": -1.5, "CH": -1.0, "CI":  1.0, 
        "CK": -1.0, "CL":  1.0, "CM":  1.0, "CN": -1.0, "CP": -0.5, "CQ": -1.0, "CR": -1.0, "CS": -1.0, 
        "CT": -1.0, "CV":  1.0, "CW":  1.0, "CY": -1.0,
        "DA":  1.0, "DC": -1.5, "DD":  3.0, "DE":  2.0, "DF":  1.0, "DG":  0.5, "DH": -1.0, "DI":  1.0, 
        "DK":  0.5, "DL":  1.0, "DM":  1.0, "DN":  0.5, "DP":  0.0, "DQ":  0.0, "DR":  0.0, "DS":  0.5, 
        "DT":  0.5, "DV":  1.0, "DW":  1.0, "DY":  0.5,
        "EA": -1.0, "EC": -1.5, "ED":  3.0, "EE":  3.0, "EF":  1.0, "EG":  0.0,"EH": -1.0, "EI":  1.0, 
        "EK":  0.0, "EL":  1.0, "EM":  1.0, "EN":  0.5, "EP":  0.0, "EQ":  0.0, "ER":  0.0, "ES":  0.0, 
        "ET":  0.0, "EV":  1.0, "EW":  1.0, "EY":  0.0,
        "FA":  1.0, "FC":  1.0, "FD":  1.0, "FE":  1.0, "FF":  1.0, "FG":  1.0, "FH":  1.0, "FI":  1.0, 
        "FK":  1.0, "FL":  1.0, "FM":  1.0, "FN":  1.0, "FP":  1.0, "FQ":  1.0, "FR":  1.0, "FS":  1.0, 
        "FT":  1.0, "FV":  1.0, "FW":  1.0, "FY":  1.0,
        "GA":  1.0, "GC": -1.5, "GD":  2.0, "GE":  1.0, "GF":  1.0, "GG":  0.0, "GH": -1.0, "GI":  1.0, 
        "GK":  0.0, "GL":  1.0, "GM":  1.0, "GN":  0.5, "GP": -0.5, "GQ":  0.0, "GR": -1.0, "GS":  0.5, 
        "GT":  0.5, "GV":  1.0, "GW":  1.0, "GY":  0.0,
        "HA": -1.0, "HC": -1.0, "HD":  0.5, "HE": -0.5, "HF":  1.0, "HG": -0.5, "HH": -1.0, "HI":  1.0, 
        "HK": -1.0, "HL":  1.0, "HM":  1.0, "HN": -1.0, "HP": -0.5, "HQ": -0.5, "HR": -1.0, "HS": -1.0, 
        "HT": -1.0, "HV":  1.0, "HW":  1.0, "HY": -1.0,
        "IA":  1.0, "IC": -1.5, "ID":  1.0, "IE":  0.0, "IF":  1.0, "IG":  0.0, "IH": -1.0, "II":  1.0, 
        "IK":  0.0, "IL":  1.0, "IM":  1.0, "IN":  0.5, "IP":  0.0, "IQ":  0.0, "IR": -1.0, "IS":  0.0, 
        "IT":  0.0, "IV":  1.0, "IW":  1.0, "IY":  0.0,
        "KA": -1.0, "KC": -1.5, "KD":  0.0, "KE":  0.0, "KF":  0.5, "KG": -0.5, "KH": -1.0, "KI":  0.0, 
        "KK":  0.0, "KL":  0.0, "KM":  0.0, "KN": -0.5, "KP": -0.5, "KQ":  0.0, "KR":  0.0, "KS": -0.5, 
        "KT": -0.5, "KV":  0.0, "KW":  0.5, "KY": -0.5,
        "LA":  1.0, "LC": -1.5, "LD":  1.0, "LE":  1.0, "LF":  1.0, "LG":  0.0, "LH": -1.0, "LI":  1.0, 
        "LK":  0.0, "LL":  1.0, "LM":  1.0, "LN":  0.5, "LP":  0.0, "LQ":  0.0, "LR": -1.0, "LS":  0.0, 
        "LT":  0.0, "LV":  1.0, "LW":  1.0, "LY":  0.0,
        "MA":  1.0, "MC":  2.0, "MD":  1.0, "ME":  1.0, "MF":  1.0, "MG":  0.0, "MH": -1.0, "MI":  1.0, 
        "MK":  0.0, "ML":  1.0, "MM":  1.0, "MN":  0.5, "MP":  0.0, "MQ":  0.0, "MR": -1.0, "MS":  0.0, 
        "MT":  0.0, "MV":  1.0, "MW":  1.0, "MY":  0.0,
        "NA":  1.0, "NC": -1.5, "ND":  2.0, "NE":  1.0, "NF":  1.0, "NG":  0.5, "NH": -1.0, "NI":  1.0, 
        "NK":  0.0, "NL":  1.0, "NM":  1.0, "NN":  3.0, "NP":  0.0, "NQ":  0.0, "NR": -1.0, "NS":  0.5, 
        "NT":  0.5, "NV":  1.0, "NW":  1.0, "NY":  0.5,
        "PA":  1.0, "PC": -1.5, "PD":  0.0, "PE":  0.0, "PF":  1.0, "PG": -0.5, "PH": -1.0, "PI":  1.0, 
        "PK": -0.5, "PL":  1.0, "PM":  1.0, "PN":  0.0, "PP": -1.0, "PQ": -0.5, "PR": -1.0, "PS": -0.5, 
        "PT": -0.5, "PV":  1.0, "PW":  1.0, "PY":  0.0,
        "QA":  1.0, "QC": -1.5, "QD":  2.0, "QE":  1.0, "QF":  1.0, "QG":  0.5, "QH": -1.0, "QI":  1.0, 
        "QK":  0.0, "QL":  1.0, "QM":  1.0, "QN":  0.5, "QP":  0.0, "QQ":  0.5, "QR":  0.0, "QS":  0.0, 
        "QT":  0.0, "QV":  1.0, "QW":  1.0, "QY":  0.0,    
        "RA": -1.0, "RC": -1.5, "RD":  0.0, "RE":  0.0, "RF":  1.0, "RG": -0.5, "RH": -1.0, "RI":  1.0, 
        "RK":  0.0, "RL":  1.0, "RM":  1.0, "RN": -0.5, "RP": -0.5, "RQ":  0.0, "RR":  0.0, "RS": -0.5, 
        "RT": -0.5, "RV":  1.0, "RW":  1.0, "RY": -0.5,
        "SA": -1.0, "SC": -1.5, "SD":  0.5, "SE":  0.0, "SF":  1.0, "SG": -0.5, "SH": -1.0, "SI":  1.0, 
        "SK": -0.5, "SL":  1.0, "SM":  1.0, "SN":  0.0, "SP": -0.5, "SQ": -0.5, "SR": -1.0, "SS":  0.5, 
        "ST":  0.5, "SV":  1.0, "SW":  1.0, "SY":  0.0,
        "TA": -1.0, "TC": -1.5, "TD":  0.5, "TE":  0.0, "TF":  1.0, "TG": -0.5, "TH": -1.0, "TI":  1.0, 
        "TK": -0.5, "TL":  1.0, "TM":  1.0, "TN":  0.0, "TP": -0.5, "TQ": -0.5, "TR": -1.0, "TS":  0.5, 
        "TT":  0.5, "TV":  1.0, "TW":  1.0, "TY":  0.0,
        "VA":  1.0, "VC": -1.5, "VD":  1.0, "VE":  1.0, "VF":  1.0, "VG":  0.0, "VH": -1.0, "VI":  1.0, 
        "VK":  0.0, "VL":  1.0, "VM":  1.0, "VN":  0.5, "VP":  0.0, "VQ":  0.0, "VR": -1.0, "VS":  0.0, 
        "VT":  0.0, "VV":  1.0, "VW":  1.0, "VY":  0.0,
        "WA":  1.0, "WC": -1.5, "WD":  1.0, "WE":  1.0, "WF":  1.0, "WG":  0.0, "WH": -1.0, "WI":  1.0, 
        "WK":  0.0, "WL":  1.0, "WM":  1.0, "WN":  0.5, "WP":  0.0, "WQ":  0.0, "WR": -1.0, "WS":  0.0, 
        "WT":  0.0, "WV":  1.0, "WW":  1.0, "WY":  0.0,
        "YA":  1.0, "YC": -1.5, "YD":  1.0, "YE":  1.0, "YF":  1.0, "YG":  0.0, "YH": -1.0, "YI":  1.0, 
        "YK":  0.0, "YL":  1.0, "YM":  1.0, "YN":  0.5, "YP":  0.0, "YQ":  0.0, "YR": -1.0, "YS":  0.0, 
        "YT":  0.0, "YV":  1.0, "YW":  1.0, "YY":  0.5}
    
    seq = seq.upper()
    if len(seq) < 2:
        raise ValueError("Sequence must contain at least 2 residues.")

    score = 0.0
    for a, b in zip(seq, seq[1:]):
        dipep = a + b
        if dipep not in instability_weights:
            raise ValueError(f"Unknown dipeptide '{dipep}'")
        score += instability_weights[dipep]

    return (10.0 * score) / (len(seq) - 1)


def create_instability_df(app, POIs_PSM_df, group_dict):
    instability_calc_df = POIs_PSM_df
    # Apply a lambda function to the 'numbers' column to square each value
    instability_calc_df['Instability Index'] = \
        instability_calc_df['seq'].apply(lambda x: calc_instability_index_per_seq(app, x))
    sample_to_group = {sample: group
        for group, sample_list in group_dict.items()
        for sample in sample_list}
    instability_calc_df['group'] = instability_calc_df['sample'].map(sample_to_group)
            
    return (instability_calc_df)



def calc_helix_propensity_per_seq(app,seq):
    """
    Docstring for calc_helix_propensity_per_seq
    
    :param app: Description
    :param seq: Description
    :return: Description
    :rtype: float
    # Values for alpha-helix propensity index from -Pace CN, Scholtz JM. A helix propensity scale based on 
    #    experimental studies of peptides and proteins. Biophys J. 1998 Jul;75(1):422-7. 
    #    doi: 10.1016/s0006-3495(98)77529-0. PMID: 9649402; PMCID: PMC1299714.
    # M. N. Fodje & S. Al‑Karadaghi (2002). Occurrence, conformational features and amino acid propensities for the 
    #    π‑helix. Protein Engineering, Design & Selection 15(5):353–358.
    # values from pi-helix propensity index from Chat-GPT using method of Sen,et al.:
    #   Sen C, Logashree V, Makde RD, Ghosh B. Amino acid propensities for secondary structures and its variation 
    #    across protein structures using exhaustive PDB data. Comput Biol Chem. 2024 Jun;110:108083. 
    # doi: 10.1016/j.compbiolchem.2024.108083. Epub 2024 Apr 24. PMID: 38691894.

    """
    pi_helix_propensity = {'A': 0.85, 'C': 0.75, 'D': 0.60, 'E': 1.10, 'F': 1.20, 'G': 0.55,
                           'H': 0.90, 'I': 1.15, 'K': 1.05, 'L': 1.25, 'M': 1.00, 'N': 0.70, 'P': 0.30,
                           'Q': 0.95, 'R': 0.98, 'S': 0.80, 'T': 0.88, 'V': 1.10, 'W': 1.30, 'Y': 1.18}
    alpha_helix_propensity = {'A': 0.00, 'R': 0.21, 'N': 0.65, 'D': 0.69, 'C': 0.68, 'E': 0.40, 'Q': 0.39,
                              'G': 1.00, 'H': 0.61, 'I': 0.41, 'L': 0.21, 'K': 0.26, 'M': 0.24, 'F': 0.54, 
                              'P': 3.16, 'S': 0.50, 'T': 0.66, 'W': 0.49, 'Y': 0.53, 'V': 0.61}
    if app.helix_type == "pi":
        propensity_dict = pi_helix_propensity
    else:
        propensity_dict = alpha_helix_propensity
    len_seq = len(seq)
    value_x = 0
    for i in range(len_seq):
        value_x += propensity_dict[seq[i]]
    return(value_x/len_seq)


def create_helix_propensity_df(app, POIs_PSM_df, group_dict):
    propensity_calc_df = POIs_PSM_df
    # Apply a lambda function to the 'numbers' column to square each value
    propensity_calc_df['propensity'] = propensity_calc_df['seq'].apply(lambda x: calc_helix_propensity_per_seq(app, x))
    # Add group column
    # Create a group column to dataframe
    sample_to_group = {sample: group
        for group, sample_list in group_dict.items()
        for sample in sample_list}
    propensity_calc_df['group'] = propensity_calc_df['sample'].map(sample_to_group)
            
    return (propensity_calc_df)

    




#########################################################################

def top_N_processing_pipeline(app, fasta_df, peptide_input_df, 
                protein_input_df, group_dict, protein_PSM_df,
                summary_info_df, summary_stats_df, protein_and_peptide_figure, PSM_figure):
    
    # This allows you to look at all peptides by putting 0 as N value
    if app.processing_method == "top_N" and int(app.N_value) == 0:
        app.N_value = len(set(protein_input_df['accession']))
        print(f"app.N_value: {app.N_value}")


    # Create dataframe with only peptides from top-N proteins
    POIs_peptide_df = create_POIs_peptide_df(app, protein_PSM_df, peptide_input_df)
    print(f"Initial shape of protein_PSM_df: {protein_PSM_df.shape}")
    print(f"Initial shape of POIs_peptide_df: {POIs_peptide_df.shape}")
    print("POIs_peptide_df:")
    print(POIs_peptide_df)
    # Save the DataFrame
    POIs_peptide_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_POIs_peptide_df"
    save_dataframe_with_counter(POIs_peptide_df, POIs_peptide_df_path_name, file_format='csv')

    # New part
    # Repeat rows based on PSM
    POIs_PSM_df = POIs_peptide_df.loc[
        POIs_peptide_df.index.repeat(POIs_peptide_df['PSM'])
    ].reset_index(drop=True)

    # Drop columns
    POIs_PSM_df.drop(['PSM', 'intensity'], axis=1, inplace=True)

    print(POIs_PSM_df)  






    # Save the DataFrame
    POIs_PSM_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_POIs_PSM_df"
    save_dataframe_with_counter(POIs_PSM_df, POIs_PSM_df_path_name, file_format='csv')
    


    # Calculate overall similarity dataframe
    if app.similarity_var == 1:
        similarity_df = calculate_overall_similarity(app, fasta_df, peptide_input_df, protein_PSM_df)
        print("similarity dataframe has been created")
        print(similarity_df)

        # Save the DataFrame
        similarity_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_similarity_dataframe"
        save_dataframe_with_counter(similarity_df, similarity_df_path_name, file_format='csv')

        similarity_table = plot_similarity_table(similarity_df, app, 
                      column_widths=[0.25] * (similarity_df.shape[1] + 1))
    

        # Calculating stats for similarity dataframe
        similarity_row_for_calc = similarity_df.head(1)

        similarity_stats_dict = summarize_by_group(similarity_row_for_calc, group_dict,
                                               assume_equal_variance=True,
                                               exclude_columns= [similarity_row_for_calc.columns[0]])
        print("Similarity_stats_dict:", similarity_stats_dict)

        # convert the dictionary into a dataframe
        similarity_stats_df = pd.concat(similarity_stats_dict, names=["property", "group"]).reset_index()
        print("similarity stats df: ", similarity_stats_df)

        similarity_stats_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_similarity_stats"
        save_dataframe_with_counter(similarity_stats_df, similarity_stats_path_name, file_format='csv')

        # Create similarity figure   
        # only use first column for plotting
        similarity_df_1st_column = similarity_df.iloc[:, [0]]
        # Get the current name of the first column
        current_first_column_name = similarity_df_1st_column.columns[0]
        # Rename the column using a dictionary
        similarity_df_1st_column = similarity_df_1st_column.rename(columns={current_first_column_name: 'Similarity Index'}).reset_index()
        print(similarity_df_1st_column)
        similarity_plot = plot_runs_and_stats(df = similarity_df_1st_column, df_2 = similarity_stats_df, 
                                          property = 'Similarity Index')
        # Save the similarity_figure
        save_figure_unique(similarity_plot, app.results_file_name + "_similarity_figure", 
                       directory=app.results_subfolder_path, format="png")  
    
    else:
        similarity_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        similarity_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        similarity_table = pd
        similarity_plot = make_blank_plot()
    
    
    print("****************************************")

    if app.AA_composition_var == 1:
        print("Calculating AA Composition Differences for top-N proteins")
        AA_df = create_AA_differences_df(POIs_PSM_df, app)
        normalized_AA_df = normalize_AA_differences(AA_df)
        AA_df_for_figure = AA_df.T
        normalized_AA_df_for_figure = normalized_AA_df.T
        print("AA_df_for_figure:")
        print(AA_df_for_figure)
        print("normalized_AA_df_for_figure:")
        print(normalized_AA_df_for_figure)

    

        # Save the AA_df and normalized_AA_df
        AA_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_AA_df"
        normalized_AA_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_normalized_AA_df"
        save_dataframe_with_counter(AA_df, AA_df_path_name, file_format='csv')
        save_dataframe_with_counter(normalized_AA_df, normalized_AA_df_path_name, file_format='csv')

        app.results_subfolder_path = "./results/" + app.results_file_name

        # Create and save AA_figure
        AA_figure = barplot_for_AA_figures(AA_df, normalized_AA_df,
        save_path= app.results_subfolder_path + "/" + app.results_file_name + "_AA_plot.png",
        dpi=600)


        # Create AA_stats_dict and normalized_AA_stats_dict
        AA_stats_dict = summarize_by_group(AA_df, group_dict, assume_equal_variance=True)
        normalized_AA_stats_dict = summarize_by_group(normalized_AA_df, group_dict, assume_equal_variance=True)
        # Convert AA_stats_dict to dataframe
        AA_stats_df = (pd.concat(AA_stats_dict, names=["AA", "group"]).reset_index())
        normalized_AA_stats_df = (pd.concat(normalized_AA_stats_dict, names=["AA", "group"]).reset_index())
        print("AA stats df: ", AA_stats_df)
        print("normalized AA stats_df: ", normalized_AA_stats_df)

        # Save the AA_stats_df and normalized_AA_stats_df
        AA_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_AA_stats_df"
        normalized_AA_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_normalized_AA_stats_df"
        save_dataframe_with_counter(AA_stats_df, AA_stats_df_path_name, file_format='csv')
        save_dataframe_with_counter(normalized_AA_stats_df, normalized_AA_stats_df_path_name, file_format='csv')

        # Create AA_stats_figure
        print("working on Stat figure for AA compostion")
        AA_stats_figure = AA_stats_bar_chart(AA_stats_df, normalized_AA_stats_df)

        # Save the AA_stats_figure
        save_figure_unique(AA_stats_figure, app.results_file_name + "_AA_stats_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("AA_statsfigure has been saved.")

    else:
        AA_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        AA_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        AA_figure = make_blank_plot()
        AA_stats_figure = make_blank_plot()
        print("AA comparison not requested.")


    
    print("****************************************")

    if app.AA_type_var == 1:
        print("Calculating AA Type Differences for top-N proteins")
        try:
            # Attempt to access the variable
            # If the variable 'my_variable' is not defined, a NameError is raised
            print(f"The AA_df exists and its shape is: {AA_df.shape}")

        except NameError:
            # This code block runs only if the variable was not defined (NameError was caught)
            print("AA_df does not exist. Creating it")
            AA_df = create_AA_differences_df(POIs_PSM_df, app)
    
        AA_type_df = create_AA_type_df(AA_df)
        normalized_AA_type_df = normalize_AA_differences(AA_type_df)
        print("AA_type_df:")
        print(AA_type_df)

        print("Normalized AA type dataframe:")
        print(normalized_AA_type_df)
    
        # Save the AA_type_df and normalized_AA_type_df
        AA_type_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_AA_type_df"
        save_dataframe_with_counter(AA_type_df, AA_type_df_path_name, file_format='csv')
        normalized_AA_type_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_normalized_AA_type_df"
        save_dataframe_with_counter(normalized_AA_type_df, normalized_AA_type_df_path_name, file_format='csv') 

    
        # Calculate AA type stats and normalized AA type stats:
        AA_type_stats_dict = summarize_by_group(AA_type_df, group_dict, assume_equal_variance=True)
        normalized_AA_type_stats_dict = summarize_by_group(normalized_AA_type_df, group_dict, assume_equal_variance=True)    
    
        # convert the dictionaries into a dataframes
        AA_type_stats_df = (pd.concat(AA_type_stats_dict, names=["AA", "group"]).reset_index())
        normalized_AA_type_stats_df = (pd.concat(normalized_AA_type_stats_dict, names=["AA", "group"]).reset_index())
        print("AA type stats df: ", AA_type_stats_df)
        print("Normalized AA type stats df: ", normalized_AA_type_stats_df)

        # Save the AA_type_stats_df and the normalized_AA_type_stats_df
        AA_type_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_AA_type_stats_df"
        normalized_AA_type_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_normalized_AA_type_stats_df"
        save_dataframe_with_counter(normalized_AA_type_stats_df, normalized_AA_type_stats_df_path_name, file_format='csv')    
        
        # Create and save  AA_type_figure
        AA_type_figure = barplot_for_AA_figures(AA_type_df, normalized_AA_type_df, save_path= app.results_subfolder_path + 
                                               "/" + app.results_file_name + "_AA_type_plot.png", dpi=600)

        # Create AA_stats_dict and the normalized_AA_stats_dict
        AA_stats_dict = summarize_by_group(AA_df, group_dict, assume_equal_variance=True)
        normalized_AA_stats_dict = summarize_by_group(normalized_AA_df, group_dict, assume_equal_variance=True)
    
        # Convert AA_stats_dict and the normalized_AA_stats_dict to dataframes
        AA_stats_df = (pd.concat(AA_stats_dict, names=["AA", "group"]).reset_index(level=0))
        print("AA stats df: ", AA_stats_df)
        
        # Create AA_type_stats_figure
        AA_type_stats_figure = AA_stats_bar_chart(AA_type_stats_df, normalized_AA_type_stats_df)

        # Save the AA_type_stats_figure
        save_figure_unique(AA_type_stats_figure, app.results_file_name + "_AA_type_stats_figure", 
                    directory=app.results_subfolder_path, format="png")
        print("AA_statsfigure has been saved.")

    else:
        AA_type_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        AA_type_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        AA_type_figure = make_blank_plot()
        AA_type_stats_figure = make_blank_plot()
        print("AA Type calculation not requested.")



    print("****************************************")

    if (app.missed_cleavages_var == True):
        print("Calculating Missed Cleavages Differences for top-N proteins")
        long_missed_cleavages_df = create_missed_cleavages_df(app, POIs_PSM_df)
        # coverting to wide format
        missed_cleavages_df = long_missed_cleavages_df.pivot(index='missed cleavages', columns='sample', values='count')

        # convert the dictionary into a dataframe
        #missed_cleavages_stats_df = (pd.concat(missed_cleavages_stats_dict, names=["missed cleavages", "group"]).reset_index(level=0))
        #print("missed_cleavages_stats_df: ", missed_cleavages_stats_df)
        
    
        # Save the missed_cleavages_df
        missed_cleavages_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_missed_cleavages_df" 
        save_dataframe_with_counter(missed_cleavages_df, missed_cleavages_df_path_name, 
                                    file_format='csv') 


        #Stats for Missed Cleavages
        # Reshape missed_cleavages_df using pivot
        #missed_cleavages_df = long_missed_cleavages_df.pivot(index='missed cleavages', columns='sample', 
        #                                                     values='count')

        missed_cleavages_stats_dict = summarize_by_group(missed_cleavages_df, group_dict,
                                                         assume_equal_variance=True)
        # convert the dictionary into a dataframe
        missed_cleavages_stats_df = (pd.concat(missed_cleavages_stats_dict, names=["missed cleavages", "group"]).reset_index(level=0))
        print("missed_cleavages_stats_df: ", missed_cleavages_stats_df)
        # Save the missed_cleavages_stats_df
        missed_cleavages_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_missed_cleavages_stats_df"
        save_dataframe_with_counter(missed_cleavages_stats_df, missed_cleavages_stats_df_path_name, 
                                    file_format='csv') 

        missed_cleavages_figure = create_missed_cleavages_plot(long_missed_cleavages_df, 
                                                               missed_cleavages_stats_df , group_dict)  
        # Save the missed_cleavages_figure
        save_figure_unique(missed_cleavages_figure, app.results_file_name + "_missed_cleavages_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("missed_cleavages_figure has been saved.")
    else:
        missed_cleavages_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        missed_cleavages_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        missed_cleavages_figure = make_blank_plot()
        print("Missed Cleavages calculation not requested.")

    print("****************************************")
    



    # Create dataframe with peptide lengths

    if app.peptide_length_var == True:
        peptide_length_df = create_peptide_length_df(POIs_PSM_df, group_dict)

        # Save the peptide_length_df
        peptide_length_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_peptide_length_df" 
        save_dataframe_with_counter(peptide_length_df, peptide_length_df_path_name, 
                                file_format='csv') 

        # Calculate stats for peptide length
        peptide_length_stats_df = (peptide_length_df.groupby(["group", "pept_len"])["count"]
            .agg(mean="mean", std="std").reset_index())


        # Save the peptide_length_stats_df
        peptide_length_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_peptide_length_stats_df" 
        save_dataframe_with_counter(peptide_length_stats_df, peptide_length_stats_df_path_name, 
                                file_format='csv')
        
        # Plot peptide length by sample/run
        peptide_length_by_run_figure = plot_peptide_length_by_run(peptide_length_df, group_dict)
        

        # Save the peptide_length_by_run_figure
        save_figure_unique(peptide_length_by_run_figure, app.results_file_name + "_peptide_length_by_run_figure", 
                    directory=app.results_subfolder_path, format="png")
        print("peptide_length_by_run_figure has been saved.")

        # Plot peptide length by group
        peptide_length_by_group_figure = plot_peptide_length_by_group(peptide_length_stats_df, group_dict)

        # Save the peptide_length_by_group_figure
        save_figure_unique(peptide_length_by_group_figure, app.results_file_name + "_peptide_length_by_group_figure", 
                    directory=app.results_subfolder_path, format="png")
        print("peptide_length_by_group_figure has been saved.")

    else:
        peptide_length_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        peptide_length_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        peptide_length_by_run_figure = make_blank_plot()
        peptide_length_by_group_figure = make_blank_plot()
        print("Peptide Length calculation not requested.")

    print("****************************************")

    if app.MW_var == 1:
        print("Calculating Molecular Weight Differences for top-N proteins")
        MW_df = calc_MW_df(app, POIs_PSM_df, group_dict)
        print("MW_df:") 
        print(MW_df)    

        # Save the MW_df
        MW_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_MW_df" 
        save_dataframe_with_counter(MW_df, MW_df_path_name, file_format='csv')  

        # Calculate stats for MW
        MW_stats_df = calc_stats_by_group(MW_df, "MW")
        print("MW_stats_df: ", MW_stats_df)

        # Save the MW_stats_df
        MW_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_MW_stats_df" 
        save_dataframe_with_counter(MW_stats_df, MW_stats_df_path_name, 
                                file_format='csv')
        
        # Plot MW by sample/run
        MW_by_run_figure = plot_property_by_run(df = MW_df, 
            property = 'MW', group_dict=group_dict)
        
        # Save the MW_by_run_figure
        save_figure_unique(MW_by_run_figure, app.results_file_name + "_MW_by_run_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("MW_by_run_figure has been saved.")

        # Plot MW by group
        MW_by_group_figure = plot_property_by_group(df = MW_df, 
                    property = 'MW', group_dict=group_dict)

        # Save the MW_by_group_figure
        save_figure_unique(MW_by_group_figure, app.results_file_name + "_MW_by_group_figure", 
                    directory=app.results_subfolder_path, format="png")
        print("MW_by_group_figure has been saved.")
    else:
        MW_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        MW_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        MW_by_run_figure = make_blank_plot()
        MW_by_group_figure = make_blank_plot()
        print("Molecular Weight calculation not requested.")    
      

    
    print("****************************************")
    # Create dataframe with hydrophobicity calculations
    if app.hydrophobicity_var == 1:
        print("Calculating Hydrophobicity Differences for top-N proteins")
        hydrophobicity_df = calc_gravy_index(app, POIs_PSM_df, group_dict)
        print("hydrophobicity_df:") 
        print(hydrophobicity_df)    

        # Save the hydrophobicity_df
        hydrophobicity_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_hydrophobicity_df" 
        save_dataframe_with_counter(hydrophobicity_df, hydrophobicity_df_path_name, 
                                    file_format='csv')  

        # Calculate stats for hydrophobicity
        hydrophobicity_stats_df = calc_stats_by_group(hydrophobicity_df, "Gravy Score")
        hydrophobicity_stats_df

        # Save the hydrophobicity_stats_df
        hydrophobicity_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_hydrophobicity_stats_df" 
        save_dataframe_with_counter(hydrophobicity_stats_df, hydrophobicity_stats_df_path_name, 
                                    file_format='csv')
        
        # Plot hydrophobicity by sample/run
        hydrophobicity_by_run_figure = plot_property_by_run(df = hydrophobicity_df, 
                property = 'Gravy Score', group_dict=group_dict)
        
        # Save the hydrophobicity_by_run_figure
        save_figure_unique(hydrophobicity_by_run_figure, app.results_file_name + "_hydrophobicity_by_run_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("hydrophobicity_by_run_figure has been saved.")

        # Plot hydrophobicity by group
        hydrophobicity_by_group_figure = plot_property_by_group(df = hydrophobicity_df, 
                     property = 'Gravy Score', group_dict=group_dict)

        # Save the hydrophobicity_by_group_figure
        save_figure_unique(hydrophobicity_by_group_figure, app.results_file_name + "_hydrophobicity_by_group_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("hydrophobicity_by_group_figure has been saved.")
    else:
        hydrophobicity_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        hydrophobicity_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        hydrophobicity_by_run_figure = make_blank_plot()
        hydrophobicity_by_group_figure = make_blank_plot()
        print("Hydrophobicity calculation not requested.")
      

      

    print("****************************************")
    # Create dataframe with pI calculations

    app.pI_method = "method 2"

    if app.pI_var == 1:
        print("Calculating pI Differences for top-N proteins")
        pI_df = create_pI_df(app, POIs_PSM_df, group_dict)
        print("pI_df:") 
        print(pI_df)
        print("pI_df shape: ", pI_df.shape)
              
        # Save the pI_df
        pI_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_pI_df" 
        save_dataframe_with_counter(pI_df, pI_df_path_name, 
                                    file_format='csv')  

        
        # Calculate stats for pI
        pI_stats_df = calc_stats_by_group(pI_df, "pI")
        print(pI_stats_df)

        # Save the pI_stats_df
        pI_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_pI_stats_df" 
        save_dataframe_with_counter(pI_stats_df, pI_stats_df_path_name, 
                                    file_format='csv')



        # Plot pI by sample/run
        pI_by_run_figure = plot_property_by_run(df = pI_df, property = 'pI', 
                                                group_dict=group_dict)
        
        # Save the PI_by_run_figure
        save_figure_unique(pI_by_run_figure, app.results_file_name + "_pI_by_run_figure", 
                        directory=app.results_subfolder_path, format="png")
        print("pI_by_run_figure has been saved.")

        # Plot pI by group
        pI_by_group_figure = plot_property_by_group(df = pI_df, property = 'pI', group_dict=group_dict)
        
        # Save the PI_by_group_figure
        save_figure_unique(pI_by_group_figure, app.results_file_name + "_pI_by_group_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("pI_by_group_figure has been saved.")
    else:
        pI_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        pI_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        pI_by_run_figure = make_blank_plot()
        pI_by_group_figure = make_blank_plot()
        print("pI calculation not requested.")    

    print("****************************************")
    # Create dataframe with instability index calculations
    if app.instability_var == 1:
        print("Calculating Instability Index Differences for top-N proteins")
        instability_df = create_instability_df(app, POIs_PSM_df, group_dict)
        
        print(instability_df.shape)
        print("instability_df head", instability_df.head())

        # Save the instability_df
        instability_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_instability_df"   
        save_dataframe_with_counter(instability_df, instability_df_path_name, 
                                    file_format='csv')  

        # Calculate stats for Instability Index
        instability_stats_df = calc_stats_by_group(instability_df, "Instability Index")
        instability_stats_df

        # Save the instability_stats_df
        instability_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_instability_stats_df" 
        save_dataframe_with_counter(instability_stats_df, instability_stats_df_path_name, 
                                    file_format='csv')
        print("Instability Stats df:", instability_stats_df)
        
        # Plot Instability Index by sample/run
        instability_by_run_figure = plot_property_by_run(df = instability_df, property = 'Instability Index', 
                                                group_dict=group_dict)
        
        # Save the instability_by_run_figure
        save_figure_unique(instability_by_run_figure, app.results_file_name + "_instability_by_run_figure", 
                        directory=app.results_subfolder_path, format="png")
        print("Instability_by_run_figure has been saved.")

        # Plot Instability by group
        instability_by_group_figure = plot_property_by_group(df = instability_df, 
                    property = 'Instability Index', group_dict=group_dict)
        
        # Save the Instability_by_group_figure
        save_figure_unique(instability_by_group_figure, app.results_file_name + "_instability_by_group_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("instability_by_group_figure has been saved.")
    else:
        instability_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        instability_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        instability_by_run_figure = make_blank_plot()
        instability_by_group_figure = make_blank_plot()
        print("Instability Index calculation not requested.")

    print("****************************************")

    if app.helix_propensity_var == True:
        print(f"Calculating {app.helix_type}-helix propensity for top-N proteins")
        propensity_df = create_helix_propensity_df(app, POIs_PSM_df, group_dict)
        print(propensity_df.shape)
        print("propensity_df head", propensity_df.head())

        # Save the propensity_df
        propensity_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_helix_propensity_df"   
        save_dataframe_with_counter(propensity_df, propensity_df_path_name, 
                                    file_format='csv')  

        # Calculate stats for helix propensity Index
        propensity_stats_df = calc_stats_by_group(propensity_df, "propensity")
        propensity_stats_df

        # Save the propensity_stats_df
        propensity_stats_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + \
            "_helix_propensity_stats_df" 
        save_dataframe_with_counter(propensity_stats_df, propensity_stats_df_path_name, 
                                    file_format='csv')
        print("Helix Propensity Stats df:", propensity_stats_df)
        
        # Plot propensity to form helices by sample/run
        propensity_by_run_figure = plot_property_by_run(df = propensity_df, property = 'propensity', 
                                                group_dict=group_dict)
        
        # Save the propensity_by_run_figure
        save_figure_unique(propensity_by_run_figure, app.results_file_name + "_helix_propensity_by_run_figure", 
                        directory=app.results_subfolder_path, format="png")
        print("Helix_propensity_by_run_figure has been saved.")

        # Plot helix propensity by group
        propensity_by_group_figure = plot_property_by_group(df = propensity_df, 
                    property = 'propensity', group_dict=group_dict)
        
        # Save the helix propensity_by_group_figure
        save_figure_unique(propensity_by_group_figure, app.results_file_name + "_helix_propensity_by_group_figure", 
                       directory=app.results_subfolder_path, format="png")
        print("propensity_by_group_figure has been saved.")
    else:
        propensity_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        propensity_stats_df = pd.DataFrame({'0': ["This dataframe was not requested"]})
        propensity_by_run_figure = make_blank_plot()
        propensity_by_group_figure = make_blank_plot()
        print("Helix Propensity calculation not requested.")

    print("****************************************")

    AA_index_reference = 'ANDN920101'
    print("AA_index_001_var")

    AA_index_001_peptide_df, AA_index_001_stats_df, AA_index_001_plot_by_run, \
        AA_index_001_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_001_var", POIs_PSM_df, group_dict)
    
    print("****************************************")
    print("AA_index_002")

    AA_index_reference = 'ARGP820103'

    AA_index_002_peptide_df, AA_index_002_stats_df, AA_index_002_plot_by_run, \
        AA_index_002_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_002_var", POIs_PSM_df, group_dict)
    

    print("****************************************")
    print("AA_index_003")

    AA_index_reference = 'BIGC670101'

    AA_index_003_peptide_df, AA_index_003_stats_df, AA_index_003_plot_by_run, \
        AA_index_003_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_003_var", POIs_PSM_df, group_dict)
    
    print("****************************************")
    AA_index_reference = 'BROC820101'
    print("AA_index_004")

    AA_index_004_peptide_df, AA_index_004_stats_df, AA_index_004_plot_by_run, \
        AA_index_004_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_004_var", POIs_PSM_df, group_dict)

    print("****************************************")
    AA_index_reference = 'BULH740101'
    print("AA_index_005")

    AA_index_005_peptide_df, AA_index_005_stats_df, AA_index_005_plot_by_run, \
        AA_index_005_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_005_var", POIs_PSM_df, group_dict)

    print("****************************************")
    AA_index_reference = 'CHAM820102'
    print("AA_index_006")

    AA_index_006_peptide_df, AA_index_006_stats_df, AA_index_006_plot_by_run, \
        AA_index_006_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_006_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'HUTJ700103'
    print("AA_index_007")

    AA_index_007_peptide_df, AA_index_007_stats_df, AA_index_007_plot_by_run, \
        AA_index_007_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_007_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'JANJ780101'
    print("AA_index_008")

    AA_index_008_peptide_df, AA_index_008_stats_df, AA_index_008_plot_by_run, \
        AA_index_008_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_008_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'KANM800102'
    print("AA_index_009")

    AA_index_009_peptide_df, AA_index_009_stats_df, AA_index_009_plot_by_run, \
        AA_index_009_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_009_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'KANM800103'
    print("AA_index_010")

    AA_index_010_peptide_df, AA_index_010_stats_df, AA_index_010_plot_by_run, \
        AA_index_010_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_010_var", POIs_PSM_df, group_dict)
    

    print("****************************************")

    AA_index_reference = 'SUEM840101'
    print("AA_index_011_var")

    AA_index_011_peptide_df, AA_index_011_stats_df, AA_index_011_plot_by_run, \
        AA_index_011_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_011_var", POIs_PSM_df, group_dict)

    print("****************************************")
    print("AA_index_012")

    AA_index_reference = 'VELV850101'

    AA_index_012_peptide_df, AA_index_012_stats_df, AA_index_012_plot_by_run, \
        AA_index_012_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_012_var", POIs_PSM_df, group_dict)
    

    print("****************************************")
    print("AA_index_013")

    AA_index_reference = 'YUTK870101'

    AA_index_013_peptide_df, AA_index_013_stats_df, AA_index_013_plot_by_run, \
        AA_index_013_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_013_var", POIs_PSM_df, group_dict)

    print("****************************************")
    AA_index_reference = 'MITS020101'
    print("AA_index_014")

    AA_index_014_peptide_df, AA_index_014_stats_df, AA_index_014_plot_by_run, \
        AA_index_014_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_014_var", POIs_PSM_df, group_dict)
    
    print("****************************************")
    AA_index_reference = 'GUOD860101'
    print("AA_index_015")

    AA_index_015_peptide_df, AA_index_015_stats_df, AA_index_015_plot_by_run, \
        AA_index_015_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_015_var", POIs_PSM_df, group_dict)

    print("****************************************")
    AA_index_reference = 'ZHOH040103'
    print("AA_index_016")

    AA_index_016_peptide_df, AA_index_016_stats_df, AA_index_016_plot_by_run, \
        AA_index_016_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_016_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'OLSK800101'
    print("AA_index_017")

    AA_index_017_peptide_df, AA_index_017_stats_df, AA_index_017_plot_by_run, \
        AA_index_017_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_017_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'KARS160105'
    print("AA_index_018")

    AA_index_018_peptide_df, AA_index_018_stats_df, AA_index_018_plot_by_run, \
        AA_index_018_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_018_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'GUYH850102'
    print("AA_index_019")

    AA_index_019_peptide_df, AA_index_019_stats_df, AA_index_019_plot_by_run, \
        AA_index_019_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_019_var", POIs_PSM_df, group_dict)
    
    print("****************************************")

    AA_index_reference = 'ZIMJ680103'
    print("AA_index_020")

    AA_index_020_peptide_df, AA_index_020_stats_df, AA_index_020_plot_by_run, \
        AA_index_020_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_020_var", POIs_PSM_df, group_dict)
    
    print("****************************************")
    

    AA_index_reference = app.AA_index_number_var
    print("AA_index_user", AA_index_reference)

    try:
        AA_index_user_peptide_df, AA_index_user_stats_df, AA_index_user_plot_by_run, \
            AA_index_user_plot_by_group = AA_index_calc_function(app, AA_index_reference, 
                    "AA_index_999_var", POIs_PSM_df, group_dict)
        user_requested_description = aaindex1['CHOP780206'].description
        print("User requested AA index description: ", user_requested_description)
    except:
        print("User requested AA index not found.")
        user_requested_description = "Description not found."
        messagebox.showerror("Error", "AA Index reference not found in AAindex1 database.")

    
    
    print("Processing Finished for top-N proteins.")
    
    # Create Results Window
    
    results_window_multiprotein = ResultsViewer(summary_info_df, 
                                   summary_stats_df, 
                                   protein_and_peptide_figure,
                                   PSM_figure,
                                   similarity_df,
                                   similarity_stats_df,
                                   similarity_table,
                                   similarity_plot,
                                   AA_df,
                                   AA_figure,
                                   AA_stats_df,
                                   AA_stats_figure, 
                                   AA_type_df,
                                   AA_type_figure,
                                   AA_type_stats_df,
                                   AA_type_stats_figure,
                                   missed_cleavages_df,
                                   missed_cleavages_stats_df,
                                   missed_cleavages_figure,
                                   peptide_length_df,
                                   peptide_length_stats_df,
                                   peptide_length_by_run_figure,
                                   peptide_length_by_group_figure,
                                   MW_df,
                                   MW_stats_df,
                                   MW_by_run_figure,
                                   MW_by_group_figure,
                                   hydrophobicity_df,
                                   hydrophobicity_stats_df,
                                   hydrophobicity_by_run_figure,
                                   hydrophobicity_by_group_figure,
                                   pI_df,
                                    pI_stats_df,
                                    pI_by_run_figure,
                                    pI_by_group_figure,
                                    instability_df,
                                    instability_stats_df,
                                    instability_by_run_figure,
                                    instability_by_group_figure,
                                    propensity_df,
                                    propensity_stats_df,   
                                    propensity_by_run_figure,
                                    propensity_by_group_figure,
                                    AA_index_001_peptide_df, 
                                    AA_index_001_stats_df, 
                                    AA_index_001_plot_by_run, 
                                    AA_index_001_plot_by_group,
                                    AA_index_002_peptide_df, 
                                    AA_index_002_stats_df, 
                                    AA_index_002_plot_by_run, 
                                    AA_index_002_plot_by_group,
                                    AA_index_003_peptide_df, 
                                    AA_index_003_stats_df, 
                                    AA_index_003_plot_by_run, 
                                    AA_index_003_plot_by_group,
                                    AA_index_004_peptide_df, 
                                    AA_index_004_stats_df, 
                                    AA_index_004_plot_by_run, 
                                    AA_index_004_plot_by_group,
                                    AA_index_005_peptide_df, 
                                    AA_index_005_stats_df, 
                                    AA_index_005_plot_by_run, 
                                    AA_index_005_plot_by_group,
                                    AA_index_006_peptide_df, 
                                    AA_index_006_stats_df, 
                                    AA_index_006_plot_by_run, 
                                    AA_index_006_plot_by_group,
                                    AA_index_007_peptide_df, 
                                    AA_index_007_stats_df, 
                                    AA_index_007_plot_by_run, 
                                    AA_index_007_plot_by_group,
                                    AA_index_008_peptide_df, 
                                    AA_index_008_stats_df, 
                                    AA_index_008_plot_by_run, 
                                    AA_index_008_plot_by_group,
                                    AA_index_009_peptide_df, 
                                    AA_index_009_stats_df, 
                                    AA_index_009_plot_by_run, 
                                    AA_index_009_plot_by_group,
                                    AA_index_010_peptide_df, 
                                    AA_index_010_stats_df, 
                                    AA_index_010_plot_by_run, 
                                    AA_index_010_plot_by_group,
                                    AA_index_011_peptide_df, 
                                    AA_index_011_stats_df, 
                                    AA_index_011_plot_by_run, 
                                    AA_index_011_plot_by_group,
                                    AA_index_012_peptide_df, 
                                    AA_index_012_stats_df, 
                                    AA_index_012_plot_by_run, 
                                    AA_index_012_plot_by_group,
                                    AA_index_013_peptide_df, 
                                    AA_index_013_stats_df, 
                                    AA_index_013_plot_by_run, 
                                    AA_index_013_plot_by_group,
                                    AA_index_014_peptide_df, 
                                    AA_index_014_stats_df, 
                                    AA_index_014_plot_by_run, 
                                    AA_index_014_plot_by_group,
                                    AA_index_015_peptide_df, 
                                    AA_index_015_stats_df, 
                                    AA_index_015_plot_by_run, 
                                    AA_index_015_plot_by_group,
                                    AA_index_016_peptide_df, 
                                    AA_index_016_stats_df, 
                                    AA_index_016_plot_by_run, 
                                    AA_index_016_plot_by_group,
                                    AA_index_017_peptide_df, 
                                    AA_index_017_stats_df, 
                                    AA_index_017_plot_by_run, 
                                    AA_index_017_plot_by_group,
                                    AA_index_018_peptide_df, 
                                    AA_index_018_stats_df, 
                                    AA_index_018_plot_by_run, 
                                    AA_index_018_plot_by_group,
                                    AA_index_019_peptide_df, 
                                    AA_index_019_stats_df, 
                                    AA_index_019_plot_by_run, 
                                    AA_index_019_plot_by_group,
                                    AA_index_020_peptide_df, 
                                    AA_index_020_stats_df, 
                                    AA_index_020_plot_by_run, 
                                    AA_index_020_plot_by_group,
                                    AA_index_user_peptide_df, 
                                    AA_index_user_stats_df, 
                                    AA_index_user_plot_by_run, 
                                    AA_index_user_plot_by_group,
                                    user_requested_description,

                                    )
    
    results_window_multiprotein.run()

