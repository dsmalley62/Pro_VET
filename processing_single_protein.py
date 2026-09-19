import random
import re
from tkinter import messagebox
from utils import calculation_alternatives, create_similarity_matrix, getting_row_means_and_sds, make_coverage_plots, \
    save_dataframe_with_counter, coverage_df_by_AA, save_figure_unique, plot_runs_and_stats
from results_viewer import ResultsViewer_Single_Protein


def get_basic_info_about_POI(app, protein_PSM_df, fasta_df):
        POI_fasta = fasta_df[fasta_df['accession'] == app.accession]
        POI_seq = POI_fasta['prot_seq'].item()
        POI_description = protein_PSM_df.loc[protein_PSM_df['accession'] == app.accession, 'description'].item() 
        total_PSM_from_POI = int(protein_PSM_df.loc[protein_PSM_df['accession'] == app.accession, 'PSM'].item())
        POI_len = int(protein_PSM_df.loc[protein_PSM_df['accession'] == app.accession, 'prot_len'].item())
        POI_dict = {'accession': 'POI',
                'seq': POI_seq,
                'description': POI_description,
                'total_PSM': total_PSM_from_POI,
                'POI_len' : POI_len}
        return(POI_dict)


def QC_check_1(df):
    # Drops rows with no accession.  (This is typically occurs for the reverse (incorrect) sequences)
    print(f"peptide_input_df shape prior to removing Nas: {df.shape}")
    df = df.dropna(subset=['accession'])
    print(f"peptide_input_df shape after removing Nas in accession number: {df.shape}")
    print("Nas (no value) in accession number usually is due to incorrect (reverse or scrambled) sequence")

    # Count NaN values in each column
    nan_counts = df.isnull().sum()
    df_nas_in_A = df[df['accession'].isna()]
    return df


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






def single_protein_processing_pipeline(app, fasta_df, peptide_input_df, 
                                       protein_input_df, group_dict, protein_PSM_df, 
                                       summary_info_df, summary_stats_df, 
                                       protein_and_peptide_figure, PSM_figure,
                                       met_df, met_stats_df, Met_Ox_fig): 
    print("*********************************************************")
    print("Starting single protein processing pipeline...")
    print(f"Single Protein Processing Method is: {app.processing_method}")
    POI = processing_alternatives(protein_PSM_df, app)
    app.accession = POI
    print(f"POI is {app.accession}")
    
    POI_dict = get_basic_info_about_POI(app, protein_PSM_df, fasta_df)
    print(f" Print POI_dict is: {POI_dict}")

    # Quality Control Check for NaN values in peptide_input_df
    peptide_input_df = QC_check_1(peptide_input_df)
    

    ### Create POI_df
    # Drops rows with no accession.  (This is typically occurs for the reverse (incorrect) sequences)
    peptides_input_df = peptide_input_df.dropna(subset=['accession'])  
    # I know this is done above, but just to be sure

    POI_df = peptide_input_df[peptide_input_df['accession'].str.contains(app.accession, na=False)]
    Total_number_of_PMS_from_peptides_df = POI_df['PSM'].sum()
    print(f"PSMs based on peptides dataframe: {Total_number_of_PMS_from_peptides_df}")
    # Save the DataFrame
    POI_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_POI_df"
    save_dataframe_with_counter(POI_df, POI_df_path_name, 
                                    file_format='csv')
    print("POI_df has been saved in results folder")


    # create calculations_df
    print("creating Calculations dataframe...")
    calculations_df = calculation_alternatives(POI_df, app)

    # Save the DataFrame
    calculations_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_calculations_df"
    save_dataframe_with_counter(calculations_df, calculations_df_path_name, file_format='csv')

    print("calculations_df has been saved in results folder")
    print(calculations_df)

    # Create Coverage save_dataframe_with_counter
    coverage_df = coverage_df_by_AA(POI_dict['seq'], calculations_df)

    # Save the DataFrame
    coverage_df_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_coverage_df"
    save_dataframe_with_counter(coverage_df, coverage_df_path_name, file_format='csv')

    print("coverage_df has been saved in results folder")
    print(coverage_df)


    similarity_matrix = create_similarity_matrix(coverage_df)
    print("similarity matrix has been created")
    print(similarity_matrix)
    similarity_matrix_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_similarity_matrix"
    save_dataframe_with_counter(similarity_matrix, similarity_matrix_path_name, file_format='csv')

    # Calculating stats for similarity matrix
    similarity_calculations_row = similarity_matrix.iloc[0:1, :]
    print("similarity calculations row:", similarity_calculations_row)
    similarity_stats = getting_row_means_and_sds(similarity_calculations_row, 
                                               group_dict = group_dict, use_base_group = False)
    print("similarity stats:")
    print(similarity_stats)
    similarity_stats_path_name = app.results_subfolder_path + "/" + app.results_file_name + "_similarity_stats"
    save_dataframe_with_counter(similarity_stats, similarity_stats_path_name, file_format='csv')    



    # Create coverage plots
    print("Creating coverage plots...")
    similarity_figure = make_coverage_plots(coverage_df, app, POI_dict)

    # Save the figure, ensuring a unique filename
    save_figure_unique(similarity_figure, app.results_file_name + "_similarity_plot", directory =app.results_subfolder_path, format = 'pdf')

    # If you run this multiple times, it will create my_plot.png, my_plot_1.png, my_plot_2.png, etc.


    results_window = ResultsViewer_Single_Protein(app,
                                                  POI_dict,
                                                  summary_info_df, 
                                   summary_stats_df, 
                                   protein_and_peptide_figure = protein_and_peptide_figure,
                                   PSM_figure = PSM_figure,
                                   calculations_df=calculations_df, 
                                   POI_df=POI_df,
                                   met_df = met_df,
                                   met_stats_df = met_stats_df,
                                   Met_Ox_fig = Met_Ox_fig,
                                   similarity_matrix=similarity_matrix,
                                   similarity_stats=similarity_stats,
                                   coverage_df=coverage_df,
                                   similarity_figure=similarity_figure)

    results_window.run()



