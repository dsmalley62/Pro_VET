import tkinter as tk
from tkinter import ttk
from tkinter import font
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pandastable import Table
from tkinter import filedialog, ttk, messagebox


import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 150        # display quality
mpl.rcParams['savefig.dpi'] = 300       # saved quality
mpl.rcParams['figure.facecolor'] = 'white'



class Results_View:

    def __init__(self, results_window_root):
        self.root = results_window_root

        self.root.title("Results Viewer")
        self.root.geometry("1200x700")

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(expand=1, fill="both")

        # Tabs
        self.tab_results_1 = ttk.Frame(self.notebook)
        self.tab_results_2 = ttk.Frame(self.notebook)


        self.notebook.add(self.tab_results_1, text="Results (Page 1)")
        self.notebook.add(self.tab_results_2, text="Results (Page 2)")


        # Build UI sections
        self.build_tab_results_1()
        self.build_tab_results_2()

    ###############################################################
    # TAB 1 — Results Page 1
    ###############################################################

    def build_tab_results_1(self):
        top = tk.Frame(self.tab_results_1, bd=2, relief="groove", bg="darkblue", height=180)
        top.pack(fill="x", padx=5, pady=5)

        self.build_tab1_top(top)

        bottom = tk.Frame(self.tab_input, bd=2, relief="groove", bg="lightgreen")
        bottom.pack(fill="both", expand=True, padx=5, pady=5)

        self.build_file_inputs(bottom)



class ResultsViewer_Single_Protein:
    def __init__(
        self,
        app,
        POI_dict,
        summary_info_df,
        summary_stats_df,
        protein_and_peptide_figure,
        PSM_figure,
        calculations_df,
        POI_df,
        similarity_matrix,
        similarity_stats,
        coverage_df,
        similarity_figure
    ):
        """
        Parameters
        ----------
        app: object
        POI_dict : dict
        summary_info_df : pd.DataFrame
        summary_stats_df : pd.DataFrame
        protein_and_peptide_figure : matplotlib.figure.Figure  
        PSM_figure : matplotlib.figure.Figure
        calculations_df : pd.DataFrame
        POI_df : pd.DataFrame
        similarity_matrix : pd.DataFrame
        similarity_stats : pd.DataFrame
        coverage_df : pd.DataFrame
        similarity_figure : matplotlib.figure.Figure
        """
        self.bold_label_font = ("TkDefaultFont", 10, "bold")
        self.app = app
        self.summary_info_df = summary_info_df
        self.summary_stats_df = summary_stats_df
        self.protein_and_peptide_figure = protein_and_peptide_figure
        self.PSM_figure = PSM_figure
        self.calculations_df = calculations_df
        self.POI_df = POI_df
        self.similarity_matrix = similarity_matrix
        self.similarity_stats = similarity_stats
        self.coverage_df = coverage_df
        self.similarity_figure = similarity_figure

        # Create the main application window
        self.results = tk.Tk()
        self.results.title("Four Frames Example")
        self.results.geometry("1200x600") # Set an initial size for visibility
        self.results.configure(background="lightgrey")

        # Iterate to create four frames
        title_frame = tk.Frame(self.results, bg="lightgrey", height=10)
        title_frame.pack(fill=tk.X, expand=True, anchor=tk.N)
        
        label = tk.Label(title_frame, text="Results for Single Protein Analysis", bg="lightgrey", font=("Arial", 16))
        label.pack(pady=5)

        general_info_frame = tk.Frame(self.results, bg="lightgrey", height = 10)
        general_info_frame.pack(fill=tk.X, expand=True, anchor = tk.N)


        # Explicitly place labels 
        label_0_0 = tk.Label(general_info_frame, text="Method for Protein Selection:", borderwidth=0, relief="solid", 
                             anchor = "w", bg = "lightgrey",font=("Helvetica", 12, "bold"))
        label_0_0.grid(row=0, column=0, padx=5, pady=5, sticky="w")
        label_0_1 = tk.Label(general_info_frame, text= app.processing_method, borderwidth=0, bg = "lightgrey", 
                             relief="solid", anchor = "w", font=("Helvetica", 12))
        label_0_1.grid(row=0, column=1, padx=15, pady=5, sticky="w")
        
        label_1_0 = tk.Label(general_info_frame, text="Accession Number: ", borderwidth=0, bg = "lightgrey", 
                             relief="solid", anchor = "w", font=("Helvetica", 12, "bold"))
        label_1_0.grid(row=1, column=0, padx=5, pady=5, sticky="w")
        label_1_1 = tk.Label(general_info_frame, text= app.accession, borderwidth=0, relief="solid", anchor = "w", 
                             bg = "lightgrey", font=("Helvetica", 12))
        label_1_1.grid(row=1, column=1, padx=15, pady=5, sticky="w")
        
        protein_name = POI_dict['description'].split(" OS=")[0]
        label_2_0 = tk.Label(general_info_frame, text="Protein: ", borderwidth=0, bg = "lightgrey", relief="solid", 
                             anchor = "w", font=("Helvetica", 12, "bold"))
        label_2_0.grid(row=2, column=0, padx=5, pady=5, sticky="w")
        label_2_1 = tk.Label(general_info_frame, text=protein_name, borderwidth=0, bg = "lightgrey", relief="solid", 
                             anchor ="w", font=("Helvetica", 12))
        label_2_1.grid(row=2, column=1, padx=15, pady=5, sticky="w")
        
        label_3_0 = tk.Label(general_info_frame, text="Total number of PSMs:", borderwidth=0, bg = "lightgrey", 
                             relief="solid", anchor = "w",font=("Helvetica", 12, "bold"))
        label_3_0.grid(row=3, column=0, padx=5, pady=5, sticky="w")
        label_3_1 = tk.Label(general_info_frame, text=str(POI_dict['total_PSM']), borderwidth=0, bg = "lightgrey", 
                             relief="solid", anchor ="w",font=("Helvetica", 12))
        label_3_1.grid(row=3, column=1, padx=15, pady=5, sticky="w")

        label_4_0 = tk.Label(general_info_frame, text="Protein Length (AAs):", borderwidth=0, bg = "lightgrey", 
                             relief="solid", anchor = "w", font=("Helvetica", 12, "bold"))
        label_4_0.grid(row=4, column=0, padx=5, pady=5, sticky="w")
        label_4_1 = tk.Label(general_info_frame, text=str(POI_dict['POI_len']), borderwidth=0, bg = "lightgrey", 
                             relief="solid", anchor ="w", font=("Helvetica", 12))
        label_4_1.grid(row=4, column=1, padx=15, pady=5, sticky="w")
        
        
        files_frame = tk.Frame(self.results, bg="lightgrey", height=10)
        files_frame.pack(fill=tk.X, expand=True, anchor = tk.N)

        tk.Label(files_frame,
            text="Raw Data", width=25, height = 2, font=self.bold_label_font, bg = "lightgrey"
            ).grid(row=0, column=1, padx=4, pady=4) 
        
        tk.Label(files_frame,
            text="Means, Stds, and P-Values", width=25, height = 2, font=self.bold_label_font, bg = "lightgrey"
            ).grid(row=0, column=2, padx=4, pady=4)   
        
        tk.Label(files_frame,
            text="First Figure", width=25, height = 2, font=self.bold_label_font, bg = "lightgrey"
            ).grid(row=0, column=3, padx=4, pady=4)   
        
        tk.Label(files_frame,
            text="Second Figure", width=25, height = 2, font=self.bold_label_font, bg = "lightgrey"
            ).grid(row=0, column=4, padx=4, pady=4) 
        
        tk.Label(files_frame,
            text="Proteins/Peptides/PSMs", width=25, height = 2, font=self.bold_label_font, bg = "lightgrey"
            ).grid(row=1, column=0, padx=4, pady=4)  
        
        tk.Button(files_frame,
            text="Summary Information",
            command=lambda: open_dataframe_window(self.summary_info_df),
            width=25, height = 2
        ).grid(row=1, column=1, padx=4, pady=4) 
        
        tk.Button(files_frame,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.summary_stats_df),
            width=25, height = 2
        ).grid(row=1, column=2, padx=4, pady=4)   
        
        tk.Button(files_frame,
            text="Protein and Peptide Figure",
            command=lambda: show_figure_in_tk(self.results, self.protein_and_peptide_figure),
            width=25, height = 2
        ).grid(row=1, column=3, padx=4, pady=4) 

        tk.Button(files_frame,
            text="PSM Figure",
            command=lambda: show_figure_in_tk(self.results, self.PSM_figure),
            width=25, height = 2
        ).grid(row=1, column=4, padx=4, pady=4) 

        ##################################################################
        
        tk.Label(files_frame,
            text="Similarity Index", width=25, height = 2, font=self.bold_label_font, bg ="lightgrey"
            ).grid(row=2, column=0, padx=4, pady=4) 
        
        
        tk.Button(files_frame,
            text="Similarity Index Dataframe",
            command=lambda: open_dataframe_window(self.similarity_matrix),
            width=25, height = 2
        ).grid(row=2, column=1, padx=4, pady=4)   
        
        tk.Button(files_frame,
            text="Similarity Index Stats",
            command=lambda: open_dataframe_window(self.similarity_stats),
            width=25, height = 2
        ).grid(row=2, column=2, padx=4, pady=4)
         
        tk.Button(files_frame,
            text="Similarity Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.results, self.similarity_figure)
        ).grid(row=2, column=3, padx=4, pady=4 )
        
        ###############################################################
        tk.Label(files_frame,
            text="Coverage Information", width=25, height = 2, font=self.bold_label_font, bg ="lightgrey"
            ).grid(row=3, column=0, padx=4, pady=4) 
        
        
        tk.Button(files_frame,
            text="Coverage Dataframe",
            command=lambda: open_dataframe_window(self.coverage_df),
            width=25, height = 2
        ).grid(row=3, column=1, padx=4, pady=4)  

        tk.Label(files_frame,
            text="Peptide Information", width=25, height = 2, font=self.bold_label_font, bg="lightgrey"
            ).grid(row=4, column=0, padx=4, pady=4) 
        
        
        tk.Button(files_frame,
            text="Peptide Dataframe",
            command=lambda: open_dataframe_window(self.POI_df),
            width=25, height = 2
        ).grid(row=4, column=1, padx=4, pady=4)  

        ###############################################

        statement_frame = tk.Frame(self.results, bg="white", height=15)
        statement_frame.pack(fill=tk.X, expand=True, anchor = tk.N)

        statement = "All files and data should have been saved to results folder.  Image quality is better in the saved files."
        label_statement = tk.Label(statement_frame, text= statement, wraplength = 1000, borderwidth=0, relief="solid", 
                                   font = ("Helvetica", 14), anchor = "w", justify=tk.LEFT, bg = "lightgrey")
        label_statement.pack(padx=5, pady=5, fill=tk.X)

    def run(self):
        #Start the Tkinter event loop.
        self.root.mainloop()
   


        
        
  



def open_dataframe_window(dataframe):
        new_window = tk.Toplevel()
        new_window.title("Pandas DataFrame Viewer")

        frame = tk.Frame(new_window)
        frame.pack(fill=tk.BOTH, expand=True)

        table = Table(frame, dataframe=dataframe, showtoolbar=True, showstatusbar=True, showindex=True)
        table.show()




def show_figure_in_tk(parent, fig, window_title="Visualization of Data"):
    win = tk.Toplevel(parent)
    win.title(window_title)

    canvas = FigureCanvasTkAgg(fig, master=win)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)




class ResultsViewer:

    def __init__(self, 
                 summary_df, summary_stats_df, protein_and_peptide_figure, PSM_figure,
                 similarity_df, similarity_stats, similarity_table, similarity_plot,
                 AA_df,AA_figure, AA_stats_df, AA_stats_figure,
                 AA_type_df, AA_type_figure, AA_type_stats_df, AA_type_stats_figure,
                 missed_cleavages_df, missed_cleavages_stats_df, missed_cleavages_figure,
                 peptide_length_df, peptide_length_stats_df, peptide_length_by_run_figure,
                 peptide_length_by_group_figure,
                 MW_df, MW_stats_df, MW_by_run_figure, MW_by_group_figure,
                hydrophobicity_df, hydrophobicity_stats_df, hydrophobicity_by_run_figure,
                hydrophobicity_by_group_figure,
                pI_df, pI_stats_df, pI_by_run_figure, pI_by_group_figure,
                instability_df, instability_stats_df, instability_by_run_figure, instability_by_group_figure,
                propensity_df, propensity_stats_df, propensity_by_run_figure, propensity_by_group_figure,
                AA_index_001_peptide_df, AA_index_001_stats_df, AA_index_001_plot_by_run, AA_index_001_plot_by_group,
                AA_index_002_peptide_df, AA_index_002_stats_df, AA_index_002_plot_by_run, AA_index_002_plot_by_group,
                AA_index_003_peptide_df, AA_index_003_stats_df, AA_index_003_plot_by_run, AA_index_003_plot_by_group,
                AA_index_004_peptide_df, AA_index_004_stats_df, AA_index_004_plot_by_run, AA_index_004_plot_by_group,
                AA_index_005_peptide_df, AA_index_005_stats_df, AA_index_005_plot_by_run, AA_index_005_plot_by_group,
                AA_index_006_peptide_df, AA_index_006_stats_df, AA_index_006_plot_by_run, AA_index_006_plot_by_group,
                AA_index_007_peptide_df, AA_index_007_stats_df, AA_index_007_plot_by_run, AA_index_007_plot_by_group,   
                AA_index_008_peptide_df, AA_index_008_stats_df, AA_index_008_plot_by_run, AA_index_008_plot_by_group,
                AA_index_009_peptide_df, AA_index_009_stats_df, AA_index_009_plot_by_run, AA_index_009_plot_by_group,
                AA_index_010_peptide_df, AA_index_010_stats_df, AA_index_010_plot_by_run, AA_index_010_plot_by_group,
                AA_index_011_peptide_df, AA_index_011_stats_df, AA_index_011_plot_by_run, AA_index_011_plot_by_group,
                AA_index_012_peptide_df, AA_index_012_stats_df, AA_index_012_plot_by_run, AA_index_012_plot_by_group,
                AA_index_013_peptide_df, AA_index_013_stats_df, AA_index_013_plot_by_run, AA_index_013_plot_by_group,
                AA_index_014_peptide_df, AA_index_014_stats_df, AA_index_014_plot_by_run, AA_index_014_plot_by_group,
                AA_index_015_peptide_df, AA_index_015_stats_df, AA_index_015_plot_by_run, AA_index_015_plot_by_group,
                AA_index_016_peptide_df, AA_index_016_stats_df, AA_index_016_plot_by_run, AA_index_016_plot_by_group,
                AA_index_017_peptide_df, AA_index_017_stats_df, AA_index_017_plot_by_run, AA_index_017_plot_by_group,
                AA_index_018_peptide_df, AA_index_018_stats_df, AA_index_018_plot_by_run, AA_index_018_plot_by_group,   
                AA_index_019_peptide_df, AA_index_019_stats_df, AA_index_019_plot_by_run, AA_index_019_plot_by_group, 
                AA_index_020_peptide_df, AA_index_020_stats_df, AA_index_020_plot_by_run, AA_index_020_plot_by_group,
                AA_index_user_peptide_df, AA_index_user_stats_df, AA_index_user_plot_by_run, AA_index_user_plot_by_group,
                user_requested_description
                ):
        
        
        self.bold_label_font = ("TkDefaultFont", 10, "bold")
        self.summary_df = summary_df
        self.summary_stats_df = summary_stats_df
        self.protein_and_peptide_figure = protein_and_peptide_figure
        self.PSM_figure = PSM_figure
        self.similarity_df = similarity_df
        self.similarity_stats = similarity_stats
        self.similarity_table = similarity_table
        self.similarity_plot = similarity_plot
        self.AA_df = AA_df
        self.AA_figure = AA_figure
        self.AA_stats_df = AA_stats_df
        self.AA_stats_figure = AA_stats_figure
        self.AA_type_df = AA_type_df
        self.AA_type_figure = AA_type_figure
        self.AA_type_stats_df = AA_type_stats_df
        self.AA_type_stats_figure = AA_type_stats_figure
        self.missed_cleavages_df = missed_cleavages_df
        self.missed_cleavages_stats_df = missed_cleavages_stats_df
        self.missed_cleavages_figure = missed_cleavages_figure
        self.peptide_length_df = peptide_length_df
        self.peptide_length_stats_df = peptide_length_stats_df
        self.peptide_length_by_run_figure = peptide_length_by_run_figure
        self.peptide_length_by_group_figure = peptide_length_by_group_figure
        self.MW_df = MW_df
        self.MW_stats_df = MW_stats_df
        self.MW_by_run_figure = MW_by_run_figure
        self.MW_by_group_figure = MW_by_group_figure
        self.hydrophobicity_df = hydrophobicity_df
        self.hydrophobicity_stats_df = hydrophobicity_stats_df
        self.hydrophobicity_by_run_figure = hydrophobicity_by_run_figure
        self.hydrophobicity_by_group_figure = hydrophobicity_by_group_figure
        self.pI_df = pI_df
        self.pI_stats_df = pI_stats_df
        self.pI_by_run_figure = pI_by_run_figure
        self.pI_by_group_figure = pI_by_group_figure
        self.instability_df = instability_df
        self.instability_stats_df = instability_stats_df
        self.instability_by_run_figure = instability_by_run_figure
        self.instability_by_group_figure = instability_by_group_figure
        self.propensity_df = propensity_df
        self.propensity_stats_df = propensity_stats_df
        self.propensity_by_run_figure = propensity_by_run_figure
        self.propensity_by_group_figure = propensity_by_group_figure
        
        self.AA_index_001_peptide_df = AA_index_001_peptide_df 
        self.AA_index_001_stats_df = AA_index_001_stats_df 
        self.AA_index_001_plot_by_run = AA_index_001_plot_by_run
        self.AA_index_001_plot_by_group = AA_index_001_plot_by_group
        
        self.AA_index_002_peptide_df = AA_index_002_peptide_df 
        self.AA_index_002_stats_df = AA_index_002_stats_df 
        self.AA_index_002_plot_by_run = AA_index_002_plot_by_run
        self.AA_index_002_plot_by_group = AA_index_002_plot_by_group
        
        self.AA_index_003_peptide_df = AA_index_003_peptide_df
        self.AA_index_003_stats_df = AA_index_003_stats_df
        self.AA_index_003_plot_by_run = AA_index_003_plot_by_run
        self.AA_index_003_plot_by_group = AA_index_003_plot_by_group
        
        self.AA_index_004_peptide_df = AA_index_004_peptide_df
        self.AA_index_004_stats_df = AA_index_004_stats_df
        self.AA_index_004_plot_by_run = AA_index_004_plot_by_run
        self.AA_index_004_plot_by_group = AA_index_004_plot_by_group
        
        self.AA_index_005_peptide_df = AA_index_005_peptide_df
        self.AA_index_005_stats_df = AA_index_005_stats_df
        self.AA_index_005_plot_by_run = AA_index_005_plot_by_run
        self.AA_index_005_plot_by_group = AA_index_005_plot_by_group
        
        self.AA_index_006_peptide_df = AA_index_006_peptide_df
        self.AA_index_006_stats_df = AA_index_006_stats_df
        self.AA_index_006_plot_by_run = AA_index_006_plot_by_run
        self.AA_index_006_plot_by_group = AA_index_006_plot_by_group
        
        self.AA_index_007_peptide_df = AA_index_007_peptide_df
        self.AA_index_007_stats_df = AA_index_007_stats_df  
        self.AA_index_007_plot_by_run = AA_index_007_plot_by_run
        self.AA_index_007_plot_by_group = AA_index_007_plot_by_group
        
        self.AA_index_008_peptide_df = AA_index_008_peptide_df  
        self.AA_index_008_stats_df = AA_index_008_stats_df  
        self.AA_index_008_plot_by_run = AA_index_008_plot_by_run
        self.AA_index_008_plot_by_group = AA_index_008_plot_by_group
        
        self.AA_index_009_peptide_df = AA_index_009_peptide_df
        self.AA_index_009_stats_df = AA_index_009_stats_df  
        self.AA_index_009_plot_by_run = AA_index_009_plot_by_run
        self.AA_index_009_plot_by_group = AA_index_009_plot_by_group
        
        self.AA_index_010_peptide_df = AA_index_010_peptide_df
        self.AA_index_010_stats_df = AA_index_010_stats_df
        self.AA_index_010_plot_by_run = AA_index_010_plot_by_run
        self.AA_index_010_plot_by_group = AA_index_010_plot_by_group

        self.AA_index_011_peptide_df = AA_index_011_peptide_df
        self.AA_index_011_stats_df = AA_index_011_stats_df
        self.AA_index_011_plot_by_run = AA_index_011_plot_by_run
        self.AA_index_011_plot_by_group = AA_index_011_plot_by_group

        self.AA_index_012_peptide_df = AA_index_012_peptide_df
        self.AA_index_012_stats_df = AA_index_012_stats_df
        self.AA_index_012_plot_by_run = AA_index_012_plot_by_run
        self.AA_index_012_plot_by_group = AA_index_012_plot_by_group

        self.AA_index_013_peptide_df = AA_index_013_peptide_df
        self.AA_index_013_stats_df = AA_index_013_stats_df
        self.AA_index_013_plot_by_run = AA_index_013_plot_by_run
        self.AA_index_013_plot_by_group = AA_index_013_plot_by_group

        self.AA_index_014_peptide_df = AA_index_014_peptide_df
        self.AA_index_014_stats_df = AA_index_014_stats_df
        self.AA_index_014_plot_by_run = AA_index_014_plot_by_run
        self.AA_index_014_plot_by_group = AA_index_014_plot_by_group

        self.AA_index_015_peptide_df = AA_index_015_peptide_df
        self.AA_index_015_stats_df = AA_index_015_stats_df
        self.AA_index_015_plot_by_run = AA_index_015_plot_by_run
        self.AA_index_015_plot_by_group = AA_index_015_plot_by_group

        self.AA_index_016_peptide_df = AA_index_016_peptide_df
        self.AA_index_016_stats_df = AA_index_016_stats_df
        self.AA_index_016_plot_by_run = AA_index_016_plot_by_run
        self.AA_index_016_plot_by_group = AA_index_016_plot_by_group
        
        self.AA_index_017_peptide_df = AA_index_017_peptide_df
        self.AA_index_017_stats_df = AA_index_017_stats_df
        self.AA_index_017_plot_by_run = AA_index_017_plot_by_run
        self.AA_index_017_plot_by_group = AA_index_017_plot_by_group
        
        self.AA_index_018_peptide_df = AA_index_018_peptide_df
        self.AA_index_018_stats_df = AA_index_018_stats_df
        self.AA_index_018_plot_by_run = AA_index_018_plot_by_run
        self.AA_index_018_plot_by_group = AA_index_018_plot_by_group

        self.AA_index_019_peptide_df = AA_index_019_peptide_df
        self.AA_index_019_stats_df = AA_index_019_stats_df
        self.AA_index_019_plot_by_run = AA_index_019_plot_by_run
        self.AA_index_019_plot_by_group = AA_index_019_plot_by_group

        self.AA_index_020_peptide_df = AA_index_020_peptide_df
        self.AA_index_020_stats_df = AA_index_020_stats_df
        self.AA_index_020_plot_by_run = AA_index_020_plot_by_run
        self.AA_index_020_plot_by_group = AA_index_020_plot_by_group
        
        self.AA_index_user_peptide_df = AA_index_user_peptide_df
        self.AA_index_user_stats_df = AA_index_user_stats_df
        self.AA_index_user_plot_by_run = AA_index_user_plot_by_run
        self.AA_index_user_plot_by_group = AA_index_user_plot_by_group
        self.user_requested_description = user_requested_description
    

        self.root = tk.Tk()
        self.root.title("Results Viewer")
        self.root.title("PRO-VET Data Results")
        self.root.geometry("1200x700")

        notebook = ttk.Notebook(self.root)
        notebook.pack(fill=tk.BOTH, expand=True)

        tab_1 = ttk.Frame(notebook)
        notebook.add(tab_1, text="Results (Page 1)")


        tk.Label(tab_1,
            text="Raw Data", width=25, height = 2, font=self.bold_label_font
            ).grid(row=0, column=1, padx=4, pady=4) 
        
        tk.Label(tab_1,
            text="Means, Stds, and P-Values", width=25, height = 2, font=self.bold_label_font
            ).grid(row=0, column=2, padx=4, pady=4)   
        
        tk.Label(tab_1,
            text="First Figure", width=25, height = 2, font=self.bold_label_font
            ).grid(row=0, column=3, padx=4, pady=4)   
        
        tk.Label(tab_1,
            text="Second Figure", width=25, height = 2, font=self.bold_label_font
            ).grid(row=0, column=4, padx=4, pady=4) 
        
        tk.Label(tab_1,
            text="Proteins/Peptides/PSMs", width=25, height = 2, font=self.bold_label_font
            ).grid(row=1, column=0, padx=4, pady=4)  
        
        tk.Button(tab_1,
            text="Summary Information",
            command=lambda: open_dataframe_window(self.summary_df),
            width=25, height = 2
        ).grid(row=1, column=1, padx=4, pady=4) 
        
        tk.Button(tab_1,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.summary_stats_df),
            width=25, height = 2
        ).grid(row=1, column=2, padx=4, pady=4)   
        
        tk.Button(tab_1,
            text="Protein and Peptide Figure",
            command=lambda: show_figure_in_tk(self.root, self.protein_and_peptide_figure),
            width=25, height = 2
        ).grid(row=1, column=3, padx=4, pady=4) 

        tk.Button(tab_1,
            text="PSM Figure",
            command=lambda: show_figure_in_tk(self.root, self.PSM_figure),
            width=25, height = 2
        ).grid(row=1, column=4, padx=4, pady=4) 

        ##################################################################
        
        tk.Label(tab_1,
            text="Similarity Index", width=25, height = 2, font=self.bold_label_font
            ).grid(row=2, column=0, padx=4, pady=4) 
        
        
        tk.Button(tab_1,
            text="Similarity Index Dataframe",
            command=lambda: open_dataframe_window(self.similarity_df),
            width=25, height = 2
        ).grid(row=2, column=1, padx=4, pady=4)   
        
        tk.Button(tab_1,
            text="Similarity Index Stats",
            command=lambda: open_dataframe_window(self.similarity_stats),
            width=25, height = 2
        ).grid(row=2, column=2, padx=4, pady=4)
         
        tk.Button(tab_1,
            text="Similarity Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.similarity_plot)
        ).grid(row=2, column=3, padx=4, pady=4 )
        
        ###############################################################
        
        tk.Label(tab_1,
            text="AA Composition", width=25, height = 2, font=self.bold_label_font
            ).grid(row=3, column=0, padx=4, pady=4) 
        
        tk.Button(tab_1,
            text="AA Composition Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.AA_df)
        ).grid(row=3, column=1, padx=4, pady=4 ) 

        tk.Button(tab_1,
            text="AA Stats Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.AA_stats_df)
        ).grid(row=3, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="AA Composition Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.AA_figure)
        ).grid(row=3, column=3, padx=4, pady=4 )

        tk.Button(tab_1,
            text="AA Stats Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.AA_stats_figure)
        ).grid(row=3, column=4, padx=4, pady=4 )

        ################################################################
        
        tk.Label(tab_1,
            text="AA Type", width=25, height = 2, font=self.bold_label_font
            ).grid(row=4, column=0, padx=4, pady=4)        
        
        tk.Button(tab_1,
            text="Open AA Type Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.AA_type_df)
        ).grid(row=4, column=1, padx=4, pady=4 )

        tk.Button(tab_1,
            text="AA Type Stats Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.AA_type_stats_df)
        ).grid(row=4, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="AA Type Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.AA_type_figure)
        ).grid(row=4, column=3, padx=4, pady=4 )


        tk.Button(tab_1,
            text="AA Type Stats Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.AA_type_stats_figure)
        ).grid(row=4, column=4, padx=4, pady=4 )
        
        ##############################################

        tk.Label(tab_1,
            text="Missed Cleavages", width=25, height = 2, font=self.bold_label_font
            ).grid(row=5, column=0, padx=4, pady=4)        
        

        tk.Button(tab_1,
            text="Open Missed Cleavages Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.missed_cleavages_df)
        ).grid(row=5, column=1, padx=4, pady=4 )
        
        tk.Button(tab_1,
            text="Missed Cleavages Stats Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.missed_cleavages_stats_df)
        ).grid(row=5, column=2, padx=4, pady=4)

        tk.Button(tab_1,
            text="Missed Cleavages Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.missed_cleavages_figure)
        ).grid(row=5, column=3, padx=4, pady=4)
        
        ##############################################
        tk.Label(tab_1,
            text="Molecular Weight", width=25, height = 2, font=self.bold_label_font,
            ).grid(row=6, column=0, padx=4, pady=4)  

        
        tk.Button(tab_1,
            text="MW Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.MW_df)
        ).grid(row=6, column=1, padx=4, pady=4 )

        tk.Button(tab_1,
            text="MW Stats", width=25, height = 2,
            command=lambda: open_dataframe_window(self.MW_stats_df)
        ).grid(row=6, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="MW by Run Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.MW_by_run_figure)
        ).grid(row=6, column=3, padx=4, pady=4)

        tk.Button(tab_1,
            text="MW by Group Figure", width = 25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.MW_by_group_figure)
        ).grid(row=6, column=4, padx=4, pady=4)

        ##############################################
        
        tk.Label(tab_1,
            text="Peptide Length", width=25, height = 2, font=self.bold_label_font
            ).grid(row=7, column=0, padx=4, pady=4)  


        tk.Button(tab_1,
            text="Peptide Length Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.peptide_length_df)
        ).grid(row=7, column=1, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Peptide Length Stats", width=25, height = 2,
            command=lambda: open_dataframe_window(self.peptide_length_stats_df)
        ).grid(row=7, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Peptide Length by Run Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.peptide_length_by_run_figure)
        ).grid(row=7, column=3, padx=4, pady=4)

        tk.Button(tab_1,
            text="Peptide Length by Group Figure", width = 25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.peptide_length_by_group_figure)
        ).grid(row=7, column=4, padx=4, pady=4)

        ##############################################
        tk.Label(tab_1,
            text="Hydrophobicity", width=25, height = 2, font=self.bold_label_font 
            ).grid(row=8, column=0, padx=4, pady=4)  


        tk.Button(tab_1,
            text="Hydrophobicity Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.hydrophobicity_df)
        ).grid(row=8, column=1, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Hydrophobicity Stats", width=25, height = 2,
            command=lambda: open_dataframe_window(self.hydrophobicity_stats_df)
        ).grid(row=8, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Hydrophobicity by Run Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.hydrophobicity_by_run_figure)
        ).grid(row=8, column=3, padx=4, pady=4)

        tk.Button(tab_1,
            text="Hydrophobicity by Group Figure", width = 25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.hydrophobicity_by_group_figure)
        ).grid(row=8, column=4, padx=4, pady=4)


        ##############################################
        tk.Label(tab_1,
            text="pI (Isoelectric Point)", width=25, height = 2, font=self.bold_label_font 
            ).grid(row=9, column=0, padx=4, pady=4)  

        
        tk.Button(tab_1,
            text="pI Dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.pI_df)
        ).grid(row=9, column=1, padx=4, pady=4 )

        tk.Button(tab_1,
            text="pI Stats", width=25, height = 2,
            command=lambda: open_dataframe_window(self.pI_stats_df)
        ).grid(row=9, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="pI by Run Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.pI_by_run_figure)
        ).grid(row=9, column=3, padx=4, pady=4)

        tk.Button(tab_1,
            text="pI by Group Figure", width = 25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.pI_by_group_figure)
        ).grid(row=9, column=4, padx=4, pady=4)

        ##############################################
        tk.Label(tab_1,
            text="Instability Index", width=25, height = 2, font=self.bold_label_font 
            ).grid(row=10, column=0, padx=4, pady=4)  

        
        tk.Button(tab_1,
            text="Instability Index dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.instability_df)
        ).grid(row=10, column=1, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Instability Index Stats", width=25, height = 2,
            command=lambda: open_dataframe_window(self.instability_stats_df)
        ).grid(row=10, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Instability Index by Run Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.instability_by_run_figure)
        ).grid(row=10, column=3, padx=4, pady=4)

        tk.Button(tab_1,
            text="Instability Index by Group Figure", width = 25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.instability_by_group_figure)
        ).grid(row=10, column=4, padx=4, pady=4)

        ##############################################
        tk.Label(tab_1,
            text="Propensity", width=25, height = 2, font=self.bold_label_font 
            ).grid(row=11, column=0, padx=4, pady=4)  

        
        tk.Button(tab_1,
            text="Propensity dataframe", width=25, height = 2,
            command=lambda: open_dataframe_window(self.propensity_df)
        ).grid(row=11, column=1, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Propensity Stats", width=25, height = 2,
            command=lambda: open_dataframe_window(self.propensity_stats_df)
        ).grid(row=11, column=2, padx=4, pady=4 )

        tk.Button(tab_1,
            text="Propensity by Run Figure", width=25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.propensity_by_run_figure)
        ).grid(row=11, column=3, padx=4, pady=4)

        tk.Button(tab_1,
            text="Propensity by Group Figure", width = 25, height = 2,
            command=lambda: show_figure_in_tk(self.root, self.propensity_by_group_figure)
        ).grid(row=11, column=4, padx=4, pady=4)

        tk.Label(tab_1,
            text="See Tab 2 for more\n data and plots", width=25, height = 2, font=self.bold_label_font 
        ).grid(row=12, column=0, padx=4, pady=4) 




 
        tab_2 = ttk.Frame(notebook)
        notebook.add(tab_2, text="Results (Page 2)")

        tk.Label(tab_2,
            text="All Data on this page Derived from AA index", width=30, height = 2, 
            wraplength = 200, font=self.bold_label_font
            ).grid(row=0, column=0, padx=0, pady=0)


        tk.Label(tab_2,
            text="Raw Data", width=20, height = 2, font=self.bold_label_font
            ).grid(row=0, column=1, padx=0, pady=0) 
        
        tk.Label(tab_2,
            text="Means, Stds, P-Values", width=20, height = 2, font=self.bold_label_font
            ).grid(row=0, column=2, padx=0, pady=0)   
        
        tk.Label(tab_2,
            text="Figure by Sample/Run", width=20, height = 2, font=self.bold_label_font
            ).grid(row=0, column=3, padx=0, pady=0)   
        
        tk.Label(tab_2,
            text="Figure by Group", width=20, height = 2, font=self.bold_label_font
            ).grid(row=0, column=4, padx=0, pady=0) 
        
        tk.Label(tab_2,
            text="Source", width=30, height = 2, font=self.bold_label_font
            ).grid(row=0, column=5, padx=0, pady=0) 
        
        ###############################################  
        tk.Label(tab_2,
            text="alpha-CH chemical shifts", width=30, height = 1, font=self.bold_label_font
            ).grid(row=1, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_001_peptide_df),
            width=20, height = 1
        ).grid(row=1, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_001_stats_df),
            width=20, height = 1
        ).grid(row=1, column=2, padx=0, pady=0)   
        
        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_001_plot_by_run),
            width=20, height = 1
        ).grid(row=1, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_001_plot_by_group),
            width=20, height = 1
        ).grid(row=1, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="ANDN920101 (Andersen et al., 1992)", width=30, height = 1, 
            ).grid(row=1, column=5, padx=0, pady=0)  
        
        ###############################################  
        tk.Label(tab_2,
            text="Membrane-buried pref.", width=30, height = 1, font=self.bold_label_font
            ).grid(row=2, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_002_peptide_df),
            width=20, height = 1
        ).grid(row=2, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_002_stats_df),
            width=20, height = 1
        ).grid(row=2, column=2, padx=0, pady=0)   
        
        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_002_plot_by_run),
            width=20, height = 1
        ).grid(row=2, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_002_plot_by_group),
            width=20, height = 1
        ).grid(row=2, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="ARGP820103 (Argos et al., 1982)", width=30, height = 1, 
            ).grid(row=2, column=5, padx=0, pady=0)  
        
        ###############################################  
        tk.Label(tab_2,
            text="Residue volume", width=30, height = 1, font=self.bold_label_font
            ).grid(row=3, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_003_peptide_df),
            width=20, height = 1
        ).grid(row=3, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_003_stats_df),
            width=20, height = 1
        ).grid(row=3, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_003_plot_by_run),
            width=20, height = 1
        ).grid(row=3, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_003_plot_by_group),
            width=20, height = 1
        ).grid(row=3, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="BIGC670101 (Bigelow, 1967)", width=30, height = 1, 
            ).grid(row=3, column=5, padx=0, pady=0)  
        
        ###############################################  
        tk.Label(tab_2,
            text="Retention coeff in TFA", width=30, height = 1, font=self.bold_label_font
            ).grid(row=4, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_004_peptide_df),
            width=20, height = 1
        ).grid(row=4, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_004_stats_df),
            width=20, height = 1
        ).grid(row=4, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_004_plot_by_run),
            width=20, height = 1
        ).grid(row=4, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_004_plot_by_group),
            width=20, height = 1
        ).grid(row=4, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="BROC820101 (Browne et al., 1982)", width=30, height = 1, 
            ).grid(row=4, column=5, padx=0, pady=0)  
        
        ###############################################  
        tk.Label(tab_2,
            text="Transfer free energy to surface", width=30, height = 1, font=self.bold_label_font
            ).grid(row=5, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_005_peptide_df),
            width=20, height = 1
        ).grid(row=5, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_005_stats_df),
            width=20, height = 1
        ).grid(row=5, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_005_plot_by_run),
            width=20, height = 1
        ).grid(row=5, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_005_plot_by_group),
            width=20, height = 1
        ).grid(row=5, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="BULH740101  (Bull-Breese, 1974)", width=30, height = 1, 
            ).grid(row=5, column=5, padx=0, pady=0)  
        
        ###############################################  
        tk.Label(tab_2,
            text="Free energy of solution in water", width=30, height = 1, font=self.bold_label_font
            ).grid(row=6, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_006_peptide_df),
            width=20, height = 1
        ).grid(row=6, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_006_stats_df),
            width=20, height = 1
        ).grid(row=6, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_006_plot_by_run),
            width=20, height = 1
        ).grid(row=6, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_006_plot_by_group),
            width=20, height = 1
        ).grid(row=6, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="CHAM820102  (Charton-Charton, 1982)", width=30, height = 1, 
            ).grid(row=6, column=5, padx=0, pady=0)  
        ###############################################

        tk.Label(tab_2,
            text="Entropy of formation", width=30, height = 1, font=self.bold_label_font
            ).grid(row=7, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_007_peptide_df),
            width=20, height = 1
        ).grid(row=7, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_007_stats_df),
            width=20, height = 1
        ).grid(row=7, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_007_plot_by_run),
            width=20, height = 1
        ).grid(row=7, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_007_plot_by_group),
            width=20, height = 1
        ).grid(row=7, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="HUTJ700103 (Hutchens, 1970)", width=30, height = 1, 
            ).grid(row=7, column=5, padx=0, pady=0)  

        ###############################################  
        tk.Label(tab_2,
            text="Average accessible surface area", width=30, height = 1, font=self.bold_label_font
            ).grid(row=8, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_008_peptide_df),
            width=20, height = 1
        ).grid(row=8, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_008_stats_df),
            width=20, height = 1
        ).grid(row=8, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_008_plot_by_run),
            width=20, height = 1
        ).grid(row=8, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_008_plot_by_group),
            width=20, height = 1
        ).grid(row=8, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="JANJ780101  (Janin et al., 1978)", width=30, height = 1, 
            ).grid(row=8, column=5, padx=0, pady=0)  

        ###############################################  
        tk.Label(tab_2,
            text="Probability of beta-sheet", width=30, height = 1, font=self.bold_label_font
            ).grid(row=9, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_009_peptide_df),
            width=20, height = 1
        ).grid(row=9, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_009_stats_df),
            width=20, height = 1
        ).grid(row=9, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_009_plot_by_run),
            width=20, height = 1
        ).grid(row=9, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_009_plot_by_group),
            width=20, height = 1
        ).grid(row=9, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="KANM800102 (Kanehisa-Tsong, 1980)", width=30, height = 1, 
            ).grid(row=9, column=5, padx=0, pady=0)  
        
        ###############################################  
        tk.Label(tab_2,
            text="Probability of inner helix", width=30, height = 1, font=self.bold_label_font
            ).grid(row=10, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_010_peptide_df),
            width=20, height = 1
        ).grid(row=10, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_010_stats_df),
            width=20, height = 1
        ).grid(row=10, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_010_plot_by_run),
            width=20, height = 1
        ).grid(row=10, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_010_plot_by_group),
            width=20, height = 1
        ).grid(row=10, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="KANM800103  (Kanehisa-Tsong, 1980)", width=30, height = 1, 
            ).grid(row=10, column=5, padx=0, pady=0)  
        
    ###############################################  
        tk.Label(tab_2,
            text="Zimm-Bragg parameter s at 20 C", width=30, height = 1, font=self.bold_label_font
            ).grid(row=11, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_011_peptide_df),
            width=20, height = 1
        ).grid(row=11, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_011_stats_df),
            width=20, height = 1
        ).grid(row=11, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_011_plot_by_run),
            width=20, height = 1
        ).grid(row=11, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_011_plot_by_group),
            width=20, height = 1
        ).grid(row=11, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="SUEM840101 (Sueki et al., 1984)", width=30, height = 1, 
            ).grid(row=11, column=5, padx=0, pady=0)  
        ###############################################  
        tk.Label(tab_2,
            text="Electron-ion interaction potential", width=30, height = 1, font=self.bold_label_font
            ).grid(row=12, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_012_peptide_df),
            width=20, height = 1
        ).grid(row=12, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_012_stats_df),
            width=20, height = 1
        ).grid(row=12, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_012_plot_by_run),
            width=20, height = 1
        ).grid(row=12, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_012_plot_by_group),
            width=20, height = 1
        ).grid(row=12, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="VELV850101 (Veljkovic et al., 1985)", width=30, height = 1, 
            ).grid(row=12, column=5, padx=0, pady=0)  
        ###############################################  
        tk.Label(tab_2,
            text="Unfolding Gibbs energy in water, pH7.0 ", width=30, height = 1, font=self.bold_label_font
            ).grid(row=13, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_013_peptide_df),
            width=20, height = 1
        ).grid(row=13, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_013_stats_df),
            width=20, height = 1
        ).grid(row=13, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_013_plot_by_run),
            width=20, height = 1
        ).grid(row=13, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_013_plot_by_group),
            width=20, height = 1
        ).grid(row=13, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="YUTK870101 (Yutani et al., 1987)", width=30, height = 1, 
            ).grid(row=13, column=5, padx=0, pady=0)  
        ###############################################  
        tk.Label(tab_2,
            text="Amphiphilicity index", width=30, height = 1, font=self.bold_label_font
            ).grid(row=14, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_014_peptide_df),
            width=20, height = 1
        ).grid(row=14, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_014_stats_df),
            width=20, height = 1
        ).grid(row=14, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_014_plot_by_run),
            width=20, height = 1
        ).grid(row=14, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_014_plot_by_group),
            width=20, height = 1
        ).grid(row=14, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="MITS020101 (Mitaku et al., 2002)", width=30, height = 1, 
            ).grid(row=14, column=5, padx=0, pady=0)  
        ###############################################  
        tk.Label(tab_2,
            text="Retention coefficient at pH 2", width=30, height = 1, font=self.bold_label_font
            ).grid(row=15, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_015_peptide_df),
            width=20, height = 1
        ).grid(row=15, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_015_stats_df),
            width=20, height = 1
        ).grid(row=15, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_015_plot_by_run),
            width=20, height = 1
        ).grid(row=15, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_015_plot_by_group),
            width=20, height = 1
        ).grid(row=15, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="GUOD860101 (Guo et al., 1986)", width=30, height = 1, 
            ).grid(row=15, column=5, padx=0, pady=0)  
        ###############################################  
        tk.Label(tab_2,
            text="Buriability", width=30, height = 1, font=self.bold_label_font
            ).grid(row=16, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_016_peptide_df),
            width=20, height = 1
        ).grid(row=16, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_016_stats_df),
            width=20, height = 1
        ).grid(row=16, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_016_plot_by_run),
            width=20, height = 1
        ).grid(row=16, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_016_plot_by_group),
            width=20, height = 1
        ).grid(row=16, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="ZHOH040103 (Zhou-Zhou, 2004)", width=30, height = 1, 
            ).grid(row=16, column=5, padx=0, pady=0)  
        ###############################################

        tk.Label(tab_2,
            text="Average internal preferences", width=30, height = 1, font=self.bold_label_font
            ).grid(row=17, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_017_peptide_df),
            width=20, height = 1
        ).grid(row=17, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_017_stats_df),
            width=20, height = 1
        ).grid(row=17, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_017_plot_by_run),
            width=20, height = 1
        ).grid(row=17, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_017_plot_by_group),
            width=20, height = 1
        ).grid(row=17, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="OLSK800101 (Olsen, 1980)", width=30, height = 1, 
            ).grid(row=17, column=5, padx=0, pady=0)  
        ###############################################  

        tk.Label(tab_2,
            text="Average eccentricity", width=30, height = 1, font=self.bold_label_font
            ).grid(row=18, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_018_peptide_df),
            width=20, height = 1
        ).grid(row=18, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_018_stats_df),
            width=20, height = 1
        ).grid(row=18, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_018_plot_by_run),
            width=20, height = 1
        ).grid(row=18, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_018_plot_by_group),
            width=20, height = 1
        ).grid(row=18, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="KARS160105 (Karkbara-Knisley, 2016)", width=30, height = 1, 
            ).grid(row=18, column=5, padx=0, pady=0)  
        ###############################################  
        tk.Label(tab_2,
            text="Apparent partition energies calculated from Wertz-Scheraga index ", width=30, height = 1, font=self.bold_label_font
            ).grid(row=19, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_019_peptide_df),
            width=20, height = 1
        ).grid(row=19, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_019_stats_df),
            width=20, height = 1
        ).grid(row=19, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_019_plot_by_run),
            width=20, height = 1
        ).grid(row=19, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_019_plot_by_group),
            width=20, height = 1
        ).grid(row=19, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="GUYH850102 (Guy,1985)", width=30, height = 1, 
            ).grid(row=19, column=5, padx=0, pady=0)  
        ###############################################  
        
        tk.Label(tab_2,
            text="Polarity", width=30, height = 1, font=self.bold_label_font
            ).grid(row=20, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw_data",
            command=lambda: open_dataframe_window(self.AA_index_020_peptide_df),
            width=20, height = 1
        ).grid(row=20, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_020_stats_df),
            width=20, height = 1
        ).grid(row=20, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_020_plot_by_run),
            width=20, height = 1
        ).grid(row=20, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_020_plot_by_group),
            width=20, height = 1
        ).grid(row=20, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text="ZIMJ680103 (Zimmerman et al., 1968)", width=30, height = 1, 
            ).grid(row=20, column=5, padx=0, pady=0)  
             
       ###############################################  
        
        tk.Label(tab_2,
            text="User Selected", width=30, height = 1, font=self.bold_label_font
            ).grid(row=21, column=0, padx=0, pady=0)  
        
        tk.Button(tab_2,
            text="raw data",
            command=lambda: open_dataframe_window(self.AA_index_user_peptide_df),
            width=20, height = 1
        ).grid(row=21, column=1, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Summary Stats",
            command=lambda: open_dataframe_window(self.AA_index_user_stats_df),
            width=20, height = 1
        ).grid(row=21, column=2, padx=0, pady=0)   

        tk.Button(tab_2,
            text="Figure by run",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_user_plot_by_run),
            width=20, height = 1
        ).grid(row=21, column=3, padx=0, pady=0) 

        tk.Button(tab_2,
            text="Figure by Group",
            command=lambda: show_figure_in_tk(self.root, self.AA_index_user_plot_by_group),
            width=20, height = 1
        ).grid(row=21, column=4, padx=0, pady=0) 

        tk.Label(tab_2,
            text=user_requested_description, width=30, height = 1, 
            ).grid(row=21, column=5, padx=0, pady=0)  
             



        





    def run(self):
        self.root.mainloop()


