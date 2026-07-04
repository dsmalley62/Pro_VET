

__author__ = "David M. Smalley"
__email__  = "dsmalley62@yahoo.com"

import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import *
from tkinter import filedialog, ttk, messagebox
from turtle import left
import customtkinter
import os
from PIL import Image, ImageTk
import sys
from processing import *
from utils import *
from types import SimpleNamespace
import matplotlib.backends.backend_pdf

# %config InlineBackend.figure_format = 'retina'


import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
mpl.rcParams['savefig.dpi'] = 300
# Test
# Set the warning threshold to 50 figures
plt.rcParams['figure.max_open_warning'] = 50

from PyInstaller.utils.hooks import collect_data_files

datas = collect_data_files("aaindex")

import os
import sys

def resource_path(rel):
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, rel)
    return os.path.join(os.path.abspath("."), rel)

# Monkey-patch aaindex expected path
os.environ["AAINDEX_CATEGORY_FILE"] = resource_path(
    "data/aaindex_to_category.txt"
)

# ==========================================================
# MODEL
# ==========================================================

class AppModel:
    """Holds all user inputs and program state."""

    def __init__(self):
        # Tab 1 state
        self.input_type = tk.StringVar(value="Manually Created")

        # Files
        self.fasta_files = []      # list of {"path": str}
        self.peptide_files = []    # list of {"path": str, "group": StringVar}
        self.protein_files = []    # list of {"path": str}

        # Tab 2 state
        self.analysis_type = tk.StringVar(value="multiple_proteins")
        self.processing_method = tk.StringVar(value="random")
        self.calculation_type = tk.StringVar(value="PSM")

        self.N_value = tk.StringVar(value="10")
        self.rank = tk.StringVar(value="10")
        self.accession = tk.StringVar()
        self.base_run_file_path = tk.StringVar()

        self.var1 = tk.IntVar()
        self.var2 = tk.IntVar()

        # Tab 2 Multiprotein Analysis checkboxes:
        self.MW_var = tk.IntVar()
        self.peptide_length_var = tk.IntVar()
        self.similarity_var = tk.IntVar()
        self.AA_composition_var = tk.IntVar()
        self.AA_type_var = tk.IntVar()
        self.missed_cleavages_var = tk.IntVar()
        self.hydrophobicity_var = tk.IntVar()
        self.instability_var = tk.IntVar()
        
        self.pI_var = tk.IntVar()
        self.helix_propensity_var = tk.IntVar()
        self.pKa_dict = tk.StringVar(value="EMBOSS")
        self.helix_type = tk.StringVar(value="alpha")

        self.AA_index_001_var = tk.IntVar()
        self.AA_index_002_var = tk.IntVar()
        self.AA_index_003_var = tk.IntVar()
        self.AA_index_004_var = tk.IntVar()
        self.AA_index_005_var = tk.IntVar()
        self.AA_index_006_var = tk.IntVar()
        self.AA_index_007_var = tk.IntVar()
        self.AA_index_008_var = tk.IntVar()
        self.AA_index_009_var = tk.IntVar()
        self.AA_index_010_var = tk.IntVar()
        self.AA_index_011_var = tk.IntVar()
        self.AA_index_012_var = tk.IntVar()
        self.AA_index_013_var = tk.IntVar()
        self.AA_index_014_var = tk.IntVar()
        self.AA_index_015_var = tk.IntVar()
        self.AA_index_016_var = tk.IntVar()
        self.AA_index_017_var = tk.IntVar()
        self.AA_index_018_var = tk.IntVar()
        self.AA_index_019_var = tk.IntVar()
        self.AA_index_020_var = tk.IntVar()
        self.AA_index_999_var = tk.IntVar()
        self.AA_index_number_var = tk.StringVar(value="FAUJ830101")



        
        # Tab 3 state
        self.results_file_name = tk.StringVar(value="Results_File_Name")


# ==========================================================
# VIEW
# ==========================================================

class AppView:
    """Constructs UI but contains NO logic."""

    def __init__(self, root):
        self.root = root

        root.title("PRO-VET Data Processing Tool")
        root.geometry("1200x700")

        self.notebook = ttk.Notebook(root)
        self.notebook.pack(expand=1, fill="both")

        # Tabs
        self.tab_input = ttk.Frame(self.notebook)
        self.tab_general = ttk.Frame(self.notebook)
        self.tab_start = ttk.Frame(self.notebook)
        self.tab_misc = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_input, text="Data Input")
        self.notebook.add(self.tab_general, text="Processing Parameters")
        self.notebook.add(self.tab_start, text="Start")
        self.notebook.add(self.tab_misc, text="Miscellaneous Information")

        # Build UI sections
        self.build_tab1()
        self.build_tab2()
        self.build_tab3()
        self.build_tab4()

    ###############################################################
    # TAB 1 — DATA INPUT
    ###############################################################

    def build_tab1(self):
        top = tk.Frame(self.tab_input, bd=2, relief="groove", bg="darkblue", height=180)
        top.pack(fill="x", padx=5, pady=5)

        self.build_tab1_top(top)

        bottom = tk.Frame(self.tab_input, bd=2, relief="groove", bg="lightgreen")
        bottom.pack(fill="both", expand=True, padx=5, pady=5)

        self.build_file_inputs(bottom)

    def build_tab1_top(self, parent):
        # Left: Logo
        left = tk.Frame(parent, width=250, bg="white", bd=1, relief="solid")
        left.pack(side="left", padx=3, pady=3)

        try:
            img = Image.open('./Pro-Vet_image.jpg')
            img = img.resize((200, 150), Image.LANCZOS)
            tk_img = ImageTk.PhotoImage(img)

            label = tk.Label(left, image=tk_img)
            label.image = tk_img
            label.pack(pady=10)
        except:
            Label(left, text="Logo missing").pack()

        # Middle: input type selection
        middle = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        middle.pack(side="left", fill="both", expand=True, padx=3, pady=3)

        Label(middle, text="Select Source of Input Data",
              bg="lightgrey", font=("Helvetica", 14)).pack(pady=20)

        self.input_type_combobox = customtkinter.CTkComboBox(
            middle,
            values=["Manually Created", "MaxQuant", "Proteome Discoverer", "MSFragger"],
            width=250,
            height=30,
            font=("Helvetica", 14),
            state="readonly"
        )
        self.input_type_combobox.pack(pady=20)
        self.input_type_combobox.set("Manually Created") 

        # Right: Instructions button
        right = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        right.pack(side="left", fill="both", expand=True, padx=3, pady=3)

        self.instructions_button = Button(right, text="Open Instructions", font=("Helvetica", 14))
        self.instructions_button.pack(pady=20)

        Label(right, text="Instructions (.pdf format)",
              bg="lightgrey", font=("Helvetica", 14)).pack()

    ###############################################################
    # File panels
    ###############################################################

    def build_file_inputs(self, parent):
        # FASTA panel
        self.fasta_panel = self._make_file_panel(parent, "FASTA Files", ["*.fasta", "*.fa"])
        self.fasta_panel.pack(side="left", fill="both", expand=True, padx=3, pady=3)

        # Peptide panel
        self.peptide_panel = self._make_file_panel(parent, "Peptide Files", ["*.csv", "*.tsv", '*.txt'], numbered=True)
        self.peptide_panel.pack(side="left", fill="both", expand=True, padx=3, pady=3)

        # Protein panel
        self.protein_panel = self._make_file_panel(parent, "Protein Files", ["*.csv", "*.tsv", '*.txt'])
        self.protein_panel.pack(side="left", fill="both", expand=True, padx=3, pady=3)

    def _make_file_panel(self, parent, title, filetypes, numbered=False):
        panel = tk.Frame(parent, bd=1, relief="solid", bg="white")
        Label(panel, text=title, font=("Helvetica", 14), bg="white").pack(pady=5)

        browse_button = tk.Button(panel, text=f"Browse {title}")
        browse_button.pack(pady=5)

        # Scrollable area
        canvas = tk.Canvas(panel)
        scrollbar = ttk.Scrollbar(panel, orient="vertical", command=canvas.yview)

        scroll_frame = tk.Frame(canvas)
        scroll_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        panel.canvas = canvas
        panel.scroll_frame = scroll_frame
        panel.browse_button = browse_button
        panel.numbered = numbered
        panel.filetypes = filetypes

        return panel

    ###############################################################
    # TAB 2 — Processing Parameters
    ###############################################################

    def build_tab2(self):
        sect_A = tk.Frame(self.tab_general, bd=2, relief="groove", bg="darkblue", height=200)
        sect_A.pack(fill="x", padx=5, pady=0)
        self.build_tab2_sect_A(sect_A)

        sect_B = tk.Frame(self.tab_general, bd=2, relief="groove", bg="lightgrey", height=200)
        sect_B.pack(fill="x", padx=5, pady=0)
        self.build_tab2_sect_B(sect_B)

        sect_C = tk.Frame(self.tab_general, bd=2, relief="groove", bg="darkblue", height=200)
        sect_C.pack(fill="x", padx=5, pady=0)
        self.build_tab2_sect_C(sect_C)

        sect_D = tk.Frame(self.tab_general, bd=2, relief="groove", bg="lightgrey", height=200)
        sect_D.pack(fill="x", padx=5, pady=0)
        self.build_tab2_sect_D(sect_D)

        sect_E = tk.Frame(self.tab_general, bd=2, relief="groove", bg="darkblue", height=200)
        sect_E.pack(fill="x", padx=5, pady=0)
        self.build_tab2_sect_E(sect_E)

        sect_F = tk.Frame(self.tab_general, bd=2, relief="groove", bg="lightgrey", height=200)
        sect_F.pack(fill="x", padx=5, pady=0)
        self.build_tab2_sect_F(sect_F)

    def build_tab2_sect_A(self, parent):
        top = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        top.pack(side="top", padx=3, pady=3, fill="both")
        Label(top, text="Select Processing Method:",
              bg="lightgrey", font=("Helvetica", 14)).pack(anchor="n", padx=5)

    def build_tab2_sect_B(self, parent):
        # Center the four sections
        parent.grid_columnconfigure((0, 1, 2, 3), weight=1)
        parent.grid_rowconfigure(0, weight=0)
        
        
        # -------- Left Frame --------
        left = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left.grid(row=0, column=0, padx=5, pady=5)

        tk.Label(
            left,
            text="General Analysis Parameters:",
            bg="lightgrey",
            font=("Helvetica", 10)
        ).grid(row=0, column=0, sticky="n", padx=5, pady=0)
        self.radio_single_protein = tk.Radiobutton(left, text="Single Protein", bg="lightgrey")
        self.radio_multiple_proteins = tk.Radiobutton(left, text="Multiple Proteins", bg="lightgrey")

        self.radio_single_protein.grid(row=1, column=0, sticky="w", padx=10)
        self.radio_multiple_proteins.grid(row=2, column=0, sticky="w", padx=10)
        
        # -------- left Frame (2 x 2 radios) --------
        left = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left.grid(row=0, column=1, padx=5, pady=0, sticky="n")

        # -------- middle Frame --------
        middle = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        middle.grid(row=0, column=1, padx=5, pady=0)


        tk.Label(
            middle,
            text="Select Calculation Method:",
            bg="lightgrey",
            font=("Helvetica", 10)
        ).grid(row=0, column=0, columnspan=2, sticky="n", padx=5, pady=2)
        self.radio_psm = tk.Radiobutton(middle, text="PSM", bg="lightgrey")
        self.radio_log_psm = tk.Radiobutton(middle, text="log2 PSM", bg="lightgrey")
        self.radio_area = tk.Radiobutton(middle, text="Intensity", bg="lightgrey")
        self.radio_log_area = tk.Radiobutton(middle, text="log2 Intensity", bg="lightgrey")

        self.radio_psm.grid(row=1, column=0, sticky="n", padx=10, pady=2)
        self.radio_log_psm.grid(row=1, column=1, sticky="n", padx=10, pady=2)
        self.radio_area.grid(row=2, column=0, sticky="n", padx=10, pady=2)
        self.radio_log_area.grid(row=2, column=1, sticky="n", padx=10, pady=2)


        # -------- Right Frame --------
        right = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        right.grid(row=0, column=2, padx=5, pady=0)

        Label(right, text="Select Base Run:", bg="lightgrey",
              font=("Helvetica", 10)).pack(anchor="w", padx=5, pady=5)

        self.button_base_run = Button(right, text="Select Base Run")
        self.button_base_run.pack(pady=5)

        self.label_base_run = Label(right, text="No file selected", bg="white", wraplength=200)
        self.label_base_run.pack(pady=5)


        # -------- Far Right Frame --------

        far_right = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        far_right.grid(row=0, column=3, padx=5, pady=0)

        Label(far_right, text="Testing Only:", bg="lightgrey",
              font=("Helvetica", 10)).pack(pady=0)  

        # Test_Checkbox
        self.checkbox_test = tk.Checkbutton(far_right, text="Single Protein", bg="lightgrey", onvalue=1, offvalue=0)
        self.checkbox_test.pack(pady=5)

        self.checkbox_test_2 = tk.Checkbutton(far_right, text="Multiple Proteins", bg="lightgrey", onvalue=1, offvalue=0)
        self.checkbox_test_2.pack(pady=5)

    def build_tab2_sect_C(self, parent):
        
        middle = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        middle.pack(side="top", padx=3, pady=3, fill="both")
        Label(middle, text="Parameters for Single Protein Analysis:",
              bg="lightgrey", font=("Helvetica", 14)).pack(anchor="n", padx=5)

    def build_tab2_sect_D(self, parent):
        # Center the sections
        parent.grid_columnconfigure((0, 1, 2, 3), weight=1)
        parent.grid_rowconfigure(0, weight=0)
        
        
        # -------- Left Frame --------
        left = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left.grid(row=0, column=1, padx=5, pady=5)


        Label(left, text="Select Protein Selection Method:",
              bg="lightgrey", font=("Helvetica", 10)).pack(anchor="w", padx=5)

        self.radio_accession = Radiobutton(left, text="Specific Protein", bg="lightgrey")
        self.radio_random = Radiobutton(left, text="Random Protein", bg="lightgrey")
        self.radio_rank = Radiobutton(left, text="Choose by Rank", bg="lightgrey")
    

        for rb in (self.radio_accession, self.radio_random, self.radio_rank):
            rb.pack(anchor="w", padx=10)

        # Middle: text inputs
        mid = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        mid.grid(row=0, column=2, padx=5, pady=5)
        

        Label(mid, text="Additional Information:", bg="lightgrey",
              font=("Helvetica", 10)).grid(row=0, column=0, columnspan=2, pady=5)

        self.entry_accession = Entry(mid, width=10)
        self.entry_rank = Entry(mid, width=10)


        labels = ["Accession:", "Rank:"]
        entries = [self.entry_accession, self.entry_rank]

        for i, (lbl, ent) in enumerate(zip(labels, entries)):
            Label(mid, text=lbl, bg="lightgrey").grid(row=i+1, column=1, padx=5, pady=3)
            ent.grid(row=i+1, column=0, padx=5, pady=3)

        # -------- Right Frame --------
        right = tk.Frame(parent, bg="lightgrey", bd=0, relief="solid")
        right.grid(row=0, column=3, padx=5, pady=5)

        text_A = "Insert accession number if you would like to analyze a specific protein."
        tk.Label(right,text= text_A, bg="lightgrey", justify="left").pack(anchor="w", pady=5)


        text_B = "Rank refers to the relative protein abundance with 1 being the most abundant protein." 
        tk.Label(right,text= text_B, bg="lightgrey", justify="left").pack(anchor="w", pady=0)




        
    def build_tab2_sect_E(self, parent):
        
        middle = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        middle.pack(side="top", padx=3, pady=3, fill="both")
        Label(middle, text="Parameters for Multiple Protein Analysis:",
              bg="lightgrey", font=("Helvetica", 14)).pack(anchor="n", padx=5)
        
    def build_tab2_sect_F(self, parent):
        container = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        container.pack(padx=5, pady=0)

        # -----------Multiprotein Analysis Data N_value_input-----------
       
        g00 = tk.Frame(container, bg="lightgrey")
        g00.grid(row=0, column=0, padx=5, pady=5, sticky="n")

        Label(g00, text="N value:", bg="lightgrey", font=("Helvetica", 10)).pack(pady=0)
        self.entry_N = Entry(g00, width=10)
        self.entry_N.pack(pady=5)
        
        
        
        
        # ---------- Multiprotein Analysis checkbox 01 ----------
        g01 = tk.Frame(container, bg="lightgrey")
        g01.grid(row=0, column=1, padx=5, pady=0, sticky="n")

        self.checkbox_MW = tk.Checkbutton(g01, text="MW (molecular Weight)", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_MW.pack(pady=5)

        # ---------- Multiprotein Analysis checkbox 02 ----------
        g02 = tk.Frame(container, bg="lightgrey")
        g02.grid(row=0, column=2, padx=5, pady=5, sticky="n")
        self.checkbox_peptide_length = tk.Checkbutton(g02, text="Peptide Length", onvalue=1, offvalue=0,
                                                    bg="lightgrey")
        self.checkbox_peptide_length.pack(pady=0)

        # ---------- Multiprotein Analysis checkbox 03 ----------
        g03 = tk.Frame(container, bg="lightgrey")
        g03.grid(row=0, column=3, padx=5, pady=0, sticky="n")

        self.checkbox_similarity = tk.Checkbutton(g03, text="Similarity Index", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_similarity.pack(pady=0)


        # ---------- Multiprotein Analysis checkbox 04 ----------
        g04 = tk.Frame(container, bg="lightgrey")
        g04.grid(row=0, column=4, padx=5, pady=0, sticky="n")

        self.checkbox_AA_composition = tk.Checkbutton(g04, text="AA Composition", onvalue=1, offvalue=0, 
                                                      bg="lightgrey")
        self.checkbox_AA_composition.pack(pady=0)

        # ---------- Multiprotein Analysis checkbox 05 ----------
        g05 = tk.Frame(container, bg="lightgrey")
        g05.grid(row=0, column=5, padx=5, pady=0, sticky="n")

        self.checkbox_AA_type = tk.Checkbutton(g05, text="AA Type", onvalue=1, offvalue=0, 
                                              bg="lightgrey")
        self.checkbox_AA_type.pack(pady=0)

        # ---------- Multiprotein Analysis checkbox 06 ----------
        g06 = tk.Frame(container, bg="lightgrey")
        g06.grid(row=0, column=6, padx=5, pady=0, sticky="n")

        self.checkbox_missed_cleavages = tk.Checkbutton(g06, text="Missed Cleavages", onvalue=1, offvalue=0, 
                                                        bg="lightgrey")
        self.checkbox_missed_cleavages.pack(pady=0)

        # ---------- Multiprotein Analysis checkbox 07 ----------
        g07 = tk.Frame(container, bg="lightgrey")
        g07.grid(row=0, column=7, padx=5, pady=0, sticky="n")


        self.checkbox_hydrophobicity = tk.Checkbutton(g07, text="Hydrophobicity", onvalue=1, offvalue=0, 
                                                      bg="lightgrey")
        self.checkbox_hydrophobicity.pack(pady=0)

        # ---------- Multiprotein Analysis checkbox 08 ----------
        g08 = tk.Frame(container, bg="lightgrey")
        g08.grid(row=0, column=8, padx=5, pady=0, sticky="n")

        self.checkbox_instability = tk.Checkbutton(g08, text="Instability", onvalue=1, offvalue=0, 
                                                  bg="lightgrey")
        self.checkbox_instability.pack(pady=0)
        
        # ---------- Multiprotein Analysis checkbox 10 ----------
        g10 = tk.Frame(container, bg="lightgrey")
        g10.grid(row=1, column=0, padx=5, pady=0, sticky="n")

        self.checkbox_pI = tk.Checkbutton(g10, text="pI (pKa values)", onvalue=1, offvalue=0, bg="lightgrey")
        self.checkbox_pI.pack(pady=0)

        self.radio_pI_EMBOSS = Radiobutton(g10, text="EMBOSS", bg="lightgrey")
        self.radio_pI_EMBOSS.pack(anchor="w")
        self.radio_pI_Bjellqvist = Radiobutton(g10, text="Bjellqvist", bg="lightgrey")
        self.radio_pI_Bjellqvist.pack(anchor="w")
        self.radio_pI_Lehninger = Radiobutton(g10, text="Lehninger", bg="lightgrey")
        self.radio_pI_Lehninger.pack(anchor="w")
        self.radio_pI_Stryer = Radiobutton(g10, text="Stryer", bg="lightgrey")
        self.radio_pI_Stryer.pack(anchor="w")

        # ---------- Multiprotein Analysis checkbox 11 Set ----------
        g11 = tk.Frame(container, bg="lightgrey")
        g11.grid(row=1, column=1, padx=5, pady=5, sticky="n")

        self.checkbox_helix_propensity = tk.Checkbutton(g11, text="Helix Propensity", onvalue=1, offvalue=0, 
                                                        bg="lightgrey")
        self.checkbox_helix_propensity.pack(pady=0)

        self.radio_helix_alpha = Radiobutton(g11, text="alpha", bg="lightgrey")
        self.radio_helix_alpha.pack(anchor="w")
        self.radio_helix_pi = Radiobutton(g11, text="pi", bg="lightgrey")
        self.radio_helix_pi.pack(anchor="w")

        # ---------- Multiprotein Analysis checkbox 12 Set ----------
        g12 = tk.Frame(container, bg="lightgrey")
        g12.grid(row=1, column=2, padx=5, pady=0, sticky="n")

        self.checkbox_AA_index_001 = tk.Checkbutton(g12, text="AA index prop 1", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_001.pack(pady=0)

        self.checkbox_AA_index_002 = tk.Checkbutton(g12, text="AA index prop 2", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_002.pack(pady=0)

        self.checkbox_AA_index_003 = tk.Checkbutton(g12, text="AA index prop 3", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_003.pack(pady=0)

        self.checkbox_AA_index_004 = tk.Checkbutton(g12, text="AA index prop 4", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_004.pack(pady=0)

        self.checkbox_AA_index_005 = tk.Checkbutton(g12, text="AA index prop 5", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_005.pack(pady=0)



        # ---------- Multiprotein Analysis checkbox 13 ----------
        g13 = tk.Frame(container, bg="lightgrey")
        g13.grid(row=1, column=3, padx=5, pady=0, sticky="n")

        self.checkbox_AA_index_006 = tk.Checkbutton(g13, text="AA index prop 6", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_006.pack(pady=0)

        self.checkbox_AA_index_007 = tk.Checkbutton(g13, text="AA index prop 7", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_007.pack(pady=0)

        self.checkbox_AA_index_008 = tk.Checkbutton(g13, text="AA index prop 8", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_008.pack(pady=0)

        self.checkbox_AA_index_009 = tk.Checkbutton(g13, text="AA index prop 9", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_009.pack(pady=0)

        self.checkbox_AA_index_010 = tk.Checkbutton(g13, text="AA index prop 10", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_010.pack(pady=0)

        # ---------- Multiprotein Analysis checkbox 14 Set ----------
        g14 = tk.Frame(container, bg="lightgrey")
        g14.grid(row=1, column=4, padx=5, pady=0, sticky="n")

        self.checkbox_AA_index_011 = tk.Checkbutton(g14, text="AA index prop 11", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_011.pack(pady=0)

        self.checkbox_AA_index_012 = tk.Checkbutton(g14, text="AA index prop 12", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_012.pack(pady=0)

        self.checkbox_AA_index_013 = tk.Checkbutton(g14, text="AA index prop 13", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_013.pack(pady=0)

        self.checkbox_AA_index_014 = tk.Checkbutton(g14, text="AA index prop 14", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_014.pack(pady=0)

        self.checkbox_AA_index_015 = tk.Checkbutton(g14, text="AA index prop 15", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_015.pack(pady=0)



        # ---------- Multiprotein Analysis checkbox 15 ----------
        g15 = tk.Frame(container, bg="lightgrey")
        g15.grid(row=1, column=5, padx=5, pady=0, sticky="n")

        self.checkbox_AA_index_016 = tk.Checkbutton(g15, text="AA index prop 16", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_016.pack(pady=0)

        self.checkbox_AA_index_017 = tk.Checkbutton(g15, text="AA index prop 17", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_017.pack(pady=0)

        self.checkbox_AA_index_018 = tk.Checkbutton(g15, text="AA index prop 18", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_018.pack(pady=0)

        self.checkbox_AA_index_019 = tk.Checkbutton(g15, text="AA index prop 19", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_019.pack(pady=0)

        self.checkbox_AA_index_020 = tk.Checkbutton(g15, text="AA index prop 20", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_020.pack(pady=0)

        # ---------- Multiprotein Analysis checkbox 16 ----------
        g16 = tk.Frame(container, bg="lightgrey")
        g16.grid(row=1, column=6, padx=5, pady=0, sticky="n")

        self.checkbox_AA_index_999_var = tk.Checkbutton(g16, text="User Selected", onvalue=1, offvalue=0,
                                                 bg="lightgrey")
        self.checkbox_AA_index_999_var.pack(pady=0)

        Label(g16, text="Enter AA Index #", bg="lightgrey", font=("Helvetica", 10)).pack(pady=0)
        self.AA_index_number_var = Entry(g16, width=15)
        self.AA_index_number_var.pack(pady=5)

        
        # ---------- Button for AA Index Details ----------
        g17 = tk.Frame(container, bg="lightgrey")
        g17.grid(row=1, column=7, padx=5, pady=0, sticky="n")   

        self.AA_index_button = Button(g17, text="AA Index Details", font=("Helvetica", 12))
        self.AA_index_button.pack(pady=15)

        Label(g17, text="AA Index Details \n(.pdf format)",
              bg="lightgrey", font=("Helvetica", 14)).pack()










    ###############################################################
    # TAB 3 — START PROCESSING  (REBUILT CLEANLY)
    ###############################################################

    def build_tab3(self):
        top = tk.Frame(self.tab_start, bd=2, relief="groove", bg="darkblue")
        top.pack(fill="x", padx=5, pady=5)

        middle = tk.Frame(self.tab_start, bd=2, relief="groove", bg="darkblue")
        middle.pack(fill="both", expand=True, padx=5, pady=5)

        bottom = tk.Frame(self.tab_start, bd=2, relief="groove", bg="darkblue")
        bottom.pack(fill="x", padx=5, pady=5)

        self.build_tab3_top(top)
        self.build_tab3_middle(middle)
        self.build_tab3_bottom(bottom)

    ###############################################################
    # TAB 3 TOP — File naming
    ###############################################################
    def build_tab3_top(self, parent):
        left = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left.pack(side="left", padx=3, pady=3, fill="both")

        Label(left, text="Results Filename and Folder",
              bg="lightgrey", font=("Helvetica", 16)).pack(anchor="w", padx=5, pady=5)

        self.entry_results = Entry(left, width=20)
        self.entry_results.pack(pady=10)

    ###############################################################
    # TAB 3 MIDDLE — Console Output + Start Button
    ###############################################################
    def build_tab3_middle(self, parent):

        left = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left.pack(side="left", expand=True, fill="both", padx=10, pady=10)

        Label(left, text="Start Processing", bg="lightgrey",
              font=("Helvetica", 14)).pack(pady=10)

        # START BUTTON
        self.button_start = Button(left, text="Start Processing",
                                   font=("Helvetica", 16))
        self.button_start.pack(pady=10)

        # CONSOLE
        console_frame = tk.Frame(left)
        console_frame.pack(expand=True, fill="both", padx=5, pady=5)

        console_text = tk.Text(console_frame, wrap='word',
                               height=20, width=90,
                               bg='black', fg='white',
                               font=('Consolas', 10))
        console_text.pack(side="left", fill="both", expand=True)

        scrollbar = tk.Scrollbar(console_frame, command=console_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        console_text.config(yscrollcommand=scrollbar.set)

        # Redirect stdout
        sys.stdout = ConsoleRedirector(console_text)

    ###############################################################
    # TAB 3 BOTTOM — Save Log button (bottom-left)
    ###############################################################
    def build_tab3_bottom(self, parent):
        left = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left.pack(side="left", padx=3, pady=3)

        # BOTTOM-LEFT SAVE LOG BUTTON (as requested)
        self.button_save_log = Button(left, text="Save Log", font=("Helvetica", 12))
        self.button_save_log.pack(padx=5, pady=5, anchor="w")

    ###############################################################
    # TAB 4 — Miscellaneous Information
    ###############################################################

    def build_tab4(self):
        left = tk.Frame(self.tab_misc, bd=2, relief="groove", bg="darkblue")
        left.grid(row=0, column=0, padx=5, pady=5, sticky="n")

        right = tk.Frame(self.tab_misc, bd=2, relief="groove", bg="darkblue")
        right.grid(row=0, column=1, padx=5, pady=5, sticky="n")

        self.build_tab4_left(left)
        self.build_tab4_right(right)



    
    def build_tab4_left(self, parent):
        left_1 = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left_1.pack(side="top", padx=3, pady=3, fill="both")

        license_text = "MIT License\n  \nCopyright (c) 2026 David M. Smalley\n\nPermission is hereby granted, free of charge, to any person obtaining a copy\nof this software and associated documentation files (the \"Software\"), to deal\nin the Software without restriction, including without limitation the rights\nto use, copy, modify, merge, publish, distribute, sublicense, and/or sell\ncopies of the Software, and to permit persons to whom the Software is\nfurnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\nOUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\nSOFTWARE."

        Label(left_1, text=license_text,
              bg="lightgrey", font=("Helvetica", 10)).pack(anchor="w", padx=5, pady=5)
        
        left_2 = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left_2.pack(side="top", padx=3, pady=3, fill="both")

        reference_text = "Please Cite: \nSmalley, D.M., PRO-VET: A Proteomics Software tool for the \
            Visualization and Evaluation of DDA-Proteomics Experiments to Assess Overall Variability. \
                  Journal of Proteome Research 2026, Submitted" 
        Label(left_2, text= reference_text,
              bg="lightgrey", font=("Helvetica", 10),  wraplength=400).pack(anchor="n", padx=5, pady=5)
        
        left_3 = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        left_3.pack(side="top", padx=3, pady=3, fill="both")
        
        Label(left_3, text="Contact Information:\nDavid M. Smalley\nEmail:dsmalley62@yahoo.com\nORCID: 0000-0002-1234-5678",
              bg="lightgrey", font=("Helvetica", 10),  wraplength=400).pack(anchor="n", padx=5, pady=5)
        
        
        
        Label(left_3, text="Please Feel Free to Contact me with any Comments, Suggestions, Questions or Issues.",
              bg="lightgrey", font=("Helvetica", 10)).pack(anchor="n", padx=5, pady=5)
        

        

    def build_tab4_right(self, parent):
        
        right = tk.Frame(parent, bg="lightgrey", bd=1, relief="solid")
        right.pack(side="right", padx=3, pady=3, fill="both")

        try:
            img = Image.open('./images/Pro-Vet_2.jpg')
            img = img.resize((500, 500), Image.LANCZOS)
            tk_img = ImageTk.PhotoImage(img)

            label = tk.Label(right, image=tk_img)
            label.image = tk_img
            label.pack(pady=10)
        except:
            Label(right, text="Logo missing").pack()
    

    def build_tab4_bottom(self, parent):
        pass
    
      
    
    

        


# ==========================================================
# CONTROLLER
# ==========================================================

class AppController:
    def __init__(self, model: AppModel, view: AppView):
        self.model = model
        self.view = view

        ################################################
        # Wire buttons
        ################################################
        self.view.instructions_button.config(command=self.open_instructions_pdf)
        self.view.AA_index_button.config(command=self.open_AA_index_pdf)

        # File browsing
        self.view.fasta_panel.browse_button.config(
            command=lambda: self.add_files("fasta")
        )
        self.view.peptide_panel.browse_button.config(
            command=lambda: self.add_files("peptide")
        )
        self.view.protein_panel.browse_button.config(
            command=lambda: self.add_files("protein")
        )

        # Base run
        self.view.button_base_run.config(command=self.select_base_run)

        # Start processing
        self.view.button_start.config(command=self.start_processing)

        # Save log (bottom-left)
        self.view.button_save_log.config(command=self.save_log_to_results_folder)

        # Testing checkbox binding
        self.view.checkbox_test.config(variable=self.model.var1)
        self.view.checkbox_test_2.config(variable=self.model.var2)

        # Multiprotein Analysis checkboxes
        self.view.checkbox_MW.config(variable=self.model.MW_var)
        self.view.checkbox_peptide_length.config(variable=self.model.peptide_length_var)
        self.view.checkbox_similarity.config(variable=self.model.similarity_var)
        self.view.checkbox_AA_composition.config(variable=self.model.AA_composition_var)    
        self.view.checkbox_AA_type.config(variable=self.model.AA_type_var)  
        self.view.checkbox_missed_cleavages.config(variable=self.model.missed_cleavages_var)
        self.view.checkbox_hydrophobicity.config(variable=self.model.hydrophobicity_var)    
        self.view.checkbox_instability.config(variable=self.model.instability_var)  
        self.view.checkbox_pI.config(variable=self.model.pI_var)
        self.view.checkbox_helix_propensity.config(variable=self.model.helix_propensity_var)    
        self.view.checkbox_AA_index_001.config(variable=self.model.AA_index_001_var)
        self.view.checkbox_AA_index_002.config(variable=self.model.AA_index_002_var)
        self.view.checkbox_AA_index_003.config(variable=self.model.AA_index_003_var)
        self.view.checkbox_AA_index_004.config(variable=self.model.AA_index_004_var)
        self.view.checkbox_AA_index_005.config(variable=self.model.AA_index_005_var)
        self.view.checkbox_AA_index_006.config(variable=self.model.AA_index_006_var)
        self.view.checkbox_AA_index_007.config(variable=self.model.AA_index_007_var)
        self.view.checkbox_AA_index_008.config(variable=self.model.AA_index_008_var)
        self.view.checkbox_AA_index_009.config(variable=self.model.AA_index_009_var)
        self.view.checkbox_AA_index_010.config(variable=self.model.AA_index_010_var)
        self.view.checkbox_AA_index_011.config(variable=self.model.AA_index_011_var)
        self.view.checkbox_AA_index_012.config(variable=self.model.AA_index_012_var)
        self.view.checkbox_AA_index_013.config(variable=self.model.AA_index_013_var)
        self.view.checkbox_AA_index_014.config(variable=self.model.AA_index_014_var)
        self.view.checkbox_AA_index_015.config(variable=self.model.AA_index_015_var)
        self.view.checkbox_AA_index_016.config(variable=self.model.AA_index_016_var)
        self.view.checkbox_AA_index_017.config(variable=self.model.AA_index_017_var)
        self.view.checkbox_AA_index_018.config(variable=self.model.AA_index_018_var)
        self.view.checkbox_AA_index_019.config(variable=self.model.AA_index_019_var)
        self.view.checkbox_AA_index_020.config(variable=self.model.AA_index_020_var)
        self.view.checkbox_AA_index_999_var.config(variable=self.model.AA_index_999_var)



        ################################################
        # Bind model variables to view widgets
        ################################################
        
        self.view.radio_single_protein.config(variable=self.model.analysis_type, value="single_protein")
        self.view.radio_multiple_proteins.config(variable=self.model.analysis_type, value="multiple_proteins")
        
        self.view.radio_accession.config(variable=self.model.processing_method, value="accession_input")
        self.view.radio_random.config(variable=self.model.processing_method, value="random")
        self.view.radio_rank.config(variable=self.model.processing_method, value="rank_input")


        self.view.radio_psm.config(variable=self.model.calculation_type, value="PSM")
        self.view.radio_log_psm.config(variable=self.model.calculation_type, value="log_2_PSM")
        self.view.radio_area.config(variable=self.model.calculation_type, value="intensity")
        self.view.radio_log_area.config(variable=self.model.calculation_type, value="log_2_intensity")

        self.view.radio_pI_EMBOSS.config(variable=self.model.pKa_dict, value="EMBOSS")
        self.view.radio_pI_Bjellqvist.config(variable=self.model.pKa_dict, value="Bjellqvist")
        self.view.radio_pI_Lehninger.config(variable=self.model.pKa_dict, value="Lehninger")
        self.view.radio_pI_Stryer.config(variable=self.model.pKa_dict, value="Stryer")   
        
        self.view.radio_helix_alpha.config(variable=self.model.helix_type, value="alpha")
        self.view.radio_helix_pi.config(variable=self.model.helix_type, value="pi")



        # Bind entries
        self.view.entry_accession.config(textvariable=self.model.accession)
        self.view.entry_rank.config(textvariable=self.model.rank)
        self.view.entry_N.config(textvariable=self.model.N_value)
        self.view.AA_index_number_var.config(textvariable = self.model.AA_index_number_var)

        # Tab 2 results filename entry
        # self.view.entry_results_param.config(textvariable=self.model.results_file_name)

        # Tab 3 results filename entry (same source)
        self.view.entry_results.config(textvariable=self.model.results_file_name)

        # Combobox
        self.view.input_type_combobox.configure(command=lambda _: self.model.input_type.set(
            self.view.input_type_combobox.get()
        ))

    ###############################################################
    # SAVE LOG
    ###############################################################
    def save_log_to_results_folder(self):
        """Save console output to the results folder."""
        log_text = sys.stdout.get_log()

        folder_name = self.model.results_file_name.get().strip()
        if folder_name == "":
            messagebox.showerror("Error", "Please enter a results filename first.")
            return

        results_path = os.path.join("results", folder_name)
        os.makedirs(results_path, exist_ok=True)

        log_path = os.path.join(results_path, "log.txt")

        try:
            with open(log_path, "w") as f:
                f.write(log_text)
            messagebox.showinfo("Success", f"Log saved to:\n{log_path}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    
    

    ###############################################################
    # FILE BROWSING
    ###############################################################

    def add_files(self, kind):
        """Adds files to model & updates UI."""
        if kind == "fasta":
            types = [("FASTA files", "*.fasta *.fa")]
            storage = self.model.fasta_files
            panel = self.view.fasta_panel

        elif kind == "peptide":
            types = [("Peptide Tables", "*.*")]
            storage = self.model.peptide_files
            panel = self.view.peptide_panel

        else:
            types = [("Protein Tables", "*.*")]
            storage = self.model.protein_files
            panel = self.view.protein_panel

        paths = filedialog.askopenfilenames(filetypes=types)
        if not paths:
            return

        storage.clear()

        for p in paths:
            if kind == "peptide":
                storage.append({"path": p, "group": tk.StringVar()})
            else:
                storage.append({"path": p})

        self.update_file_display(panel, storage, panel.numbered)

    def update_file_display(self, panel, storage, numbered):
        """Refresh file list display."""
        frame = panel.scroll_frame

        for w in frame.winfo_children():
            w.destroy()

        if not storage:
            Label(frame, text="No files selected").pack(pady=5)
            return

        row = 0
        if numbered:
            Label(frame, text="File", font=("Arial", 10, "bold")).grid(row=row, column=0)
            Label(frame, text="Group", font=("Arial", 10, "bold")).grid(row=row, column=1)
            row += 1

        for item in storage:
            fname = os.path.basename(item["path"])
            Label(frame, text=fname).grid(row=row, column=0, sticky="w")

            if numbered:
                Entry(frame, width=5, textvariable=item["group"]).grid(row=row, column=1)

            row += 1

    ###############################################################
    # BASE RUN SELECTION
    ###############################################################

    def old_select_base_run(self):
        path = filedialog.askopenfilename()
        if path:
            self.model.base_run_file_path.set(path)
            self.view.label_base_run.config(text=path)
        else:
            self.view.label_base_run.config(text="No file selected")


    def select_base_run(self):
        path = filedialog.askopenfilename()
        if path:
            self.model.base_run_file_path.set(path)          # full path stored
            filename = os.path.basename(path)                # just the file name
            self.view.label_base_run.config(text=filename)   # display only filename
        else:
            self.view.label_base_run.config(text="No file selected")

    ###############################################################
    # INSTRUCTIONS PDF
    ###############################################################

    def open_instructions_pdf(self):
        pdf_path = "Pro-VET_instructions.pdf"
        try:
            if os.name == 'nt':
                os.startfile(pdf_path)
            else:
                os.system(f"open {pdf_path}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    ###############################################################
    # AA INDEX PDF
    ###############################################################

    def open_AA_index_pdf(self):
        pdf_path = "AA_Index_Details.pdf"
        try:
            if os.name == 'nt':
                os.startfile(pdf_path)
            else:
                os.system(f"open {pdf_path}")
        except Exception as e:
            messagebox.showerror("Error", str(e))





    ###############################################################
    # COLLECT FINAL INPUTS
    ###############################################################

    def collect_all_inputs(self):
        """Collect everything from model (clean MVC style)."""

        data = {
            "input_type": self.model.input_type.get(),
            "analysis_type": self.model.analysis_type.get(),
            "processing_method": self.model.processing_method.get(),
            "calculation_type": self.model.calculation_type.get(),
            "N_value": self.model.N_value.get(),
            "rank": self.model.rank.get(),
            "accession": self.model.accession.get(),
            "results_file_name": self.model.results_file_name.get(),
            "base_run_file_path": self.model.base_run_file_path.get(),

            "fasta_files": [f["path"] for f in self.model.fasta_files],
            "protein_files": [f["path"] for f in self.model.protein_files],
            "peptide_files": [
                {"path": f["path"], "group": f["group"].get()}
                for f in self.model.peptide_files
            ],
            "testing_single": self.model.var1.get(),
            "testing_multiple": self.model.var2.get(),
            "MW_var": self.model.MW_var.get(),
            "peptide_length_var": self.model.peptide_length_var.get(),
            "similarity_var": self.model.similarity_var.get(),  
            "AA_composition_var": self.model.AA_composition_var.get(),  
            "AA_type_var": self.model.AA_type_var.get(),
            "missed_cleavages_var": self.model.missed_cleavages_var.get(),
            "hydrophobicity_var": self.model.hydrophobicity_var.get(),
            "instability_var": self.model.instability_var.get(),
            "pI_var": self.model.pI_var.get(),
            "helix_propensity_var": self.model.helix_propensity_var.get(), 
            "pKa_dict": self.model.pKa_dict.get(),
            "helix_type": self.model.helix_type.get(),
            "AA_index_001_var": self.model.AA_index_001_var.get(),
            "AA_index_002_var": self.model.AA_index_002_var.get(),
            "AA_index_003_var": self.model.AA_index_003_var.get(),
            "AA_index_004_var": self.model.AA_index_004_var.get(),
            "AA_index_005_var": self.model.AA_index_005_var.get(),
            "AA_index_006_var": self.model.AA_index_006_var.get(),
            "AA_index_007_var": self.model.AA_index_007_var.get(),
            "AA_index_008_var": self.model.AA_index_008_var.get(),
            "AA_index_009_var": self.model.AA_index_009_var.get(),
            "AA_index_010_var": self.model.AA_index_010_var.get(),
            "AA_index_011_var": self.model.AA_index_011_var.get(),
            "AA_index_012_var": self.model.AA_index_012_var.get(),
            "AA_index_013_var": self.model.AA_index_013_var.get(),
            "AA_index_014_var": self.model.AA_index_014_var.get(),
            "AA_index_015_var": self.model.AA_index_015_var.get(),
            "AA_index_016_var": self.model.AA_index_016_var.get(),
            "AA_index_017_var": self.model.AA_index_017_var.get(),
            "AA_index_018_var": self.model.AA_index_018_var.get(),
            "AA_index_019_var": self.model.AA_index_019_var.get(),
            "AA_index_020_var": self.model.AA_index_020_var.get(),
            "AA_index_999_var": self.model.AA_index_999_var.get(),
            "AA_index_number_var": self.model.AA_index_number_var.get(),

        }

        return data

    ###############################################################
    # START PROCESSING + AUTO-SAVE LOG
    ###############################################################

    def start_processing(self):
        print("Starting processing...")
        data = self.collect_all_inputs()

        print("Converting to SimpleNamespace…")
        app = SimpleNamespace(**data)

        print(f"app = {app}")

        # Run user’s processing function
        initial_processing(app)

        print("Processing finished.")

        # === AUTO-SAVE LOG AFTER PROCESSING ===
        print("Point D: ", app.results_path, "results_file_name",app.results_file_name)
        self.save_log_to_results_folder(app)

# ==========================================================
# CONSOLE REDIRECTOR
# ==========================================================

class ConsoleRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget
        self.buffer = []  # stores log text

    def write(self, message):
        self.buffer.append(message)
        self.text_widget.insert(tk.END, message)
        self.text_widget.see(tk.END)

    def flush(self):
        pass

    def get_log(self):
        return "".join(self.buffer)


# ==========================================================
# MODAL APPLICATION DIALOG
# ==========================================================

class ProVetDialog:
    """Wraps the full MVC Tkinter application into a modal dialog."""

    def __init__(self):
        self.root = customtkinter.CTk()
        self.root.title("PRO-VET Data Processing Tool")

        # Build MVC
        self.model = AppModel()
        self.view = AppView(self.root)
        self.controller = AppController(self.model, self.view)

        # Cancel button at bottom
        self.button_cancel = Button(self.root, text="Cancel",
                                    font=("Helvetica", 14),
                                    command=self.cancel)
        self.button_cancel.pack(side="bottom", pady=10)

        self.result = None

        # Make modal
        self.root.grab_set()

    def cancel(self):
        # Restore stdout so terminal works normally
        sys.stdout = sys.__stdout__
        self.result = None
        self.root.quit()
        self.root.destroy()

    def run(self):
        self.root.mainloop()

        # Always restore stdout
        sys.stdout = sys.__stdout__

        # If data was collected, return it
        if hasattr(self.controller, "data"):
            self.result = self.controller.data

        return self.result


# ==========================================================
# ENTRY POINT FUNCTION
# ==========================================================

def run_provet_dialog():
    dialog = ProVetDialog()
    return dialog.run()


# ==========================================================
# OPTIONAL DIRECT EXECUTION
# ==========================================================

if __name__ == "__main__":
    result = run_provet_dialog()
    print("\nDialog returned:\n", result)


    # Restore stdout after window closes
    sys.stdout = sys.__stdout__




