# PRO-VET: A Proteomics Software tool for the Visualization and Evaluation of Bottom-up Proteomics Experiments to Monitor Overall Variability.   

Abstract

Proteomics experiments rely on complex, multi-step workflows that include sample solubilization, protein digestion, cysteine reduction and alkylation, peptide separation by reverse-phase nano-liquid chromatography, tandem mass spectrometry, spectral matching, and data processing. Ensuring consistent performance is critical, as variability at any stage can compromise data quality. Although extensive methodologies and software exist to monitor instrument performance, overall workflow evaluation is often limited to coarse metrics such as total number of proteins, peptides, and/or peptide spectral matches. To address this limitation, Proteomics Visualization and Evaluation Tool (Pro-VET), a software platform that provides comprehensive metrics for assessing proteomics workflow performance, was developed. Pro-VET is implemented in Python and distributed as a standalone Windows executable. The software enables comparison of a wide range of peptide properties, including hydrophobicity, missed cleavages, amino acid composition, and peptide length. In addition, Pro-VET supports visualization of any of 566 peptide properties derived from the AAindex database. A novel quantitative metric, similarity index, which enables direct comparison of proteomic runs is also introduced. It compares the peptides detected from the more abundant proteins.  Pro-VET outputs both figures and raw data files, facilitating downstream visualization and analysis. The utility of Pro-VET is demonstrated using three representative case studies.
<img width="3000" height="3000" alt="PRO-VET_image" src="https://github.com/user-attachments/assets/463bcdac-2ddf-4f08-9cbb-097ebba6dbed" />

Submitted to J. of Proteome Research

Installation:
Pro-Vet (Pro-VET_ver_00.40.exe) should be downloaded with two .pdf files and one image file.  The two .pdf files are i) this instructions file (Pro-VET_instructions.pdf) and ii) file describing the peptide properties available for the AAindex and their index numbers (AA_index_details.pdf).  In addition, an image file (PRO-VET_image.jpg) is included.   All four files should be put in a single folder and the resulting results folder and subfolders will be created in this folder.  Only the .exe file is needed for most operations.     

Steps to Use:
1)	Generation of Data
This program is designed to work with proteomics data generated using DDA workflows.  Search instrument generated “raw” files with your favorite search engine.  Each sample/run should be searched independently.   All searches should be performed using the same search engine, search parameters and fasta database(s).  Results from several of the more common search engines/proteomics analysis software suites can be directly imported into Pro-VET (see Appendix 1).  
2)	Data Import.
Prepare data for import into Pro_VET.  For any software, one can manually convert the results into .csv files and import into this software.  For each run, two files are needed.  The first is the peptide file which must contain at least 4 columns, with the headings “seq”, “accession”, “PSM”, and “Area”.  Sequence (“seq”) is the peptide sequence using all upper-case 1-letter abbreviations for the 20 canonical amino acids.  Peptide spectral matches (“PSM”) or spectral count are the number of times that peptide is matched to a spectra.  Accession (“accession”) is the accession number of the protein that matches that peptide sequence.  This must be the same accession number used in the fasta database in the search and there should only one per peptide.  While certain peptides can be derived from more than one protein, its most likely associated protein should be used.  This is termed differently for different search engines.  Area (“Area”) is the peak area for that peptide.  If only peak heights are known, these can be used in place of the areas but the heading should remain “Area”.  For If no peak areas are available, simply put 1.0 into all entries of that column.  (Do not select “area” or “log of area” options during processing but most functions will work properly).  An optional fifth column for modifications, entitled ‘mods’, can be used.  The current version of Pro-VET only examines Met oxidation and should use the format such as “1xOxidation [M5]”.  This indicates that one Met at position 5 is oxidized.  This file should be saved as a .csv (Comma-Separated Values) file with file name underscore (‘_’) then “_peptides.csv”.  For example, “sample_1A_peptides.csv”.  The second file should be the protein file with at least 4 columns with the heading “description”, “accession”, “PSM”, and “prot_len”.  The description (“description”) generally includes that protein name and is how the protein is described.  Accession number (“accession”) number and “PSM” are as described above.  The protein length (“prot_len”) should also be included.  These files should also be .csv files, labelled with sample name, followed by underscore and “_proteins.csv”.  For example, the protein file corresponding to the peptide file above should be “sample_1a_proteins.csv”.  Several manually created files should be available with the distribution of this software, entitled: “manual_input_example_1_peptides.csv”, and “manual_input_example_1_proteins.csv” for sample 1.  Similar files for example 2 and 3 should also be available.  All these column and file names are case sensitive and do NOT include the quotation marks.  These files should be relatively easy to generate in Excel for those familiar with that software and then exported as a .csv file.  

In addition to these protein and peptide files, a fasta file used for these searches with a “.fasta” extension is needed.  This should be identical to the one used during the searches and in the format from uniport.org.  More than 1 fasta file can be used but it must match the fasta file(s) used during the searches.

If the various samples are associated with groups, you can assign each sample to a group.  If groups are used, a minimum of 3 samples per group is recommended.  

The data import page should appear similar to below with your data in the associated fields:
 
 
3)	The General Analysis parameters must be selected.  These are:
a.	Select Processing Method.  
-	Single or Multiple Protein Analysis.  Indicate if you want analysis of a single protein or multiple proteins.  
-	Calculation Method.  For some analysis, one of several calculation methods are used.  The four options are PSM, log 2 of PSM, peak area (or peak height) and log 2 of peak area (or peak height).  These are used for calculating similarity index and displaying similarity values.  
-	Base Run. (Optional).  Base run, or reference sample, is the sample that is used to compare all others.  It will appear first in most analysis.
-	Testing Only.  Do NOT check.  Only for testing.
b.	Parameters for Single Protein Analysis.  (Only used for single protein Analysis.)  
-	Protein Selection Method.  For single protein analysis, there are several options available.  Protein can be selected randomly from one of the top 500 (???) most abundant proteins based on spectral count.  They can be selected by user, which requires user to put protein accession number in entry box.  They can also be selected by rank.  Rank is the rank of the protein based on the number of PSMs for all samples.  The most abundant protein is 1 and they are in decreasing order.  This also requires user to specify the rank in the proper entry box.  (10 is default).  
-	Additional Information.  Add the need information, when needed.  (See Protein Selection Method)

c.	Parameters for Multiple Protein Analysis.  (Only used for multiple protein Analysis.)  
-	N-Value.  Select the number of proteins to include in the analysis.  These will use the most abundant proteins based on overall spectral count.  For example, N-Value of 10 (default) will examine peptides from the 10 most abundant proteins.  For some analysis (e.g. Similarity Index), it is best to avoid proteins with a low spectral count.  Increasing the N-value will also require more processing time and computer resources.
-	Peptide Properties.  Select the peptide properties to examine.  These include:
a.	Peptide Molecular Weight (MW).  Monoisotopic Molecular Weight of peptides. 
b.	Peptide Length.  Number of Amino acids in the peptide.
c.	Similarity Index.  (See manuscript).
d.	AA Composition. (Prevalence of each Amino acid, weighed by number of PSMs)
e.	AA Type.  AAs can be classified in a number of different ways. For this study, the 7 categories used are i) aliphatic (A, G, I, L, P, V), aromatic (F, W, Y), sulfur containing (C, M), Hydroxyl (S, T), basic (R, H, K), acidic (D, E), and amide (N, Q).  
f.	Number of Missed Cleavages.  Missed Cleavages assumes the use of trypsin as the protease and cleavage of the peptide bonds on the C-terminal of lysine (K) and arginine (R) residues, but not before proline (Pro).  
g.	Hydrophobicity (or gravy score).  The Grand average of hydropathicity (GRAVY) of the linear polypeptide sequence is calculated as the sum of hydropathy values of all amino acids, divided by the number of residues in the sequence. Increasing positive score indicates greater hydrophobicity, but no account is taken of the way the protein folds in three dimensions or the percentage of residues buried in the hydrophobic core of the protein. The calculation is based on the Kyte-Doolittle scale [1].
h.	Instability Index.  Based on the instability of presence of certain dipeptide sequences.  Used the method of Guruprasad, et al. [2] as described by Gasteiger, et al. [3].  
i.	Isoelectric Point (pI).  The method described by Kozlowski [4] was used with an option of 4 different pKa values from EMBOSS [5], Lehninger [6], Styer [7], and Bjellqvist [8].  
j.	Helix Propensity.  Propensity of peptide to for either alpha- or pi- helix.
k.	AA index properties.  Twenty (20) are incorporated automatically.  Any other can be specified (see table below).


 
The Parameter Selection Tab:
 


4)	Start Tab.  
-	Results Filename and Folder.  All results including data tables and figures will be placed I the indicated folder in the same folder where the Pro-VET file is located.  These images should be better than those generated on-screen.  Higher quality figures can be created using the created data files.  All dataframes and figures are automatically saved following creation. 
-	Start Processing.  Make sure all other parameters are properly selected prior to clicking the “Start Processing” file.  Once clicked, it may take several minutes to generate results.  Most of the computer run time is used to upload data.  It is highly recommended that one examines a few files samples initially to become familiar with the program.  Several sample files are provided.   
-	Save Log.  If problems develop in sample processing, the log can be saved and it will help in troubleshooting.
 



5)	Miscellaneous Information Tab.  
Self Explanatory.
 



Importing Data from Search Engines.
General Comments:  Run protein search engine and obtain results.  Each raw file (.raw)  should be searched independently.  Data from most common search engines can be used but some may need to be manually manipulate prior to generate data that can be imported into Pro-Vet.  There are others that can be exported or created by search engines and used directly.  In addition, at least one fasta file is required with a “.fasta” extension.  There are many formats possible for fasta files.  This program works best for the fasta format used by uniport (www.uniprot.org).   Other fasta formats should be tested prior to use.
Possible methods to import search results data:   
a.	Manually Generated Files.  For the search results of each raw file, two Comma-Separated Values (.csv) files need to be created representing the proteins and peptides detected.  The protein file should have at least four (4) columns for i) protein (column name = “description”), ii) accession number (column name = “accession”) used during the spectral match search,  iii) the number of spectral matches, also referred to as peptide spectral matches (column name = “PSM”), and iv) protein length (column name = “prot_len”).  The protein files should be names with any identifier followed by an underscore (“_”) and “proteins”.  For example, sample 1, would have a file names of “sample_1_proteins.csv”.
Peptide files should also have at least five (5) columns.  These are for i) sequence (column = “seq”), ii) protein accession number (column name = “accession”), iii) spectral count (“PSM”), and iv) ion intensity (“intensity”) of the peptide and v) oxidation of Methionine (column name = “mods” for modifications).  They should follow the format 1xOxidation for 1 oxidized Met and 2xOxidaton for 2 oxidized Met in a peptide.  If no intensity information is available, please make the intensity equal to 1 for all peptides.  If not Oxidation of methionine is detected, the column is still needed and should be blank (or it can contain other modifications.)  The peptide files should be names with any identifier followed by an underscore (“_”) and “peptide”.  For example, sample 1, would have a file names of “sample_1_peptides.csv”.  

All these column and file names are case sensitive and do NOT include the quotation marks.  These files should be relatively easy to generate in Excel for those familiar with that software and then exported as a .csv file.  Sample data (ex. manual_input_example_peptides.csv and manual_input_example_1_protiens.csv) can be downloaded from the github page. (www.github.?????)  

b.	Maxquant is a proteomics software package developed and updated by Jurgen Cox, Matthias Mann and coworkers for more than a decade  [9,10].  It uses the Andromea search engine, is in widespread use, and is freely available.  The version used to generate the test data was 2.7.3.0.  Pro-Vet has been tested and works with this version.  Each file search generates a “peptides.txt” and “proteinGroups.txt”.  These are typically in the “combined/txt” subdirectory and can be used directly by simply adding a sample identifier followed by an underscore (‘_’) prior to these.  For example, sample 1 would create “sample_1_peptides.txt” and sample_1_proteinGroups.txt”.  
c.	Proteome DiscovererTM is also a widely used proteomics software and is available from Thermo Fisher Scientific.  While it is easy to use and vendor supported, it is not freely available and different versions and add-ons are available.  The version utilized in the development of this software was 2.5.0.400 with Amanda 2.0 search Engine [11] using IMP-apQuant [12] for peptide quantitation.  The Proteins and the Peptide Groups worksheets are exported as Excel workbooks and labeled with sample designation followed by “_Proteins” and “_ProteinGroups”, respectively.  These should have “.xlsx” extensions.  For example, sample 1 data could be denoted as “sample_1_Proteins.xlsx” and “sample_1_PeptideGroups.xlsx”.   
d.	MSFragger search engine and the FragPipe software were created by Nesvizhskii and coworkers [13].  This is also a nice proteomics software suite free to academic users.  The version used and tested for this software was generated in version 23.1 of FragPipe.  Each file generates a protein.tsv and a peptide.tsv.  These tab-separated values (.tsv) files are similar to the .csv files.  These files can also be used directly by adding an identifier followed by an underscore (“_”) and then protein.tsv and peptide.tsv for the protein and peptide file, respectively.  For example, sample 1 would be sample_1_protein.tsv and sample_2_peptide.tsv.  

 
References:
[1]	J. Kyte, R.F. Doolittle, A simple method for displaying the hydropathic character of a protein., J Mol Biol 157 (1982) 105–132. https://doi.org/10.1016/0022-2836(82)90515-0.
[2]	K. Guruprasad, B.V. Reddy, M.W. Pandit, Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence., Protein Eng 4 (1990) 155–161. https://doi.org/10.1093/protein/4.2.155.
[3]	E. Gasteiger, C. Hoogland, A. Gattiker, S. Duvaud, M.R. Wilkins, R.D. Appel, A. Bairoch, Protein Identification and Analysis Tools on the ExPASy Server, in: J.M. Walker (Ed.), The Proteomics Protocols Handbook, Humana Press, Totowa, NJ, 2005: pp. 571–607. https://doi.org/10.1385/1-59259-890-0:571.
[4]	L.P. Kozlowski, IPC - Isoelectric Point Calculator., Biol Direct 11 (2016) 55. https://doi.org/10.1186/s13062-016-0159-9.
[5]	P. Rice, I. Longden, A. Bleasby, EMBOSS: the European Molecular Biology Open Software Suite., Trends Genet 16 (2000) 276–277. https://doi.org/10.1016/s0168-9525(00)02024-2.
[6]	D.L. Nelson, M.M. Cox, Lehninger Principles of Biochemistry, 5th Edition, W.H. Freeman, New York, NY, USA, 2008.
[7]	J.M. Berg, J.L. Tymoczko, G.J. Gatto, L. Stryer, Biochemistry, W. H. Freeman and Company, New York, 2022.
[8]	B. Bjellqvist, B. Basse, E. Olsen, J.E. Celis, Reference points for comparisons of two-dimensional maps of proteins from different human cell types defined in a pH scale where isoelectric points  correlate with polypeptide compositions., Electrophoresis 15 (1994) 529–539. https://doi.org/10.1002/elps.1150150171.
[9]	J. Cox, M. Mann, MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification, Nature Biotechnology 26 (2008) 1367–1372. https://doi.org/10.1038/nbt.1511.
[10]	J. Cox, N. Neuhauser, A. Michalski, R.A. Scheltema, J.V. Olsen, M. Mann, Andromeda: A Peptide Search Engine Integrated into the MaxQuant Environment, J. Proteome Res. 10 (2011) 1794–1805. https://doi.org/10.1021/pr101065j.
[11]	V. Dorfer, P. Pichler, T. Stranzl, J. Stadlmann, T. Taus, S. Winkler, K. Mechtler, MS Amanda, a Universal Identification Algorithm Optimized for High Accuracy Tandem Mass Spectra, J. Proteome Res. 13 (2014) 3679–3684. https://doi.org/10.1021/pr500202e.
[12]	J. Doblmann, F. Dusberger, R. Imre, O. Hudecz, F. Stanek, K. Mechtler, G. Dürnberger, apQuant: Accurate Label-Free Quantification by Quality Filtering, J. Proteome Res. 18 (2019) 535–541. https://doi.org/10.1021/acs.jproteome.8b00113.
[13]	A.T. Kong, F.V. Leprevost, D.M. Avtonomov, D. Mellacheruvu, A.I. Nesvizhskii, MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry-based proteomics., Nat Methods 14 (2017) 513–520. https://doi.org/10.1038/nmeth.4256.
[14]	S. Kawashima, P. Pokarowski, M. Pokarowska, A. Kolinski, T. Katayama, M. Kanehisa, AAindex: amino acid index database, progress report 2008., Nucleic Acids Res 36 (2008) D202-205. https://doi.org/10.1093/nar/gkm998.
 
