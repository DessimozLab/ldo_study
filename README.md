# LDO project Repo
------------------
The python scripts contained in this repo were used to calculate the branch length significant differences between paralogs and to identify the asymmetric evolution. This method was applied to all gene families across the Tree of Life in the PANTHER database. There is also code used to analyse the gene structures and gene expression profiles. Each notebook contains a first explanatory markdown cell and comments in the code to help users replicate the analysis.

step1_expected_branches
-----------------------
Use to compute the expected branch lengths, using a simple evolutionary model. This gives us a method to identify the unexpectedly long branches, so that we can test the hypothesis of the Least Diverged Orthologue (LDO).

The imput data was downloaded from the Panther database:

#### data/panther-18.0/trees/ -> contains all trees from panther

1. download: wget http://data.pantherdb.org/ftp/panther_library/18.0/PANTHER18.0_hmmscoring.tgz

2. extract only the tree files: tar -zxvf PANTHER18.0_hmmscoring.tgz target/famlib/rel/PANTHER18.0_altVersion/hmmscoring/PANTHER18.0/books/PTHR*/tree.tree

3. rename to trees/PTHR*.tree

#### data/panther-18.0/species_tree.nhx -> species tree
downloaded using panther api with: scripts/panther_species_tree.py


step2_expectedness_of_duplications
----------------------------------
Use to filter the branches to only duplication events.
This is performed by computing the difference from the expected branch length and then translating this into z-scores. Branches can then be classified into six categories (p<0.05): 

i) normal-normal: both branches are not significantly different

ii) short-short: both branches are significantly shorter than expected

iii) long-long: both branches are significantly longer than expected

iv) normal-long: only one branch is significantly longer than expected

v) short-normal: only one branch is significantly shorter than expected

vi) short-long: one branch is significantly shorter than expected and the other is significantly longer than expected.

Also generates Fig 1, Supp Fig 1, Supp Fig 2

Additionally, identifies the genes for the outgroup test.


step3_structure_data
--------------------
Use to download the structure data and compute the structural alignment.


step4_loading_expression_data-rna_seq
-------------------------------------
This document contains the code that downloads and reformats the available expression data from the bgee database (https://www.bgee.org/) and reformat all the expression data used.

The plant expression data was downloaded from: https://expression.plant.tools/


step5_genome_mapping
--------------------
Generate genome mapping tables for each of the species. This is necessary as PANTHER genomes are imported from UniProt RPs, whereas bgee is using ensembl data directly.


step6_pairwise
--------------
This notebook contains the code to run the inparalogue pairwise Pearson's correlation and tissue specificity ($\tau$) tests.


step7_tissue_specificity
------------------------
This notebook contains the code to compute tissue specicity scores ($\tau$). To do this the TPM data is first transformed using the `arcsinh` function -- $\textrm{arcsinh}(x) := \ln (x + \sqrt{x^2 + 1})$ -- before taking the mean value of any replicates.


step8_outgroup
--------------
This notebook identifies relevant species to use for each branch in the species tree as outgroup species.

After, LDO / MDO are compared to the outgroup gene with both a PCC and tau analysis.


step9_plots
-----------
This notebook contains all the code used to analyze the data and generate the plots presented in the paper.


lib/
----
This folder contains the modules used to parse the panther trees in step_1.


general_scripts/
----------------
Scripts used to analyse the data and download specific datasets.


structure_scripts/
------------------
Scripts used in step_3.


A note on software and system requirements
------------------------------------------
The scripts used for this project have only been tested on Ubuntu 24.04 environment.

At various steps several external tools are called by these scripts and notebooks. 

Installation instuctions and dependencies for the software used in this project can be found at the following locations:

foldseek:[https://github.com/steineggerlab/foldseek] (v8.ef4e960)
