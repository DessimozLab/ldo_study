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
This document co

step4_loading_expression_data-rna_seq
-------------------------------------
This document contains the code that downloads and reformats the available expression data from the bgee database and reformat all the expression data used.
The plant expression data was downloaded from: https://expression.plant.tools/

step5_genome_mapping
--------------------
Generate genome mapping tables for each of the species. This is necessary as PANTHER genomes are imported from UniProt RPs, whereas bgee is using ensembl data directly.


step6_tissue_specificity
------------------------


step7_pairwise
--------------
This notebook was used to parse the results of the phmmer search against the metaclust database. The HMM used to search the DB is also included in this folder. ( Source data table 1 )

step8_outgroup
--------------
Data, scripts, and notebooks to generate IMEs clustering based on Jaccard indexes (Supplementary Fig. 8) and  synteny plots (Fig. 4b; Supplementary Fig. 9).

step9_plots
-----------
Data, scripts and structural models to produce the phylogenetic tree with environmental information and selected trimers with electrostatic surfaces to be calculated with APBS (Extended Fig. 8). 




non_rev_rooting
---------------
Non reversible rooting analysis using sequence data to place the root of Fsx1 in either archaea or eukaryota. (Supplementary Fig. 10).

time_acquisition
----------------
Time of aquisition analysis to find the approximate evolutionary period where Fsx1 emerged. ( Supplementary Fig. 15). 





Contains supplementary files produced by bioinformatics analysis, primers, synthesized sequences, movies, sequence identifiers etc.

A note on software and system requirements
----------------------
This software and the scripts used for this project have only been tested on Ubuntu 18.04 and 20.04 environments
At various steps of the project compelementary bioinformatics strategies were used to generate models or search for homology. Furthermore, several external tools are called by these scripts and notebooks. Installation instuctions and dependencies for the software used in this project can be found at the following locations:

HH-Suite:[https://github.com/soedinglab/hh-suite] ( v3.3.0 )
Emboss:[http://emboss.sourceforge.net/download/] ( v6.5.7 )
Blast suite:[https://blast.ncbi.nlm.nih.gov/Blast.cgi] ( v2.12.0 )
Hmmer:[http://hmmer.org/] ( v3.3.2 ) 
OpenMM suite and pdbfixer:[https://openmm.org/] ( v7.6.0 and 1.8.1 respectively )
ClustalO:[http://www.clustal.org/omega/] ( v1.2.2 )
iqtree:[http://www.iqtree.org/] ( v2 )
MAD root:[https://github.com/davidjamesbryant/MADroot] ( v1.0 )
AlphaFold:[https://github.com/deepmind/alphafold] ( v2.1.1 )
TMalign:[https://zhanggroup.org/TM-align/] (Version 20210224)
FATCAT:[https://fatcat.godziklab.org/] ( v2.0 )
I-tasser:[https://zhanggroup.org/I-TASSER/] ( v5.1 )
Modeller:[https://salilab.org/modeller/] ( v10.2 )
TOPCONS:[https://topcons.cbr.su.se/] ( v2.0 )
McScan:[https://github.com/tanghaibao/mcscan] ( v0.8 )
FastME:[http://www.atgc-montpellier.fr/fastme/] ( v2.0 )
rootDigger:[https://github.com/computations/root_digger] ( v1.7.0 )

The metaclust and uniclust30 databases used with Hmmer and HHblits can be found on the mmseqs web page at [https://metaclust.mmseqs.org/] and [https://uniclust.mmseqs.com/] respectively. The Pfam database is available at Pfam:[http://pfam.xfam.org/]( v33 ). 
