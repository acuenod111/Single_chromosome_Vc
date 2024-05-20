# Long-read sequencing reveals prevalent chromosome fusion in Vibrio cholerae O1
This directory contains all scripts used for the analysis in "Long-read sequencing reveals prevalent chromosome fusion in Vibrio cholerae O1". 
All input data to reproduce the figures and run the R scripts can be found in 'Input_data'. 

## Abstract 
Two circular chromosomes are a defining feature of the family Vibrionaceae, including the pathogen V. cholerae, with rare reports of isolates with a single, fused chromosome. Here we report chromosome fusions in clinical V. cholerae isolates, including several independent fusion events stable enough to be transmitted between patients within a household. Fusion occurs in a 12 kilobase pair homologous sequence (HS1) which is shared between the two chromosomes.

## Overview
### Vc_Chromosome_Fusion_Bash_Scripts.Rmd
This includes bash scripts which were run for the analysis comments allowing to follow the workflow. 

### lineage_assignment.R
To assign each of our genomes to one of the sublineages within the V. cholerae 7PET, medaka was used to could the number of high quality SNV to one reference per sublineage. This acript evaluates the output and assigns the lineage for which the lowest number of SNV was called. 

### tree_heatmap.R 
This script plots a phylogenetic tree (SNV tree), which was previously compiled using RAxML and adds the information of household membership, the 7PET sublineage designation, the number of chromosomes identified by sequencing, and the number of times HS1 was detected in the genome (Fig. 2C)

### plot_QC_ONT.R 
This script evaluates and plots the quality control data for the genomes sequenced for this study (Fig. 1C, Fig. S1)

### eval_readcount_HR.R
This script summarises how often the HS1 was found in each genome sequenced for this study (screened via blastn). 
It plots the annotation of HS1 (Fig. 2B) and evaluates and plots the number of reads which span HS1 and are flanking different chr sequences excluding the ones where a sequencing adapter was found (Fig. S2). 

### transitionfinder.R
Performs ancestral state reconstruction and identifies fusion / fission events along the phylogenetic tree (Figure S3). 

### public_long_reads.R: 
This script summarises QC, outputs from GTDB-tk, kraken / bracken and screens for assemblies with one fused chromosome. Further it summarises the comparisons between chr 1 and chr 2 (Fig. S4AB, Fig. S5). 

### public_short_read_HR_eval.R:
This plots HS1 and VSP-II copy number with core genome tree (Fig. S4C)
