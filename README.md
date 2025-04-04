# Prevalent chromosome fusion in Vibrio cholerae O1
This directory contains all scripts used for the analysis in "Prevalent chromosome fusion in Vibrio cholerae O1". 
All input data to reproduce the figures and run the R scripts can be downloaded via the OpenScienceFoundation (https://osf.io/xyfvg/). 

## Abstract 
Two circular chromosomes are a defining feature of the family Vibrionaceae, including the pathogen Vibrio cholerae, with rare reports of isolates with a single, fused chromosome. Here we report chromosome fusions in clinical V. cholerae O1 isolates, including several independent fusion events that are likely transmissible within a household. Fusions occur in a 12 kilobase-pair homologous sequence shared between the two chromosomes and are stable for 200 generations under laboratory conditions and. We found no detectable effect of fusion on V. cholerae growth, virulence factor expression, or biofilm formation. The factors promoting fusion, affecting chromosome stability, and subtle phenotypic or clinical consequences merit further investigation.

## Overview
### Vc_Chromosome_Fusion_Bash_Scripts.Rmd
This includes bash scripts which were run for the analysis comments allowing to follow the workflow. 

### lineage_assignment.R
To assign each of our genomes to one of the sublineages within the V. cholerae 7PET, medaka was used to could the number of high quality SNV to one reference per sublineage. This acript evaluates the output and assigns the lineage for which the lowest number of SNV was called. 

### tree_heatmap_iqtree_UF.R
This script plots a phylogenetic tree (SNV tree), which was previously compiled using IQtree2 and adds the information of household membership, the 7PET sublineage designation, the number of chromosomes identified by sequencing, and the number of times HS1 was detected in the genome.

### plot_QC_ONT.R 
This script evaluates and plots the quality control data for the genomes sequenced for this study. 

### eval_readcount_HR.R
This script summarises how often the HS1 was found in each genome sequenced for this study (screened via blastn). 
It plots the annotation of HS1 and evaluates and plots the number of reads which span HS1 and are flanking different chr sequences excluding the ones where a sequencing adapter was found (Fig. S2). 
Finally, this script evaluates the presence similarity of sequences which are known to play a role in chromosome replication, such as the par genes, the origins of replication (ori 1 and 2), the Chr2 replication triggering Site crtS, and the dam gene. 

### eval_island_flanking_blast.R
This script evaluates the blast of the flanking side of the four virulence associated genomic islands VPI-1, VPI-2, VSP-I and VSP-II and outputs the island positions for each sample.

### circular_chromosome_plot.R
This script plots the non-fused V. cholerae chromosomes and a fused V. cholerae chromosome with the location of the two origins of replication (ori1 and ori2), crtS, HS1 and the four virulence associated genomic islands (VPI-1, VPI-2, VSP-I and VSP-II). 

### passaging_eval_reads.R
We tested the stability of chromosome fusion state by passaging an isolate with a fused and a strain with a non-fused chromosome or 200 generations (20 days).
The two isolates were sequenced on day 0, 3, 8, 11, 16 and 20. To evaluate the fusion state of these isolates, we mapped the reads agains HS1 in its fused states (chr1-HS1-ch2/chr2-HS1-ch1) and its non-fused states (chr1-HS1-ch1/chr2-HS1-ch2). This script plots the number of reads mapping to each state from this isolates. 
We called variants from each passaging experiment isolate to its source isolate using clair3. This script plots the number of high quality SNV identified.

### mut_rec.R
This script evaluates whether there are SNVs in mut or rec genes (involved in homologous recombination) which are specific to fused chromosome genomes. SNVs were called for each genome against the reference genome V. cholerae N16961 using medaka. 

### public_long_reads.R: 
This script summarises QC, outputs from GTDB-tk, kraken / bracken and screens for assemblies with one fused chromosome. Further it summarises the comparisons between chr 1 and chr 2. 

### public_short_read_HR_eval.R:
This plots HS1 and VSP-II copy number with core genome tree 
