# load packages
library('tidyr')
library('dplyr')
library('phytools')
library('ggplot2')
library('ggtree')
library('ggnewscale')
library('RColorBrewer')
library('colorspace')

# setwd
setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# import basic metadata 
meta_basic <- read.csv2("input_data/very_basic_meta_.csv", sep = ',')

# import the phylogenetic tree (RaxML SNP tree)
tree <- read.tree("input_data/RAxML_bestTree.raxmltree_all_ref_N16961")

# extract the sames of the samples (these are 9-10 colony picks per patient, thats why they are called 'colonies' in the following)
colonies <- as.data.frame(tree$tip.label)
colnames(colonies) <- 'colony'
# extract the patient number
colonies['Individual'] <- gsub('Vc.*', '', colonies$colony)
colonies <- merge(colonies, meta_basic, by = 'Individual', all.x = T)
colonies <- colonies[,c("colony","Individual","household","ind_contact")]
rownames(colonies) <- colonies$colony
# check if all match
setdiff(tree$tip.label, colonies$colony)

colonies$household <- as.factor(colonies$household)
colourCount <- length(unique(colonies$household))
# define getpallete
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# subset household information to add as a column to tree-heatmap plot
plot <- as.data.frame(colonies[,'household'])
colnames(plot) <- c('household')
rownames(plot) <- colonies$colony

# rotate some branches (for visualisation)
tree_1 <- ape::rotate(tree, node = 643)
tree_2 <- ape::rotate(tree_1, node = 490)

# reroot on the reference
tree_b <- root(tree_2, outgroup = 'NZ_LT907989.1_Vibrio_cholerae_strain_A19')

# basic ggtree
tree.gg_1 <- ggtree(tree_b, layout = 'rectangular', ladderize = FALSE) + geom_treescale()

# add first column (household)
tree_heatmap <- gheatmap(tree.gg_1, plot, width = 0.1,
                         colnames=T, legend_title="") + scale_fill_manual(values = c(getPalette(colourCount), "#9ECAE1",  "#DEEBF7"), na.value = "grey91") + geom_rootedge()
# add new color (fill) scale, such that a new column can be added using a different color scheme
tree_heatmap_1 <- tree_heatmap + new_scale_fill()

# add lineage info
lineage <- read.table('input_data//lineage_assignemnet_from_medaka.txt', header = T)
lineage_plot <- lineage[,c('sample', 'Lineage_group')]
# samples to rownames, such that in gheatmap it can be added to tree
rownames(lineage_plot) <- lineage_plot$sample
lineage_plot$sample <- NULL
lineage_plot$Lineage_group <- factor(lineage_plot$Lineage_group, levels = c('BD2', 'IND1.2', 'IND1.3'))

tree_heatmap_1 <- gheatmap(tree_heatmap_1, lineage_plot, width = 0.1, offset = 0.000011,
                         colnames=T, legend_title="") + scale_fill_manual(values = c("#FDD58C", "#4763AB", "#B2258F"), na.value = "grey91")
# add new color (fill) scale, such that a new column can be added using a different color scheme
tree_heatmap_1 <- tree_heatmap_1 + new_scale_fill()

# import info on how many chromosomes were assembled
QC <- read.csv2('input_data/QC_extended.csv')
QC <- QC[,c('Sample', 'chr')]
rownames(QC) <- QC$Sample
QC$Sample <- NULL

QC$chr <- ifelse(QC$chr %in% c("Two circular chromosomes assembled",   "Not all chromosomes circularised"), "Two chromosomes assembled", "One chromosome assembled")

tree_heatmap_1c <- gheatmap(tree_heatmap_1, QC, offset=0.000022, width = 0.1,
                            colnames=T, legend_title="") + scale_fill_manual(values = c('lightblue', 'darkblue'), na.value = "grey91")
# add new fill scale
tree_heatmap_1cb <- tree_heatmap_1c + new_scale_fill()

# import how many copies of HS1 (aka 'overlap') have been detected per genome
overlap_sum <- read.table('input_data/blastout_HR_all_sum.tab', header = T)
rownames(overlap_sum) <- overlap_sum$Sample
overlap_sum$Sample <- NULL
overlap_sum$on_x_chr <- NULL
overlap_sum$n_HR <- factor(overlap_sum$n_HR)
tree_heatmap_1cc <- gheatmap(tree_heatmap_1cb, overlap_sum, offset=0.000032, width = 0.1,
                            colnames=T, legend_title="") + scale_fill_manual(values = c("#FCAD91", "#F34A35"))
tree_heatmap_1cc


pdf("output_figures/all_tree_heatmap_chr_HR_v2.pdf", height = 10, width = 8)
tree_heatmap_1cc
dev.off()
