library('phytools')
library('ggplot2')
library('ggtree')
library('data.table')
library('dplyr')
library('tidyr')
library('treeio')
library('phangorn')
library('stringr')
library('ggnewscale')
library('RColorBrewer')

setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# import metadata and QC for publicly available QC data
incl_meta <- read.csv2('input_data/public/short_reads/metadata_public_short_read_incl.csv', sep=';')

## import information if the read depth of the HS1 and the VSP-I region and the entire genome
read_depth_HR <- read.table('input_data/public/short_reads/all_public_short_read_HR_average_coverage_1.tab', header = T)
read_depth_HR['Sample'] <- gsub('_average_coverage.*', '', read_depth_HR$name_file)
# add read depth of HR2
read_depth_HR2 <- read.table('input_data/public/short_reads/all_public_short_read_both_HR_average_coverage_1.tab', header = T)
read_depth_HR2['Sample'] <- gsub('_average_coverage_.*', '', read_depth_HR2$name_file)
# merge
read_depth_HR <- merge(read_depth_HR[,c("Average_coverage_full_genome","Average_coverage_of_HR_region","Sample")], read_depth_HR2[,c("Average_coverage_full_genome","Average_coverage_of_HR_region","Average_coverage_of_HR2_region","Sample")], by = c("Average_coverage_full_genome","Average_coverage_of_HR_region","Sample"), all=T)
# convert to numeric
read_depth_HR$Average_coverage_full_genome <- as.numeric(as.character(read_depth_HR$Average_coverage_full_genome))
read_depth_HR$Average_coverage_of_HR_region <- as.numeric(as.character(read_depth_HR$Average_coverage_of_HR_region))
read_depth_HR$Average_coverage_of_HR2_region <- as.numeric(as.character(read_depth_HR$Average_coverage_of_HR2_region))

#calculate the ratio of the read depth in the regions of interest compared to the rest of the genome as a proxy for copy number
read_depth_HR['ratio_HR1_full'] <- read_depth_HR$Average_coverage_of_HR_region / read_depth_HR$Average_coverage_full_genome
read_depth_HR['ratio_HR2_full'] <- read_depth_HR$Average_coverage_of_HR2_region / read_depth_HR$Average_coverage_full_genome

# check before merge
setdiff(incl_meta$Sample, read_depth_HR$Sample)
#merge
read_depth_HR_meta <- merge(read_depth_HR, incl_meta, by = 'Sample', all.x=F, all.y=T)

read_depth_HR_meta$Average_coverage_full_genome <- as.numeric(as.character(read_depth_HR_meta$Average_coverage_full_genome))
read_depth_HR_meta$Read_depth <- as.numeric(as.character(read_depth_HR_meta$Read_depth))
ggplot(read_depth_HR_meta, aes(x=Read_depth, y=Average_coverage_full_genome)) + geom_point()

# I manually checked, the read depth calculation in the illumina QC script is slightly off (although strongly correlated with the one calculated here), the one I get seems correct
read_depth_HR_meta$Read_depth <- NULL
read_depth_HR_meta$Read_quality <- NULL

# in 14 samples, no HR was found, therefore set ratio to 0
read_depth_HR_meta$ratio_HR1_full <- ifelse(is.na(read_depth_HR_meta$ratio_HR1_full), 0, read_depth_HR_meta$ratio_HR1_full)
read_depth_HR_meta$ratio_HR2_full <- ifelse(is.na(read_depth_HR_meta$ratio_HR2_full), 0, read_depth_HR_meta$ratio_HR2_full)

read_depth_HR_meta['HR1_ratio_round'] <- round(read_depth_HR_meta$ratio_HR1_full)
# summarise the ones which were >2
read_depth_HR_meta['HR1_ratio_round'] <- factor(ifelse(read_depth_HR_meta$HR1_ratio_round > 2, '>2', read_depth_HR_meta$HR1_ratio_round), levels = c('0', '1', '2', '>2'))
read_depth_HR_meta['HR2_ratio_round'] <- round(read_depth_HR_meta$ratio_HR2_full)
# summarise the ones which were >2
read_depth_HR_meta['HR2_ratio_round'] <- factor(ifelse(read_depth_HR_meta$HR2_ratio_round > 2, '>2', read_depth_HR_meta$HR2_ratio_round), levels = c('0', '1', '2', '>2'))
# convert year to numeric
read_depth_HR_meta$Year <- as.numeric(as.character(read_depth_HR_meta$Year)) 

# add info whether 7PET or not
Is7PET <- read.table('input_data/public/short_reads/Is_it_7PET.txt', header = T)
setdiff(read_depth_HR_meta$Sample, Is7PET$Sample)
read_depth_HR_meta <- merge(read_depth_HR_meta, Is7PET, by = 'Sample', all.x = T, all.y = F)

# import lineage assignment
snippy_out <- read.table('input_data/public/short_reads/all_snps_1.txt', header = T, row.names = NULL)
colnames(snippy_out) <- c('variant', 'Nr', 'sample', 'reference')
# only count snv
snippy_out <- snippy_out[snippy_out$variant == 'Variant-SNP',]
snippy_out$Nr <- as.numeric(snippy_out$Nr)
snippy_out$reference <- gsub('.fna', '', snippy_out$reference)
# import which ref is which lineage
ref_lin <- read.table('input_data/sel_lineage_assignment_refs.csv', sep=';', header = T)
#merge
snippy_out <- merge(snippy_out, ref_lin[,c('Accession_merge', 'Lineage_group')], by.x='reference', by.y='Accession_merge', all.x=T, all.y=F)

# use the lineage assignment of the reference to which the least number of SNV were detected
snippy_out_sum <- snippy_out %>% 
  group_by(sample, Nr) %>% # summarise ties
  mutate(group_sum = str_c(sort(Lineage_group), collapse="_")) %>% 
  ungroup() %>%
  dplyr::select(sample, Nr, group_sum) %>%
  group_by(sample) %>%
  slice_min(Nr, with_ties =  F)

colnames(snippy_out_sum) <- c('Sample', 'Nr_SNP_to_lineage_ref', 'Lineage')
# check before merge
setdiff(read_depth_HR_meta$Sample, snippy_out_sum$Sample)
# merge
read_depth_HR_meta <- merge(read_depth_HR_meta, snippy_out_sum, by = 'Sample', all.x = T, all.y = F)

# plot on tree
# import tree
nwk <- 'input_data/public/short_reads/short_read_public_fasttree.tree'
tree <- read.tree(nwk)
# draw basic tree
p <- ggtree(tree, layout = 'rectangular') 
# add information which belong to 7PET
read_depth_HR_meta
tree_p <- p %<+% read_depth_HR_meta + 
  geom_tippoint(aes(label = Is_7PET, colour = Is_7PET)) + scale_color_brewer(palette="Accent", na.value="grey91") + geom_treescale()
# find 7PET MRCA
tree2 <- ggtree::groupOTU(tree, read_depth_HR_meta$Sample[read_depth_HR_meta$Is_7PET == 'Yes'], group_name = "group")
# color the non-7PET part of the tree grey
p <- ggtree(tree2, layout = 'rectangular', aes(color=group)) + scale_color_manual(values = c('grey', 'black'))


#read_depth_HR_meta_plot <- read_depth_HR_meta[,c('Sample', 'Is_7PET','Lineage', 'Year', 'HR_ratio_round')]
#read_depth_HR_meta_plot_7p <- read_depth_HR_meta[,c('Sample', 'Is_7PET')]
#rownames(read_depth_HR_meta_plot_7p) <- read_depth_HR_meta_plot_7p$Sample
#read_depth_HR_meta_plot_7p$Sample <- NULL

#tree_heatmap_1a <- gheatmap(p, read_depth_HR_meta_plot_7p, offset=0.00, width = 0.15,
#                            colnames=T, legend_title="") 
#library(ggnewscale)
#tree_heatmap_1b <- tree_heatmap_1a + new_scale_fill()

# add lineage
read_depth_HR_meta['lineage_7P'] <- ifelse(read_depth_HR_meta$Is_7PET != 'Yes', NA, read_depth_HR_meta$Lineage)
# simplify the ones where there are ties (are too many, gets difficult to read). if eg BD1_IND is assigned, choose BD1 (choose BD over IND and IND over T)
unique(read_depth_HR_meta$lineage_7P[grepl('_', read_depth_HR_meta$lineage_7P)])
read_depth_HR_meta$lineage_7P <- gsub("BD1.1_IND1.1", 'BD1.1', read_depth_HR_meta$lineage_7P)
read_depth_HR_meta$lineage_7P <- gsub("BD2_IND2" , 'BD2', read_depth_HR_meta$lineage_7P)
read_depth_HR_meta$lineage_7P <- gsub("BD1.2_IND1.3" , 'BD1.2', read_depth_HR_meta$lineage_7P)
read_depth_HR_meta$lineage_7P <- gsub("BD1.1_IND1" , 'BD1.1', read_depth_HR_meta$lineage_7P)
read_depth_HR_meta$lineage_7P <- gsub("IND1.3_T13" , 'IND1.3', read_depth_HR_meta$lineage_7P)

# plot lineages to tree
read_depth_HR_meta_plot_lineage <- read_depth_HR_meta[,c('Sample', 'lineage_7P')]
rownames(read_depth_HR_meta_plot_lineage) <- read_depth_HR_meta_plot_lineage$Sample
read_depth_HR_meta_plot_lineage$Sample <- NULL
# choose anough colors
lineage_palletter <- colorRampPalette(brewer.pal(8,"Accent") )(length(unique(read_depth_HR_meta_plot_lineage$lineage_7P[!is.na(read_depth_HR_meta_plot_lineage$lineage_7P)])))
#tree_heatmap_1b <- gheatmap(tree_heatmap_1b, read_depth_HR_meta_plot_lineage, offset=0.050, width = 0.15,
#                            colnames=T, legend_title="") + scale_fill_manual(values =lineage_palletter, na.value = "grey91")

# leave out 7PET as column (not included in branch color)
tree_heatmap_1b <- gheatmap(p, read_depth_HR_meta_plot_lineage, offset=0.00, width = 0.15,
                            colnames=T, legend_title="") + scale_fill_manual(values =lineage_palletter, na.value = "grey91")

tree_heatmap_1c <- tree_heatmap_1b + new_scale_fill()

# add whether the strain was isolated from Bangladesh or from another country
read_depth_HR_meta_plot_country <- read_depth_HR_meta[,c('Sample', 'Country.clean')]
read_depth_HR_meta_plot_country$Country.clean <- ifelse(read_depth_HR_meta_plot_country$Country.clean == 'Bangladesh', 'Bangladesh', 'Other')
rownames(read_depth_HR_meta_plot_country) <- read_depth_HR_meta_plot_country$Sample
read_depth_HR_meta_plot_country$Sample <- NULL
tree_heatmap_1cb <- gheatmap(tree_heatmap_1c, read_depth_HR_meta_plot_country, offset=0.050, width = 0.15,
                             colnames=T, legend_title="") + scale_fill_manual(values = c('navy', 'orange')) 
# add new fill scale
tree_heatmap_1d <- tree_heatmap_1cb + new_scale_fill()
# add number of HR
read_depth_HR_meta_plot_HR <- read_depth_HR_meta[,c('Sample', 'HR1_ratio_round', 'HR2_ratio_round')]
rownames(read_depth_HR_meta_plot_HR) <- read_depth_HR_meta_plot_HR$Sample
read_depth_HR_meta_plot_HR$Sample <- NULL
read_depth_HR_meta_plot_HR$HR1_ratio_round <- factor(read_depth_HR_meta_plot_HR$HR1_ratio_round, levels = c('0', '1', '2', '>2'))
read_depth_HR_meta_plot_HR$HR2_ratio_round <- factor(read_depth_HR_meta_plot_HR$HR2_ratio_round, levels = c('0', '1', '2', '>2'))

# although the factor level is defined, this somehow puts '>2' as lowest value, therefore, define colors manually
tree_heatmap_1d <- gheatmap(tree_heatmap_1d, read_depth_HR_meta_plot_HR, offset=0.10, width = 0.30,
                            colnames=T, legend_title="") + scale_fill_manual(values = c("#99000D","#FFF5F0","#FCAD91", "#F34A35")) 

tree_heatmap_1d

# --> check in blast out of the genomic islands if HR2 is part of any
# --> is almost entirely covering VSP-I

pdf('output_figures/public_short_read_tree_heatmap.pdf', width = 6, height =6)
tree_heatmap_1d
dev.off()

