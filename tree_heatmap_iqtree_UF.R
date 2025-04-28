# load packages
library('tidyr')
library('dplyr')
library('phytools')
library('ggplot2')
library('ggtree')
library('ggnewscale')
library('RColorBrewer')
library('colorspace')
library('treeio')
library('patchwork')
library(scales)

# setwd
setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# import basic metadata 
meta_basic <- read.csv2("input_data/very_basic_meta_2.csv", sep = ',')

# I had previsouly constructed a tree with RaXML which I realised was not the ideal appraoch for our dataset, as it assigns a minimal branch length even for the same sequences
#tree <- read.tree("input_data/RAxML_bestTree.raxmltree_all_ref_N16961")
# I therefore repeated this by using IQtree and ultrafast bootstrap approximation as input
tree <- read.tree('input_data/all_concat_highQ_variants_no_indels_concat_snp-sites.fasta.contree')

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

# color tree by ultrafast (UF) bootstrap approximation support. I am aware that these should pimarily be used on single gene trees, rather than on concatenated sequences as this one, so I will not interprete these quantitatively
tree.gg_1$data['UF_support'] <- as.numeric(as.character(ifelse(tree.gg_1$data$isTip == TRUE, '', tree.gg_1$data$label)))
tree.gg_111 <- tree.gg_1 + geom_tree(aes(color=UF_support)) + scale_colour_gradient(low='red', high='forestgreen')

# add tiplabels which are fused
# import info on how many chromosomes were assembled
QC <- read.csv2('input_data/QC_extended.csv')
QC <- QC[,c('Sample', 'chr')]
# output chromosome fusion state as binary endpoint for pyseer (GWAS)
pyseer <- QC
pyseer$chr <- ifelse(QC$chr %in% c("Two circular chromosomes assembled",   "Not all chromosomes circularised"), "0", "1")
write.table(pyseer, '../../01_data/16_check_singular_contig/pyseer/pyseer_fusion.tab', sep='\t', row.names = F, col.names = F, quote = F)

# for 8/9 where not all chromosomes were circularised, one chromosome was circularised, so count them to non-fused
QC$chr <- ifelse(QC$chr %in% c("Two circular chromosomes assembled",   "Not all chromosomes circularised"), "Two chromosomes assembled", "One chromosome assembled")
# order such that one 'two chr' is plotted fist
QC$chr <- factor(QC$chr, levels = c("Two chromosomes assembled", "One chromosome assembled"))
QC_keep <- QC

# add new color scale
tree.gg_111b <- tree.gg_111 + new_scale_color()
# add tippoints
tree.gg_2 <- tree.gg_111b %<+% QC + geom_tippoint(aes(color=as.factor(chr)), size=1) + scale_color_manual(values = c('darkblue', 'lightblue'), na.value = "grey91")
# merge
QC_meta <- merge(colonies, QC_keep, by.x = 'colony', by.y='Sample')
# find all which are in household with at least one fused
# this takes too long, do this household by household
fusion_hh <- unique(QC_meta[QC_meta$chr == 'One chromosome assembled',]$household)

## for each household where at least one fused chromosome isolate was detected, I compared all against all using clair3. therefore, I here printed all fusion-household isolate names
#QC_meta[QC_meta$household %in% fusion_hh[2],]$colony
#QC_meta[QC_meta$household %in% fusion_hh[1],]$colony
#QC_meta[QC_meta$household %in% fusion_hh[3],]$colony
#QC_meta[QC_meta$household %in% fusion_hh[4],]$colony
#QC_meta[QC_meta$household %in% fusion_hh[5],]$colony

# plot proportion of fused and non-fused strains for each of these households. 
QC_meta_fusion_households <- QC_meta[QC_meta$household %in% fusion_hh,]
#order according to tree
QC_meta_fusion_households$household <- factor(QC_meta_fusion_households$household, levels = c('E', 'P', 'I', 'F', 'H'))

plots <- list()
for (i in 1:length(levels(QC_meta_fusion_households$household))){
  plots[[i]] <- ggplot(QC_meta_fusion_households[QC_meta_fusion_households$household == levels(QC_meta_fusion_households$household)[i],]) + 
    geom_bar(aes(x=Individual, fill=chr)) + ylim(0,10) + scale_fill_manual(values = c('darkblue', 'lightblue'), na.value = "grey91") + theme_light() + theme(legend.position = 'none')  #+
    #coord_fixed()
}


# export
pdf('output_figures/proportion_fusion_household.pdf', height=8, width=2)
plots[[1]] / plots[[2]] / plots[[3]]  / plots[[4]] / plots[[5]]  
dev.off()

# highlight the branches with a fused chromoosme and label the nodes
d <- data.frame(node=c(869, 615, 534,  498, 926), chr_high=c("fusion_branch"))

tree.gg_2a <- tree.gg_2  + geom_hilight(data=d, aes(node=node, fill=chr_high), align="right",
                          extend=0.0035, #type = "gradient", gradient.direction = 'rt',
                          alpha = .8) + scale_fill_manual(values='lightblue') + 
  geom_point2(aes(subset=(node %in% c(869, 615, 534,  498, 926))), shape=23, size=3, fill='red') # add nodes to label
tree.gg_2a

# add new fill scale
tree.gg_2b <- tree.gg_2a + new_scale_fill()

# add first column (household)
tree_heatmap <- gheatmap(tree.gg_2b, plot, width = 0.08,
                         colnames=T, legend_title="") + scale_fill_manual(values = c(getPalette(colourCount), "#9ECAE1",  "#DEEBF7"), na.value = "grey91") + geom_rootedge()
# add new color (fill) scale, such that a new column can be added using a different color scheme
tree_heatmap_1 <- tree_heatmap + new_scale_fill()

# add lineage info
lineage <- read.table('input_data/lineage_assignemnet_from_medaka.txt', header = T)
lineage_plot <- lineage[,c('sample', 'Lineage_group')]
# samples to rownames, such that in gheatmap it can be added to tree
rownames(lineage_plot) <- lineage_plot$sample
lineage_plot$sample <- NULL
lineage_plot$Lineage_group <- factor(lineage_plot$Lineage_group, levels = c('BD2', 'IND1.2', 'IND1.3'))

tree_heatmap_1 <- gheatmap(tree_heatmap_1, lineage_plot, width = 0.08, offset = 0.0008,
                         colnames=T, legend_title="") + scale_fill_manual(values = c("#FDD58C", "#4763AB", "#B2258F"), na.value = "grey91")
# add new color (fill) scale, such that a new column can be added using a different color scheme
tree_heatmap_1 <- tree_heatmap_1 + new_scale_fill()
tree_heatmap_1

rownames(QC) <- QC$Sample
QC$Sample <- NULL
#QC$tick <- NULL
tree_heatmap_1c <- gheatmap(tree_heatmap_1, QC, offset=0.0016, width = 0.08,
                            colnames=T, legend_title="") + scale_fill_manual(values = c('lightblue', 'darkblue'), na.value = "grey91")
# add new fill scale
tree_heatmap_1cb <- tree_heatmap_1c + new_scale_fill()

# import how many copies of HS1 (aka 'overlap') have been detected per genome
overlap_sum <- read.table('input_data/blastout_HR_all_sum.tab', header = T)
rownames(overlap_sum) <- overlap_sum$Sample
overlap_sum$Sample <- NULL
overlap_sum$on_x_chr <- NULL
overlap_sum$n_HR <- factor(overlap_sum$n_HR)
tree_heatmap_1cc <- gheatmap(tree_heatmap_1cb, overlap_sum, offset=0.0024, width = 0.08,
                            colnames=T, legend_title="") + scale_fill_manual(values = c("#FCAD91", "#F34A35"))
tree_heatmap_1cc


pdf("output_figures/all_tree_heatmap_chr_HR_v3_iqtree_highlight_UF.pdf", height = 10, width = 8)
tree_heatmap_1cc
dev.off()


### this was done on he iqtree with no bootsrapps so the nodelabels dont match
## add chr info to meta
#meta_chr <- merge(QC, colonies, by=0)
## Find all in same household with fused, but which are unfused
#meta_chr[meta_chr$household %in% unique(meta_chr[meta_chr$chr == 'One chromosome assembled',]$household) & meta_chr$chr != 'One chromosome assembled',]$colony
#
## import called variants for 114Vc03 and 117Vc01 (have single chr and are being tested by Ana)
#variants_114Vc03 <- read.table('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/01_data/20_clair3/114Vc03_merge_output.vcf')
#variants_114Vc03['query'] <- '114Vc03'
#variants_117Vc01 <- read.table('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/01_data/20_clair3/117Vc01_merge_output.vcf')
#variants_117Vc01['query'] <- '117Vc01'
#
#variants <- merge(variants_114Vc03, variants_117Vc01, all = T)
#colnames(variants) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE', 'query')
#variants['ref'] <- gsub('\\/.*', '', variants$CHROM)
#variants$CHROM <- gsub('.*\\:', '', variants$CHROM)
#variants['Allele_Frequency_ALT'] <- gsub('.*\\:', '', variants$SAMPLE)
#variants['Estimated_read_depth'] <- gsub('(\\d+\\/\\d+\\:)(\\d+\\:)(\\d+)(\\:\\d+\\,\\d+)(\\:.*)', '\\3', variants$SAMPLE)
#
#variants_sum <- variants %>% group_by(query, ref) %>%
#  summarise(n_hq_variants = sum(Estimated_read_depth >=10 & Allele_Frequency_ALT >=0.8 & QUAL >=20))
#variants_sum[variants_sum$query == '114Vc03' & variants_sum$n_hq_variants == min(variants_sum[variants_sum$query == '114Vc03', ]$n_hq_variants),]
## for 114Vc03, either choose 114Vc01 or 114Vc02
#variants_sum[variants_sum$query == '117Vc01' & variants_sum$n_hq_variants == min(variants_sum[variants_sum$query == '117Vc01', ]$n_hq_variants),]
## for 117Vc01, either choose 117Vc07 or 117Vc10
#
## plot subtrees for each of the five branches with a fused strain
#tree.gg_2$data['label_nodes'] <- ifelse(tree.gg_2$data$isTip == T, NA, tree.gg_2$data$node)
#pdf('/Users/alinecuenod/Downloads/check_nodelabs.pdf', height=45,width = 25)
#tree.gg_2 + geom_text(aes(label=label_nodes)) + geom_tippoint(aes(color=chr))
#dev.off()
#
#meta_chr['ind_contact_label'] <- ifelse(meta_chr$ind_contact == 'index', 'i', 'c')
##meta_chr['ind_contact_label'] <- '-'
#
## nodes of interest to subset are 747, 887, 543, 491, 717
#tree_747 <- tree_subset(tree_b, 747, levels_back = 0)
#ggtree_747 <- ggtree(tree_747, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_747 +  geom_tippoint(aes(color=chr)) + geom_text(aes(label=node))
## try 767
#tree_767 <- tree_subset(tree_b, 767, levels_back = 0)
#ggtree_767 <- ggtree(tree_767, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_767 <- ggtree_767 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), align=T, linetype = 'solid') + scale_color_manual(values = c('darkblue', 'lightblue',  '#f47e71'), na.value = "grey91")+ # only index cases
#  theme(legend.position = 'none')+ geom_treescale() #7
#
## 887
#tree_887 <- tree_subset(tree_b, 887, levels_back = 0)
#ggtree_887 <- ggtree(tree_887, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_887 +  geom_tippoint(aes(color=chr)) 
## 887
#tree_895 <- tree_subset(tree_b, 895, levels_back = 0)
#ggtree_895 <- ggtree(tree_895, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_895 <- ggtree_895 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), align=T, linetype = 'solid') + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+
#  theme(legend.position = 'none') + geom_treescale() #1 at 897
#
## 543
#tree_543 <- tree_subset(tree_b, 543, levels_back = 0)
#ggtree_543 <- ggtree(tree_543, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_543 +  geom_tippoint(aes(color=chr)) 
## 584
#tree_584 <- tree_subset(tree_b, 584, levels_back = 0)
#ggtree_584 <- ggtree(tree_584, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_584 <- ggtree_584 + geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), align=T, linetype = 'solid') + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+
#  theme(legend.position = 'none') + geom_treescale()#6
#tree_584
#
## 491
#tree_491 <- tree_subset(tree_b, 491, levels_back = 0)
#ggtree_491 <- ggtree(tree_491, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_491 +  geom_tippoint(aes(color=chr)) 
## 512
#tree_512 <- tree_subset(tree_b, 512, levels_back = 0)
#ggtree_512 <- ggtree(tree_512, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_512 <- ggtree_512 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), align=T, linetype = 'solid') + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+
#  theme(legend.position = 'none')+ geom_treescale()#21
#
## 717
#tree_717 <- tree_subset(tree_b, 717, levels_back = 0) 
#ggtree_717 <- ggtree(tree_717, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_717 <- ggtree_717 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), align=T, linetype = 'solid') + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+ 
#  theme(legend.position = 'none') + geom_treescale()
#tree_717
## 23
#
#
## nodes to plot 717, 512, 584. 887, 767
## plot together 
#tree_895 + tree_767 + tree_584 + tree_512 + tree_717 +
#  plot_layout(ncol = 1)
#
#pdf("output_figures/subtrees_fusion_lines.pdf", height = 8, width = 1.1)
#tree_895 + tree_767 + tree_584 + tree_512 + tree_717 + # plotted in same order as occuring in the tree (top to bottom)
#  plot_layout(ncol = 1)
#dev.off()
#
## try offset
#
## nodes of interest to subset are 747, 887, 543, 491, 717
#tree_747 <- tree_subset(tree_b, 747, levels_back = 0)
#ggtree_747 <- ggtree(tree_747, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_747 +  geom_tippoint(aes(color=chr)) + geom_text(aes(label=node))
## try 767
#tree_767 <- tree_subset(tree_b, 767, levels_back = 0)
#ggtree_767 <- ggtree(tree_767, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_767 <- ggtree_767 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), offset=4e-06) + scale_color_manual(values = c('darkblue', 'lightblue',  '#f47e71'), na.value = "grey91")+ # only index cases
#  theme(legend.position = 'none')+ geom_treescale() #7
#tree_767
#
## 887
#tree_887 <- tree_subset(tree_b, 887, levels_back = 0)
#ggtree_887 <- ggtree(tree_887, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_887 +  geom_tippoint(aes(color=chr)) 
## 887
#tree_895 <- tree_subset(tree_b, 895, levels_back = 0)
#ggtree_895 <- ggtree(tree_895, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_895 <- ggtree_895 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), offset=4e-08) + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+
#  theme(legend.position = 'none') + geom_treescale() #1 at 897
#tree_895
## 543
#tree_543 <- tree_subset(tree_b, 543, levels_back = 0)
#ggtree_543 <- ggtree(tree_543, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_543 +  geom_tippoint(aes(color=chr)) 
## 584
#tree_584 <- tree_subset(tree_b, 584, levels_back = 0)
#ggtree_584 <- ggtree(tree_584, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_584 <- ggtree_584 + geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), offset=4e-06) + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+
#  theme(legend.position = 'none') + geom_treescale()#6
#tree_584
#
## 491
#tree_491 <- tree_subset(tree_b, 491, levels_back = 0)
#ggtree_491 <- ggtree(tree_491, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#ggtree_491 +  geom_tippoint(aes(color=chr)) 
## 512
#tree_512 <- tree_subset(tree_b, 512, levels_back = 0)
#ggtree_512 <- ggtree(tree_512, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_512 <- ggtree_512 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), offset=4e-06) + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+
#  theme(legend.position = 'none')+ geom_treescale()#21
#
## 717
#tree_717 <- tree_subset(tree_b, 717, levels_back = 0) 
#ggtree_717 <- ggtree(tree_717, layout = 'rectangular', ladderize = FALSE)  %<+% meta_chr
#tree_717 <- ggtree_717 +  geom_tippoint(aes(color=chr)) + geom_tiplab(aes(label=ind_contact_label, col=ind_contact_label), offset=4e-06) + scale_color_manual(values = c('darkblue', 'lightblue', '#beb9d8', '#f47e71'), na.value = "grey91")+ 
#  theme(legend.position = 'none') + geom_treescale()
#tree_717
## 23
#
#
## nodes to plot 717, 512, 584. 887, 767
## plot together 
#tree_895 + tree_767 + tree_584 + tree_512 + tree_717 +
#  plot_layout(ncol = 1)
#
#pdf("output_figures/subtrees_fusion_offset.pdf", height = 8, width = 1.1)
#tree_895 + tree_767 + tree_584 + tree_512 + tree_717 + # plotted in same order as occuring in the tree (top to bottom)
#  plot_layout(ncol = 1)
#dev.off()
#
#


# eval clair 3 within fusion housholds
# all
#all_household <- read.table('input_data/clair3_fusion_households_merge_output_sum_household_all.tab')

all_household <- read.table('input_data/clair3_fusion_households_merge_output_sum_fusion_households_non_vc_assemblies_removed.tab')
colnames(all_household) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE', 'query_ref', 'household')
all_household['query'] <- gsub('(\\d{3}Vc\\d{2})(\\d{3}Vc\\d{2})', '\\1', all_household$query_ref)
all_household['ref'] <- gsub('(\\d{3}Vc\\d{2})(\\d{3}Vc\\d{2})', '\\2', all_household$query_ref)
# only consider SNV > 20 QUAL
all_household$QUAL <- as.numeric(as.character(all_household$QUAL))
all_household_hq <- all_household[all_household$QUAL >= 20,]
# check if in any of the hq SNV more than one allele called
grep('\\,', all_household_hq$ALT) #--> no
# add allele freq and read depth
all_household_hq['Allele_Frequency_ALT'] <- gsub('.*\\:', '', all_household_hq$SAMPLE)
all_household_hq['Estimated_read_depth'] <- gsub('(\\d+\\/\\d+\\:)(\\d+\\:)(\\d+)(\\:\\d+\\,\\d+)(\\:.*)', '\\3', all_household_hq$SAMPLE)
# only consider >= 0.8 allele freq and >= 10 read depth
all_household_hq$Allele_Frequency_ALT <- as.numeric(as.character(all_household_hq$Allele_Frequency_ALT))
all_household_hq$Estimated_read_depth <- as.numeric(as.character(all_household_hq$Estimated_read_depth))
all_household_hq <- all_household_hq[all_household_hq$Allele_Frequency_ALT >= 0.8 & all_household_hq$Estimated_read_depth >= 10,]

#all_household_hq['check'] <- paste0(all_household_hq$CHROM, all_household_hq$POS, all_household_hq$ID, all_household_hq$REF, all_household_hq$ALT, all_household_hq$query_ref,   all_household_hq$household,   all_household_hq$query,     all_household_hq$ref)
#all_household_hq_old['check'] <- paste0(all_household_hq_old$CHROM, all_household_hq_old$POS, all_household_hq_old$ID, all_household_hq_old$REF, all_household_hq_old$ALT, all_household_hq_old$query_ref,   all_household_hq_old$household,   all_household_hq_old$query,     all_household_hq_old$ref)

# summarise for each household
E_household_hq_sum <- all_household_hq[all_household_hq$household == 'household_E',] %>% group_by(query, ref) %>% summarise(n_hq_SNV = n()) %>% ungroup %>% complete(query, ref, fill=list(n_hq_SNV= 0))
# color label my fused and unfused
x_cols <- ifelse(unique(E_household_hq_sum$query) %in% QC_meta[QC_meta$chr == 'One chromosome assembled',]$colony, 'lightblue', 'darkblue')
# plot
household_E_plot <-ggplot(E_household_hq_sum, aes(query, ref, fill= n_hq_SNV)) + 
  geom_tile() +
  geom_text(aes(label=n_hq_SNV), col='white') +
  #scale_fill_gradientn(
  #  colours = c("blue", "orange", "darkred"),
  #  values = c(0, 5, 29)) + 
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 10, l3 = 0, p3 = .8, p4 = .6, limits = c(0, 29)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=x_cols), axis.text.y = element_text(colour = x_cols), legend.position = 'none') + ggtitle('Household E (Event 1)')
household_E_plot

# household F
F_household_hq_sum <- all_household_hq[all_household_hq$household == 'household_F',] %>% group_by(query, ref) %>% summarise(n_hq_SNV = n()) %>% ungroup %>% complete(query, ref, fill=list(n_hq_SNV= 0))
# color label my fused and unfused
x_cols <- ifelse(unique(F_household_hq_sum$query) %in% QC_meta[QC_meta$chr == 'One chromosome assembled',]$colony, 'lightblue', 'darkblue')
# plot
household_F_plot <-ggplot(F_household_hq_sum, aes(query, ref, fill= n_hq_SNV)) + 
  geom_tile() +
  geom_text(aes(label=n_hq_SNV), col='white') +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 10, l3 = 0, p3 = .8, p4 = .6, limits = c(0, 29)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=x_cols), axis.text.y = element_text(colour = x_cols), legend.position = 'none')  + ggtitle('Household F (Event 5)')
range(F_household_hq_sum$n_hq_SNV)


# household H
H_household_hq_sum <- all_household_hq[all_household_hq$household == 'household_H',] %>% group_by(query, ref) %>% summarise(n_hq_SNV = n()) %>% ungroup %>% complete(query, ref, fill=list(n_hq_SNV= 0))
# color label my fused and unfused
x_cols <- ifelse(unique(H_household_hq_sum$query) %in% QC_meta[QC_meta$chr == 'One chromosome assembled',]$colony, 'lightblue', 'darkblue')
# plot
household_H_plot <-ggplot(H_household_hq_sum, aes(query, ref, fill= n_hq_SNV)) + 
  geom_tile() +
  geom_text(aes(label=n_hq_SNV), col='white') +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 10, l3 = 0, p3 = .8, p4 = .6, limits = c(0, 29)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=x_cols), axis.text.y = element_text(colour = x_cols), legend.position = 'none')  + ggtitle('Household H (Event 2)')
range(H_household_hq_sum[grepl('118',H_household_hq_sum$ref) & grepl('117',H_household_hq_sum$query),]$n_hq_SNV)
range(H_household_hq_sum[grepl('117',H_household_hq_sum$ref) & grepl('118',H_household_hq_sum$query),]$n_hq_SNV)

# household I
I_household_hq_sum <- all_household_hq[all_household_hq$household == 'household_I',] %>% group_by(query, ref) %>% summarise(n_hq_SNV = n()) %>% ungroup %>% complete(query, ref, fill=list(n_hq_SNV= 0))
# color label my fused and unfused
x_cols <- ifelse(unique(I_household_hq_sum$query) %in% QC_meta[QC_meta$chr == 'One chromosome assembled',]$colony, 'lightblue', 'darkblue')
# plot
household_I_plot <-ggplot(I_household_hq_sum, aes(query, ref, fill= n_hq_SNV)) + 
  geom_tile() +
  geom_text(aes(label=n_hq_SNV), col='white') +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 10, l3 = 0, p3 = .8, p4 = .6, limits = c(0, 29)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=x_cols), axis.text.y = element_text(colour = x_cols), legend.position = 'none')  + ggtitle('Household I (Event 4)')
range(I_household_hq_sum$n_hq_SNV)



# household P
P_household_hq_sum <- all_household_hq[all_household_hq$household == 'household_P',] %>% group_by(query, ref) %>% summarise(n_hq_SNV = n()) %>% ungroup %>% complete(query, ref, fill=list(n_hq_SNV= 0))
# color label my fused and unfused
x_cols <- ifelse(unique(P_household_hq_sum$query) %in% QC_meta[QC_meta$chr == 'One chromosome assembled',]$colony, 'lightblue', 'darkblue')
# plot
household_P_plot <-ggplot(P_household_hq_sum, aes(query, ref, fill= n_hq_SNV)) + 
  geom_tile() +
  geom_text(aes(label=n_hq_SNV), col='white') +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 10, l3 = 0, p3 = .8, p4 = .6, limits = c(0, 29)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=x_cols), axis.text.y = element_text(colour = x_cols), legend.position = 'right')  + ggtitle('Household P (Event 3)')
household_P_plot
range(P_household_hq_sum$n_hq_SNV)

pdf('output_figures/clair3_heatmap_within_fusion_household_non_Vc_removed.pdf', height=10, width = 15)
household_E_plot + household_H_plot + household_P_plot + household_I_plot+ household_F_plot
dev.off()
