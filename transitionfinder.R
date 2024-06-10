#  load packages
library('tidyr')
library('dplyr')
library('phytools')
library('ggplot2')
library('ggtree')
library('ggnewscale')
library('RColorBrewer')
library('colorspace')
library('ape')

setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

#load 'transitionfinder.R' script
source('../transitionsfinder.R') # this was downloaded from GitHUb (https://github.com/Saannah/QuidiPhydy)

# import basic metadata 
meta_basic <- read.csv2("input_data/very_basic_meta_1.csv", sep = ',')

# plot the tree similarly as for Fig 2C
# import the tree
tree <- read.tree("input_data/RAxML_bestTree.raxmltree_all_ref_N16961")
# rotate around two nodes (for visualisation)
tree_1 <- ape::rotate(tree, node = 643)
tree_2 <- ape::rotate(tree_1, node = 490)
# root on the reference
tree <- root(tree_2, outgroup = 'NZ_LT907989.1_Vibrio_cholerae_strain_A19')
# plot basic ggtree
tree.gg_1 <- ggtree(tree, layout = 'rectangular', ladderize = FALSE) + geom_treescale()


# Add info to tree for which assembly how many chromosomes were assembled
#chr_plot <- read.table('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/01_data/08_seqs_to_screen/outputs/merge_blast_plot.csv', sep=';', header = T)
QC <- read.csv2('input_data/QC_extended.csv')
# summarise the ones where two large contigs have been assembled as carrying two chromosomes (independent of whether these were asselmbled to circular chromosomes or not, as I do not have evidence for these to be fused)
QC$chr <- ifelse(QC$chr %in% c('Not all chromosomes circularised','Two circular chromosomes assembled'), 2, 1)

# add a line for the reference (has two chromosomes)
Ref<-data.frame('NZ_LT907989.1_Vibrio_cholerae_strain_A19','2')
names(Ref)<-c('Sample', 'chr')
QC <- merge(QC, Ref, by = c('Sample', 'chr'), all=T)

# supset to only include chr column
chr_plot <- QC[,c('Sample', 'chr')]
chr_plot$chr <- paste0(chr_plot$chr, '_chr')
rownames(chr_plot) <- chr_plot$strain
chr_plot$strain <- NULL

# build an array which includes the number of chromosomes for each strain in the same order as the tree tiplabels (needed for ancestral state reconstruction)
QC <- QC[match(tree$tip.label, QC$Sample),]
chr <- array(paste0(QC$chr, '_chr'))

# using transitionfinder find fusion events (transition from state 2 chr to 1 chr)
MRCAin_2to1chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("2_chr", "1_chr"))$MRCAin
MRCAout_2to1chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("2_chr", "1_chr"))$MRCAout
node_reconstruction_2to1chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("2_chr", "1_chr"))$reconstruction
all_transitions_2to1chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("2_chr", "1_chr"))$all_transitions

# using transitionfinder find fission (aka unfusion) events (transition from state 1 chr to 2 chr)
MRCAin_1to2chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("1_chr", "2_chr"))$MRCAin
MRCAout_1to2chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("1_chr", "2_chr"))$MRCAout
node_reconstruction_1to2chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("1_chr", "2_chr"))$reconstruction
all_transitions_1to2chr <- transitionsfinder(tree = tree, trait = chr, fromto = c("1_chr", "2_chr"))$all_transitions
# should be the same and are!
all(node_reconstruction_2to1chr == node_reconstruction_1to2chr)
node_reconstruction <- node_reconstruction_2to1chr

# plot chromosome state to tips
p <- ggtree(tree, layout = 'circular', branch.length="none")  %<+% chr_plot
p2 <- p + geom_tippoint(aes(col=factor(as.character(chr))), size=1, show.legend = T) + scale_color_manual(values = c('lightblue', 'darkblue'))

# label at which nodes transition happened as these whould be highlighted in the plot
node_reconstruction['transition'] <- ifelse(node_reconstruction$node %in% all_transitions_2to1chr, 'transition_2to1', 
                                            ifelse(node_reconstruction$node %in% all_transitions_1to2chr, 'transition_1to2', 'no transition'))
# add a column only including transitions
node_reconstruction['nodelab_check'] <- ifelse(node_reconstruction$transition %in% c('transition_2to1', 'transition_1to2'), node_reconstruction$node, '')

# this only adds transitions to nodes which are not tips, also mark transitions at tips (see below)
p2 %<+% node_reconstruction + 
  geom_nodepoint(aes(col=factor(labels), shape=factor(transition), size=factor(transition))) + scale_shape_manual(values = c(19,18,17)) + scale_size_manual(values = c(1,5,5)) +
  geom_text(aes(label=nodelab_check))

# add the informations whethe transition happened also for the tips
node_reconstruction_tipinfo <- merge(node_reconstruction, p$data, by = 'node', all = T)
# reorder the columns such that the 'label' column is first (chich matches the tiplabels)
node_reconstruction_tipinfo <- node_reconstruction_tipinfo[,c("label", colnames(node_reconstruction_tipinfo)[!colnames(node_reconstruction_tipinfo) == 'label'])]
# add this info to the tree
p0 <- p  %<+% node_reconstruction_tipinfo
# label chromosome state and where transition happened
p3 <- p0 + geom_tippoint(aes(col=factor(as.character(chr)), shape=factor(transition), size=factor(transition)), show.legend = T) + 
  scale_color_manual(values = c('lightblue', 'darkblue')) + scale_shape_manual(values = c(19,18,17)) + scale_size_manual(values = c(1,5,5)) + 
  geom_nodepoint(aes(col=factor(labels), shape=factor(transition), size=factor(transition))) + theme(legend.position = 'bottom')
p3 

pdf('output_figures/ASR_fusion_transitionfinder.pdf', width =7, height =7.7)
p3
dev.off()

# it looks like fusion happened at three relatively basal nodes
# for each node calculate the node to root distance
# extract all distances from tree
dist <- dist.nodes(tree)
# node 167 is the reference root tip (strain A19 and N16961 are the same)
ref_node <- node_reconstruction_tipinfo[grepl('A19',node_reconstruction_tipinfo$label),'node']

# extract distances for the different transition nodes (the rownumbers match the nodenumbers)
dist_trans <- as.data.frame(dist[c(all_transitions_1to2chr, all_transitions_2to1chr),ref_node])

dist_trans['trans'] <- ifelse(rownames(dist_trans) %in% all_transitions_1to2chr, 'Fission', 
                              ifelse(rownames(dist_trans) %in% all_transitions_2to1chr, 'Fusion', NA))
colnames(dist_trans)<- c('dist', 'trans')

#plot
comp_dist_p <- ggplot(dist_trans, aes(y=dist, x=trans)) +
   geom_jitter(width = 0.2, aes(col=trans, shape=trans), size=3) + scale_color_manual(values = c('darkblue', 'lightblue')) + scale_shape_manual(values = c(17,18)) +
  ylab('Distance from root [Substitotions per site]') + xlab('') + theme_light() + theme(legend.position = 'bottom')
# export
pdf('output_figures/comp_dist_transitionfinder_fusion.pdf', width = 3, height =3.5)
comp_dist_p  
dev.off()

# check the distance between fusion and the next downstream fusion event
dist_fusion_unfusion <- data.frame(dist[all_transitions_2to1chr, all_transitions_1to2chr])
dist_fusion_unfusion['fusion_node'] <- rownames(dist_fusion_unfusion)
dist_fusion_unfusion['min_dist'] <- apply(dist_fusion_unfusion, 1, FUN = min)
dist_fusion_unfusion['closest_unfusion_node'] = colnames(dist_fusion_unfusion)[apply(dist_fusion_unfusion, 1, which.min)]

# three relatively basal fusion nodes: 871, 901, 656
# find minimal dist
p3 +  geom_nodelab(aes(subset=!isTip & transition == 'transition_2to1',
                       label = nodelab_check))
fusion_nodes <- p3$data[p3$data$isTip == F & p3$data$transition == 'transition_2to1',]$nodelab_check

dist_fusion_unfusion_min_sel <- dist_fusion_unfusion[dist_fusion_unfusion$fusion_node %in% fusion_nodes, c('fusion_node', 'closest_unfusion_node', 'min_dist')]

# 4133751 is the mean assembly length
dist_fusion_unfusion_min_sel['Nr_SNV'] <- as.numeric(as.character(dist_fusion_unfusion_min_sel$min_dist)) * 4133751
range(dist_fusion_unfusion_min_sel['Nr_SNV'])
# In the third wave of the 7PET, the rates of SNP accumulation is estimated to 3.5 / year https://www.nature.com/articles/nature10392
dist_fusion_unfusion_min_sel['estimated_time'] <- dist_fusion_unfusion_min_sel$Nr_SNV / 3.5
range(dist_fusion_unfusion_min_sel['estimated_time'])

