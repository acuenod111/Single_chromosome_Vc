# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# did all again with new assemblies
QC <- read.csv('input_data/quality_after_non-Vc_contig_removed.csv')

# remove last row (no sample)
QC <- QC[-nrow(QC),]

# import flye assembly stats
flye_info <- read.csv2('input_data/all_fly_assembly_info.txt', sep = '\t')
# add column names
colnames(flye_info) <- c("seq_name","length","cov.","circ.","repeat.","mult.","alt_group","graph_path","sample")

# I ran kraken on the assemblies to identify the ones where there are non-Vc contigs (most probably stemming from sequencing artefacts)
Vc_kraken <- read.table('input_data/Vc_all_report_kraken2.txt',sep='\t')
Vc_kraken['sample'] <- gsub('\\_.*', '', Vc_kraken$V1)
#check if all only once
unique(table(Vc_kraken$sample)) # all exactly once
# add perc Vc
Vc_kraken['perc_Vc'] <- gsub('.*\\:', '', Vc_kraken$V1)
# for the ones which were not 100 Vc I ran SprayNPray to identify the non-Vc contigs
sprayNpray <- read.csv('input_data/all_sprayNpray_out.csv')
sprayNpray <- sprayNpray[!sprayNpray$contig == 'contig',]
# check length of all non-Vc contigs
sprayNpray$contig_length <- as.numeric(as.character(sprayNpray$contig_length))
quantile(sprayNpray[!grepl('Vibrio cholerae',sprayNpray$closest_blast_hits),]$contig_length)
#View(sprayNpray[!grepl('Vibrio cholerae',sprayNpray$closest_blast_hits),])

non_Vc_contigs <- sprayNpray[!grepl('Vibrio cholerae',sprayNpray$closest_blast_hits),]$contig # I manually check, these are all repetitive and probably sequencing artefact
#write.table(sprayNpray[!grepl('Vibrio cholerae',sprayNpray$closest_blast_hits),]$contig, row.names = F, quote = F, '/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/01_data/02_ONT_QC/non-Vc_contigs.tab')


# check circularity of contigs from the ReCycled output
# more chromosomes were identified as circularly assembled using ReCycled, which checks for reads spanning the end-start of each contig 
recycled <- read.table('input_data/all_ReCycled_logs.tab', comment.char = "", header = T)
# some contigs were identified as non-Vc after running recycled, remove these
recycled['contig_name'] <- gsub('(.*\\:)(\\d{3}Vc\\d{3,})(\\_.*$)', '\\2', recycled$X.contigName)
recycled <- recycled[!recycled$contig_name %in% non_Vc_contigs, ]

# it only sais 'Circularizable' for chromosome 1, but thats because it does not find dnaa on chromosome 2, so I just check if i find overlapping reads
table(recycled$StartOfTargetMapping)
range(recycled[recycled$OverlapCov == 0,]$ContigLength) # some small contigs are not circularly assembled
range(recycled[recycled$ContigLength > 800000,]$OverlapCov) # for all chromosomes, I have reads spanning the break --> circularly assembled
recycled['overlap_dev_total_cov'] <- recycled$OverlapCov / recycled$ContigCov
#check overlap fraction of the ones I consider circular
quantile(recycled[recycled$OverlapCov >0,]$overlap_dev_total_cov) # looks ok


recycled$X.contigName <- gsub('^.*\\:', '', recycled$X.contigName)
recycled['Isolate'] <- gsub('(\\d{3}Vc\\d{2})(.*)', '\\1', recycled$X.contigName)
recycled['chromosome'] <- ifelse(recycled$ContigLength > 2000000, 1, 
                                 ifelse(recycled$ContigLength <= 2000000 & recycled$ContigLength > 800000, 2, 
                                        ifelse(recycled$ContigLength < 800000, 'extrachromosomal element', NA)))
table(recycled$chromosome)

recycled['topology'] <- ifelse(recycled$OverlapCov > 0, 'circular', 'unknown')
recycled['completeness'] <- ifelse(recycled$OverlapCov > 0, 'complete', 'partial')

table(recycled$chromosome, recycled$topology) 
# all chromosomes were circularly assembled (409xtwo chromosomes, 58 one large fused chromosome). 
# further, we assemble 108 extrachromosomal contigs (103 circular and 5 non-circular ones)

unique(recycled[recycled$chromosome == 'extrachromosomal element',]$Isolate)
length(unique(recycled[recycled$chromosome == 'extrachromosomal element',]$Isolate))

# create table to rename contigs for NCBI submission
contig_rename <- data.frame(old=recycled$X.contigName, new = paste0('[organism=Vibrio cholerae][isolate=', recycled$Isolate,'][chromosome=', recycled$chromosome,'][topology=', recycled$topology, '][completeness=', recycled$completeness, ']'))
contig_rename$new <- gsub('chromosome=extrachromosomal contig', 'extrachromosomal contig', contig_rename$new)
#write.table(contig_rename, '/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/01_data/02_ONT_QC/rename_contigs.tab', quote = F, row.names = F, col.names = F)

recycled_sum <- recycled %>% group_by(Isolate) %>% summarise(n_chr = sum(ContigLength > 800000), 
                                                                fract_circ = mean(OverlapCov[ContigLength> 800000] >0))

# add column indicating whether (i) not all chromosomes where fully circularised (n=9) and where one fused chromosome was circularry assembled
#QC['chr'] <- ifelse(QC$Sample %in% flye_info_sum[flye_info_sum$n_chr == '1', ]$sample, 'One chromosome assembled', 
#                    ifelse(QC$Sample %in% not_all_chr_circ$sample, 'Not all chromosomes circularised', 'Two circular chromosomes assembled'))

QC['chr'] <- ifelse(QC$Sample %in% recycled_sum[recycled_sum$n_chr == '1' & recycled_sum$fract_circ == 1, ]$Isolate, 'One chromosome assembled', 
                    ifelse(QC$Sample %in% recycled_sum[recycled_sum$n_chr == '2' & recycled_sum$fract_circ == 1,]$Isolate, 'Two circular chromosomes assembled', NA))

# write 
write.table(QC, 'input_data/QC_extended.csv', row.names = F, quote = F, sep = ';')

# add meta
meta <- read.csv2('input_data/very_basic_meta_2.csv', sep=',')
QC['Individual'] <- gsub('(\\d{3})(Vc\\d{2})', '\\1', QC$Sample)
QC_meta <- merge(QC, meta, by='Individual')
write.table(QC_meta, '../../01_data/02_ONT_QC/QC_meta.csv', row.names = F, quote = F, sep = ';')


# plot some basic QC measures
# plot the N50  of the reads
hist_head_length <- ggplot(QC, aes(x=N50_reads)) + 
  geom_histogram(bins = 200) + xlab('N50 (reads)') + ylab('Number of samples') + theme_light() 

# quality score vs. coverage
qual_depth_mean <- ggplot(QC, aes(x=Read_depth, y=Mean_read_quality_nanostat)) + geom_point() +theme_light() +ylab('Mean Read quality') + xlab('Average read depth') + scale_color_manual(values = c('orange', 'blue','grey'))

# quality score vs. read depth
qual_depth <- ggplot(QC, aes(x=Read_depth, y=Median_read_quality_nanostat)) + geom_point(alpha=0.7)+ theme_light()+theme(legend.position = 'None') +ylab('Median read quality') + xlab('Average read depth') + scale_color_manual(values = c('orange', 'blue','grey'))

# assembly length vs. read depth
length_depth <- ggplot(QC, aes(y=Total_length_assembly, x=Read_depth, col = chr)) +geom_point(alpha=0.7) +theme_light()+theme(legend.position = 'None')  +xlab('Average read depth') + ylab('Assembly length')+ scale_color_manual(values = c('orange', 'blue','grey'))
length_depth <- ggplot(QC, aes(y=Total_length_assembly, x=Read_depth)) +geom_point(alpha=0.7) +theme_light()+theme(legend.position = 'None')  +xlab('Average read depth') + ylab('Assembly length')+ scale_color_manual(values = c('orange', 'blue','grey'))

# print stats on QC measures
quantile(QC$Total_length, na.rm = T)
quantile(QC$GC_percent, na.rm = T)
quantile(QC$Read_depth, na.rm = T)
quantile(QC$Mean_read_quality_nanostat, na.rm = T)
quantile(QC$Total_length_assembly, na.rm = T)
quantile(QC$N50_reads, na.rm = T)

QC_p <- cowplot::plot_grid(hist_head_length, qual_depth, length_depth, nrow = 1)

pdf('output_figures/ONT_QC_all_no_colors.pdf', width = 9, height =3.5)
QC_p
dev.off()


# import basic metadata
meta_basic <- read.csv2("input_data/very_basic_meta_2.csv", sep = ',')

# merge to QC
QC['Individual'] <- gsub('Vc.*', '', QC$Sample)
QC_household <- merge(meta_basic, QC[,c('Individual', 'Sample', 'chr')])
QC_household['chr_num'] <- ifelse(QC_household$chr == 'One chromosome assembled', 1, 2)

# add in which households a fusion event was detected
QC_household_sum <- QC_household %>% group_by(household) %>%
  reframe(mean_chr_nr = mean(chr_num), 
          chr_state_sum = ifelse(mean(chr_num) == 2, 'No fusion in household', 'Fusion detected'),# %>%
          chr_fusion_ind_contact = toString(unique(ind_contact[chr_num == 1])))
QC_household_sum['transmission'] <- ifelse(QC_household_sum$chr_fusion_ind_contact == 'index, contact', 'FusionInIndexAndContact', 
                                           ifelse(QC_household_sum$chr_fusion_ind_contact == 'contact', 'FusionInContactOnly',
                                                  ifelse(QC_household_sum$chr_fusion_ind_contact == 'index', 'FusionInIndexOnly', 'NoFusionDetected')))

QC_household_sum$transmission <- factor(QC_household_sum$transmission, levels = c('NoFusionDetected', 'FusionInIndexAndContact', 'FusionInContactOnly', 'FusionInIndexOnly'))
QC_household_sum$chr_state_sum <- factor(QC_household_sum$chr_state_sum, levels = c('No fusion in household', 'Fusion detected'))

fusion_sum_p <- ggplot(QC_household_sum, aes(x= chr_state_sum, fill = factor(transmission), col = factor(transmission))) +
  geom_bar(position="stack",stat="count")  + scale_fill_manual(values = c('darkblue', 'lightblue', 'lightblue', 'lightblue')) +
  scale_color_manual(values = c('darkgrey','darkgrey','darkgrey','darkgrey')) +
  ylab('Nr. of households') + xlab('') + theme_light()
fusion_sum_p

pdf('output_figures/fusion_sum_p_no_pattern.pdf', width = 3.0, height =4.5)
fusion_sum_p + theme(legend.position = 'bottom')
dev.off()

