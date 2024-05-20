# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# did all again with new assemblies
QC <- read.csv('input_data/quality.csv')

# remove 'sample' 1139Vc03 (wrong sample name in QC script, therefore no data)
QC <- QC[QC$Sample != '1139Vc03',]

# remove last row (no sample)
QC <- QC[-nrow(QC),]

# import flye assembly stats
flye_info <- read.csv2('input_data/all_fly_assembly_info.txt', sep = '\t')
# add column names
colnames(flye_info) <- c("seq_name","length","cov.","circ.","repeat.","mult.","alt_group","graph_path","sample")

# check how many on the chromosomes are not circular
flye_info$length <- as.numeric(as.character(flye_info$length))
flye_info$cov. <- as.numeric(as.character(flye_info$cov.))
ggplot(flye_info) + geom_point( aes(x=length, y = cov., col=circ.)) 

p <- ggplot(flye_info, aes(x=length/1000000)) + 
  geom_histogram(bins = 50, aes(fill=circ.)) + xlab('length [MB]') + ylab('Number of contigs') + theme_light() + scale_fill_manual(values = c('grey77', 'grey22'))

sum(flye_info$length > 500000 & flye_info$length < 2000000) #409 chromosome 1 assembled
table(flye_info[flye_info$length > 500000 & flye_info$length < 2000000,]$circ.) # 408/409 chromosome 1 assembled circularly
sum(flye_info$length > 2000000 & flye_info$length < 3700000) #409 chromosome 2 assembled
table(flye_info[flye_info$length > 2000000 & flye_info$length < 3700000,]$circ.) # 400/409 chromosome 1 assembled circularly

sum(flye_info$length < 50000 & flye_info$circ. == 'Y') #126 non-chromosomal small elements were circularly assembled
sum(flye_info$length > 3500000) # in 28 samples one large chromosome was assembled (instead of chr1 and chr2)

# 10 chromosomes of 9 isolates did not circularly assemble
not_all_chr_circ <- flye_info[flye_info$length > 500000 & flye_info$length < 3700000 & flye_info$circ. == 'N',]

# summarise to how many chromosomes were assembled
flye_info_sum <- flye_info %>% group_by(sample) %>% summarise(n_chr = sum(length > 800000), 
                                                              fract_circ = mean(circ.[length> 800000] == 'Y'))
flye_info_sum$n_chr <- factor(flye_info_sum$n_chr, levels = c(2,1))
flye_info_sum$fract_circ <- ifelse(flye_info_sum$fract_circ == 1, 'all circular', 
                                   ifelse(flye_info_sum$fract_circ == 0.5, '1/2 circular', 
                                          ifelse(flye_info_sum$fract_circ == 0, '0/2 circular', NA)))
flye_info_sum$fract_circ <- factor(flye_info_sum$fract_circ, levels = c('all circular', '1/2 circular', '0/2 circular'))
flye_info_sum_p <- ggplot(flye_info_sum, aes(x= n_chr, fill = factor(n_chr), alpha = factor(fract_circ))) +
  geom_bar() + scale_alpha_manual(values = c(1,0.6,0.2)) + scale_fill_manual(values = c('darkblue', 'lightblue')) +
  ylab('Nr. of assemblies') + xlab('Nr. of chromosomes') + theme_light()
flye_info_sum_p


# add column indicating whether (i) not all chromosomes where fully circularised (n=9) and where one fused chromosome was circularry assembled
QC['chr'] <- ifelse(QC$Sample %in% flye_info_sum[flye_info_sum$n_chr == '1', ]$sample, 'One chromosome assembled', 
                    ifelse(QC$Sample %in% not_all_chr_circ$sample, 'Not all chromosomes circularised', 'Two circular chromosomes assembled'))
# write 
write.table(QC, 'input_data/QC_extended.csv', row.names = F, quote = F, sep = ';')

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
meta_basic <- read.csv2("input_data/very_basic_meta_.csv", sep = ',')

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

