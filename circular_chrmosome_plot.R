library('circlize')
library('ggplot2')
library('dplyr')
library('tidyr')

# setwd
setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')


# import locations of ori, ctrS
# import oriC positions
oriC <- read.table('input_data//oriC.tsv', sep='\t')
colnames(oriC) <- c('SequenceId','Type','Start','Stop','Strand','LocusTag','Gene','Product')
oriC['Sample'] <- gsub('\\/.*', '', oriC$SequenceId)
oriC['Chr'] <- gsub('(.*\\:)(.*)', '\\2', oriC$SequenceId)
oriC <- oriC %>% 
  mutate(n_ori_per_chr = n(), 
         length_ori = Stop - Start) %>%
  group_by(Sample) %>% 
  mutate(n_ori_per_sample = n()) %>% # add how often found per sample and per chr
  group_by(Sample, Chr) %>%
  mutate(n_ori_per_chr = n()) 
oriC['ori'] <- ifelse(oriC$length_ori < 500, 'ori1', 'ori2')

oriC_plot <- oriC[oriC$Sample %in% c('114Vc01', '114Vc03'), c('Sample','Chr', 'Start','Stop', 'ori')]
colnames(oriC_plot) <- c('Sample','Chr', 'Start','Stop', 'Gene')

# import crtS 
blast_out_crtS <- read.table('input_data//blastout_crtS_site_household_genomes.tab')
colnames(blast_out_crtS) <-  c('qseqid', 'sseqid', 'bitscore', 'pident', 'nident', 'mismatch', 'length', 'qcovs', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send')
table(blast_out_crtS$pident, blast_out_crtS$qcovs) #
blast_out_crtS['Sample'] <- gsub('(\\d{3}Vc\\d{2})(\\d{1})', '\\1', blast_out_crtS$sseqid)
blast_out_crtS['Chr'] <- gsub('(\\d{3}Vc\\d{2})(\\d{1})', '\\2', blast_out_crtS$sseqid)

crtS_plot <- blast_out_crtS[blast_out_crtS$Sample %in% c('114Vc01', '114Vc03'), c('Sample','sseqid', 'sstart','send', 'qseqid')]
colnames(crtS_plot) <- c('Sample','Chr', 'Start','Stop', 'Gene')

# import HS1
# import blastout
HS1 <- read.table('input_data/blastout_HR_all.tab')
colnames(HS1) <-  c('qseqid', 'sseqid', 'bitscore', 'pident', 'nident', 'mismatch', 'length', 'qcovs', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send')
# only consider hits when > 95% of the query was found
HS1 <- HS1[HS1$length >= HS1$qlen*0.95,]
HS1['Sample'] <- gsub('(\\d{3}Vc\\d{2})(\\d{1})','\\1',HS1$sseqid)
# add orientation
HS1['orientation'] <- ifelse(HS1$sstart < HS1$send, 'F', 'R')

HS1_plot <- HS1[HS1$Sample %in% c('114Vc01', '114Vc03'), c('Sample','sseqid', 'sstart','send', 'qseqid')]
colnames(HS1_plot) <- c('Sample','Chr', 'Start','Stop', 'Gene')

# import positions of virulence isolands
PI <- read.csv2('input_data/island_positions_all.csv')
PI_plot <- PI[PI$strain %in% c('114Vc01', '114Vc03'), c('strain','sseqid', 's_min_loc','s_max_loc', 'island')]
colnames(PI_plot) <- c('Sample','Chr', 'Start','Stop', 'Gene')

# merge all
loci <- rbind(crtS_plot, HS1_plot, oriC_plot, PI_plot)

# in fused chromosome VPI-1 is broken through chr2. Therefore adjust boundaries
loci_VPI2 <- data.frame(Sample = '114Vc03', Chr = '114Vc031', Start = c(loci[loci$Chr == '114Vc031' & loci$Gene == 'VPI-2', 'Start'], max(loci[loci$Chr == '114Vc031' & loci$Gene == 'overlap_134Vc041_134Vc042', 'Stop'])),
                                                              Stop = c(min(loci[loci$Chr == '114Vc031' & loci$Gene == 'overlap_134Vc041_134Vc042', 'Start']), loci[loci$Chr == '114Vc031' & loci$Gene == 'VPI-2', 'Stop']), 
                        Gene =c('VPI-2_part1', 'VPI-2_part2'))

# replace VPI-2 coordinated with the 2 part coordinates
loci <- rbind(loci[!(loci$Chr == '114Vc031' & loci$Gene == 'VPI-2'),], loci_VPI2)
loci['Gene_label'] <- ifelse(loci$Gene == 'overlap_134Vc041_134Vc042', 'HS1', 
                             ifelse(loci$Gene == 'crtS_NSCV1', 'crtS', 
                                    ifelse(grepl('VPI-2', loci$Gene), 'VPI-2',loci$Gene)))

# Add chr identity
Chr_ID <- data.frame(Sample = c('114Vc01', '114Vc01', '114Vc03', '114Vc03', '114Vc03'), 
                     Chr = c('114Vc011', '114Vc012', '114Vc031', '114Vc031', '114Vc031'), 
                     Start = c(0,0,0,min(loci[loci$Chr == '114Vc031' & loci$Gene == "overlap_134Vc041_134Vc042",]$Stop) + 1, 
                                     max(loci[loci$Chr == '114Vc031' & loci$Gene == "overlap_134Vc041_134Vc042",]$Start)),
                     Stop = c(loci[loci$Chr == '114Vc011' & loci$Gene == 'ori1',]$Stop,
                              loci[loci$Chr == '114Vc012' & loci$Gene == 'ori2',]$Stop,
                              min(loci[loci$Chr == '114Vc031' & loci$Gene == "overlap_134Vc041_134Vc042",]$Stop),
                              max(loci[loci$Chr == '114Vc031' & loci$Gene == "overlap_134Vc041_134Vc042",]$Start) - 1,
                              loci[loci$Chr == '114Vc031' & loci$Gene == 'ori1',]$Stop), 
                     Gene = "",
                     Chr_ID = c('Chr1', 'Chr2', 'Chr1', 'Chr2', 'Chr1'))


# Set up the circular plot
# # # run the following chunk of code for each of the three chromosomes (fused 114Vc031, unfused chr1 114Vc011, unfused chr2 114Vc012)
# fused 
#Chromosome = '114Vc031'

# unfused chr 1 
#Chromosome = '114Vc011'

# unfused chr2
Chromosome = '114Vc012'

loci_plot  <- loci[loci$Chr == Chromosome,]
chromosome_length <- max(loci_plot$Stop)
loci_plot['col'] <- ifelse(grepl('^V', loci_plot$Gene), 'red',
                           ifelse(grepl('overlap', loci_plot$Gene), 'blue', 'lightgrey'))
# order loci plot such that HS1 is plotted last (as its part of VPI-2)
loci_plot <- loci_plot[c(which(loci_plot$Gene_label!='HS1'), which(loci_plot$Gene_label =='HS1')),]

chr_id_plot <- Chr_ID[Chr_ID$Chr == Chromosome,]
chr_id_plot['col'] <- ifelse(chr_id_plot$Chr_ID == 'Chr1', 'purple', 'orange')

circos.clear()  # Reset circos plot if re-running
circos.par("start.degree" = 90)  # Rotate so the start is at the top
circos.initialize(
  factors = "chromosome",
  xlim = c(0, chromosome_length)
)

# Add the chromosome track
circos.trackPlotRegion(
  factors = "chromosome",
  ylim = c(0, 1),
  track.height = 0.1,
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ycenter + 0.1,
      "",
      facing = "bending.outside",
      cex = 0.8
    )
  }
)

# Add Chromosome ID as segments
for (i in 1:nrow(chr_id_plot)) {
  circos.segments(
    sector.index = "chromosome",
    x0 = chr_id_plot$Start[i],
    x1 = chr_id_plot$Stop[i],
    y0 = 0.5,
    y1 = 0.5,
    col = chr_id_plot$col[i],
    lwd = 5
  )
}

# Add labels for chromosome ID
for (i in 1:nrow(chr_id_plot)) {
  circos.text(
    x = (chr_id_plot$Start[i] + chr_id_plot$Stop[i]) / 2,
    y = -2,
    labels = chr_id_plot$Chr_ID[i],
    sector.index = "chromosome",
    facing = "clockwise",
    niceFacing = TRUE,
    col = chr_id_plot$col[i],
    cex = 0.7
  )
}


# Add loci as segments
for (i in 1:nrow(loci_plot)) {
  circos.segments(
    sector.index = "chromosome",
    x0 = loci_plot$Start[i],
    x1 = loci_plot$Stop[i],
    y0 = 0.5,
    y1 = 0.5,
    col = loci_plot$col[i],
    lwd = 6
  )
}

# Add labels for loci
for (i in 1:nrow(loci_plot)) {
  circos.text(
    x = (loci_plot$Start[i] + loci_plot$Stop[i]) / 2,
    y = 2,
    labels = loci_plot$Gene_label[i],
    sector.index = "chromosome",
    facing = "clockwise",
    niceFacing = TRUE,
    col = loci_plot$col[i],
    cex = 0.7
  )
}

#dev.print(pdf, 'output_figures/fused_chr_scheme_plot.pdf', width = 2.4, height = 2.4)
#dev.print(pdf, 'output_figures/non-fused_chr1_scheme_plot.pdf', width = 1.6, height = 1.6)
dev.print(pdf, 'output_figures/non-fused_chr2_scheme_plot.pdf', width = 0.8, height = 0.8)

# # #

# Add an axis to represent the genome position
#circos.axis(
#  h = "top",
#  labels.cex = 0.6,
#  major.at = seq(0, chromosome_length, by = 1e6),
#  labels = paste0(seq(0, chromosome_length, by = 1e6) / 1e6, " Mb")
#)


# export positions of HS1 flanking sides as bed files
HS1_chr_1 <- HS1[HS1$slen > 2000000 & HS1$slen < 3500000,]
HS1_flanking_1 <- data.frame(chr = HS1_chr_1$sseqid, start = HS1_chr_1$sstart - 500, stop = HS1_chr_1$sstart, id='HS1_flanking_1')
HS1_flanking_2 <- data.frame(chr = HS1_chr_1$sseqid, start = HS1_chr_1$send, stop = HS1_chr_1$send + 500, id='HS1_flanking_2')
HS1_chr_2 <- HS1[HS1$slen < 2000000,]
HS1_flanking_3 <- data.frame(chr = HS1_chr_2$sseqid, start = HS1_chr_2$sstart - 500, stop = HS1_chr_2$sstart, id='HS1_flanking_3')
HS1_flanking_4 <- data.frame(chr = HS1_chr_2$sseqid, start = HS1_chr_2$send, stop = HS1_chr_2$send + 500, id='HS1_flanking_4')

write.table(HS1_flanking_1, 'input_data/check_gene_transfer/HS1_flanking_1.bed', quote = F, sep = '\t', col.names = F, row.names = F)
write.table(HS1_flanking_2, 'input_data/check_gene_transfer/HS1_flanking_2.bed', quote = F, sep = '\t', col.names = F, row.names = F)
write.table(HS1_flanking_3, 'input_data/check_gene_transfer/HS1_flanking_3.bed', quote = F, sep = '\t', col.names = F, row.names = F)
write.table(HS1_flanking_4, 'input_data/check_gene_transfer/HS1_flanking_4.bed', quote = F, sep = '\t', col.names = F, row.names = F)

# I checked, for HS1_flanking_1, HS1_flanking_2 and HS1_flanking_4, the values are all > 99.9. --> can be missing values if below 70 --> eval all

ANI_1 <- read.table('input_data/check_gene_transfer/ANI_HS1_flanking_1.tab')
colnames(ANI_1) <- c('query', 'ref', 'ANI', 'nr_mapped_fragments', 'nr_total_fragments')
ANI_1$query <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_1$query)
ANI_1$ref <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_1$ref)
length(unique(ANI_1$query))
length(unique(ANI_1$ref))
nrow(ANI_1) == 407*407
range(ANI_1$ANI) #--> are all the same
ANI_1['seq'] <- 'flanking_1'

# create df with all combinations
all_combi_chr1 <- expand.grid(unique(ANI_1$query),unique(ANI_1$ref))
colnames(all_combi_chr1) <- c('query', 'ref')

ANI_2 <- read.table('input_data/check_gene_transfer/ANI_HS1_flanking_2.tab')
nrow(ANI_2) == 407*407
# merge to all combi df
colnames(ANI_2) <- c('query', 'ref', 'ANI', 'nr_mapped_fragments', 'nr_total_fragments')
ANI_2$query <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_2$query)
ANI_2$ref <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_2$ref)
ANI_2 <- merge(ANI_2, all_combi_chr1, by=c('query', 'ref'), all = T)
ANI_2['seq'] <- 'flanking_2'
#plot 
ggplot(ANI_2, aes(ref, query, fill= ANI)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", na.value = 'grey') 

ANI_2$ANI <- ifelse(is.na(ANI_2$ANI), 69, ANI_2$ANI)
sort(table(ANI_2[ANI_2$ANI < 99.5,]$ref)) #--> 112Vc05 against all which is below 99.9 


ANI_3 <- read.table('input_data/check_gene_transfer/ANI_HS1_flanking_3.tab')
colnames(ANI_3) <- c('query', 'ref', 'ANI', 'nr_mapped_fragments', 'nr_total_fragments')
ANI_3$query <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_3$query)
ANI_3$ref <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_3$ref)
# create df with all combinations chr2 
all_combi_chr2 <- expand.grid(unique(ANI_3$query),unique(ANI_3$ref))
colnames(all_combi_chr2) <- c('query', 'ref')
nrow(ANI_3) == nrow(all_combi_chr2)
ANI_3 <- merge(ANI_3, all_combi_chr2, by=c('query', 'ref'), all = T)
#plot 
ggplot(ANI_3, aes(ref, query, fill= ANI)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", na.value = 'grey') 
ANI_3$ANI <- ifelse(is.na(ANI_3$ANI), 69, ANI_3$ANI)
sort(table(ANI_3[ANI_3$ANI < 99.5,]$ref)) 
#--> 103Vc01 and 112Vc05 are the same
#--> 108Vc03 and 109Vc05 are the same
# both are distant to all other and betwen


ANI_4 <- read.table('input_data//check_gene_transfer/ANI_HS1_flanking_4.tab')
colnames(ANI_4) <- c('query', 'ref', 'ANI', 'nr_mapped_fragments', 'nr_total_fragments')
ANI_4$query <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_4$query)
ANI_4$ref <- gsub('(/lustre04/scratch/acuenod/household/01_data/21_check_chr_dna_transfer/flanking_HS1/seqs/HS1_flanking_./)(.*)(.fasta)', '\\2', ANI_4$ref)
nrow(ANI_4) == nrow(all_combi_chr2)
ANI_4 <- merge(ANI_4, all_combi_chr2, by=c('query', 'ref'), all = T)
#plot 
ggplot(ANI_4, aes(ref, query, fill= ANI)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", na.value = 'grey') 
ANI_4$ANI <- ifelse(is.na(ANI_4$ANI), 69, ANI_3$ANI)
sort(table(ANI_4[ANI_4$ANI < 99.5,]$ref)) 
#--> 103Vc01 and 112Vc05 are the same
#--> 108Vc03 and 109Vc05 are very similar
# both are distant to all other and betwen

# check 103Vc01, 112Vc05, 108Vc03 and 109Vc05 in more detail
# compare them each to a closely related genome (with 2 chr) sampled from the same household
# 103Vc01 to 103Vc05 inversion of seq on chr2 which includes HS1
# 112Vc05 to 114Vc09 large rearrangement with big exchanges between chr1 and chr2
# 108Vc03 to 109Vc02 --> rearrangement just upstream of HS1
# 109Vc05 to 109Vc02 --> rearrangement just upstream of HS1
# compared with minlength 1,000


