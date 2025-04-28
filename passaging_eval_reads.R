library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# import the reads spanning the HS1 region
reads_spanning <- read.table('input_data/passaging/count_reads_HS1_passaging.txt')
# remove last line, is total count
reads_spanning <- reads_spanning[-nrow(reads_spanning),]
colnames(reads_spanning) <- c('nr_reads','path')
# add if HS1 was in fused (chr1-HS1-chr2/chr2-HS1-chr1) or non-fused state (chr1-HS1-chr1/chr2-HS1-chr2)
reads_spanning['mapping_HS1'] <- gsub('(readmapping_)(.*fused)(_HR.*)', '\\2', reads_spanning$path)
# add sample, isolate and date (from the sequence name)
reads_spanning['sample'] <- gsub('(.*\\/)([[:alnum:]]*\\_*[[:alnum:]]*)(_reads_spanning_HR(1|2).txt)', '\\2', reads_spanning$path)
reads_spanning['isolate'] <- gsub('(.*)(D\\d{2})$', '\\1', reads_spanning$sample)
reads_spanning['day'] <- gsub('(.*)(D)(\\d{2})$', '\\3', reads_spanning$sample)
# add fusion state
reads_spanning['fusion_state'] <- ifelse(reads_spanning$isolate %in% c('117Vc01','114Vc03'), 'fused', 'non-fused')
# add the number of estimated generations
# discussed this again, the more conservative estimate is 10 generations per day (as in Bruhn et al.)
reads_spanning['generations_estimate'] <- ifelse(reads_spanning$day == '00', '00', 
                                                 ifelse(reads_spanning$day == '03', '30', 
                                                        ifelse(reads_spanning$day == '08', '80', 
                                                               ifelse(reads_spanning$day == '11', '110', 
                                                                      ifelse(reads_spanning$day == '16', '160', 
                                                                             ifelse(reads_spanning$day == '20', '200', reads_spanning$day))))))


# I checked if any of these reads contain an adaptor (fake fused reads) (by converting them to fasta seqs and blasting the adaptor seqs). No blasthit was found, so these are all genuine reads
# import total read length to normalise (from flye output)
total_read_length <- read.table('input_data/passaging/total_read_length_info_sample.txt')
colnames(total_read_length) <- c('barcode', 'total_read_length', 'sample')
# merge
reads_spanning <- merge(reads_spanning, total_read_length, by='sample')
reads_spanning['nr_reads_norm'] <- reads_spanning$nr_reads/reads_spanning$total_read_length

# remove TE Buffer reads (none mapped against HS1)
reads_spanning <- reads_spanning[reads_spanning$sample != 'TE_Buffer',]

# convert generation times into factors
reads_spanning$generations_estimate <- factor(reads_spanning$generations_estimate, levels=c( "00","30","80","110","160","200","BS01","N16961"))
# plot
p <- ggplot(reads_spanning) +
  geom_boxplot(aes(x=generations_estimate, y=nr_reads_norm, fill=mapping_HS1, col=mapping_HS1), alpha= 0.5) + scale_fill_manual(values = c('lightblue', 'darkblue'), na.value = "grey91") + scale_color_manual(values = c('lightblue', 'darkblue'), na.value = "grey91") + 
  geom_point(aes(x=generations_estimate, y=nr_reads_norm, col=mapping_HS1), position = position_dodge(width = .75), shape = 15) +
  facet_grid(~fusion_state) + 
  ylab('# reads spanning HS1 / total read length') +
  xlab('Estimated # generations') + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 60)) 
p  

pdf("output_figures/passaging_read.pdf", height = 3.5, width = 7)
p
dev.off()
#

# I first called these with medaka which gave me many SNV (which I did not trust, because there were so many). I therefore repeated the variant calling with clair3, which also gives me information on the allele frequency which I want to check 
clair3 <- read.table('input_data/passaging/clair3_25_sum.tab', sep='\t')
# add col names
colnames(clair3) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE', 'sample')
clair3['ref'] <- gsub('(\\d{3}Vc\\d{2})(\\1)', '\\1', clair3$CHROM)
clair3['Allele_Frequency_ALT'] <- gsub('.*\\:', '', clair3$SAMPLE)
clair3['Estimated_read_depth'] <- gsub('(\\d+\\/\\d+\\:)(\\d+\\:)(\\d+)(\\:\\d+\\,\\d+)(\\:.*)', '\\3', clair3$SAMPLE)
# add isolate and day
clair3['isolate'] <- gsub('(.*)(D\\d{2})$', '\\1', clair3$sample)
clair3['day'] <- gsub('(.*)(D)(\\d{2})$', '\\3', clair3$sample)
# add fusion state
clair3['fusion_state'] <- ifelse(clair3$isolate %in% c('117Vc01','114Vc03'), 'fused', 'non-fused')
# remove the low qual ones
clair3 <- clair3[clair3$FILTER !='LowQual',]

# plot 
# check which SNV appear more than once
clair3 <- clair3 %>% group_by(CHROM, POS, REF, ALT) %>% 
  mutate(n_occurence_per_snv = n(), #check how many SNV per sample 
         same_SNV = cur_group_id()) %>% #check how many SNV per sample 
  ungroup() %>%
  group_by(sample) %>%
  mutate(n_snv_per_sample= n())

clair3$Allele_Frequency_ALT <- as.numeric(as.character(clair3$Allele_Frequency_ALT))

# plot the fraction of each SNV
#ggplot(clair3, aes(y=Allele_Frequency_ALT, x=day, colour = same_SNV, group=1)) +
#  geom_point() + geom_line() + theme(legend.position = 'none') + facet_grid(isolate~.)

clair3$isolate <- factor(clair3$isolate, levels = c('114Vc03', '117Vc01', '114Vc01', '117Vc10'))

snv_phase <- ggplot(clair3, aes(y=Allele_Frequency_ALT, x=day, group = same_SNV)) +
  geom_line(aes(col=fusion_state)) + 
  geom_point(alpha=0.2, aes(col=fusion_state)) +
  theme_light() +
  theme(legend.position = 'none') + facet_grid(isolate~.) + scale_color_manual(values = c('lightblue', 'darkblue'), na.value = "grey91")
snv_phase

n_snv <- ggplot(clair3, aes(x=day, y=n_snv_per_sample, col=fusion_state)) +
  geom_point() +
  #geom_text(aes(label = isolate), hjust=0, vjust=0) + 
  theme_light() + facet_grid(isolate~.) + scale_color_manual(values = c('lightblue', 'darkblue'), na.value = "grey91") + ggtitle('Clair3 - all')

n_snv + snv_phase

# all come from a single col isolate, so I will apply the following cut-offs to confidently call a SNV: read depth >= 10; allele frequency >= 0.8 and quality >= 20
clair3_hq_col <- clair3[clair3$Estimated_read_depth >=10 & clair3$Allele_Frequency_ALT >=0.8 & clair3$QUAL >=20,]

# check which SNV appear more than once
clair3_hq_col <- clair3_hq_col %>% group_by(CHROM, POS, REF, ALT) %>% 
  mutate(n_occurence_per_snv = n(), #check how many SNV per sample 
         same_SNV = cur_group_id()) %>% #check how many SNV per sample 
  ungroup() %>%
  group_by(sample) %>%
  mutate(n_snv_per_sample= n())

#summarise
clair3_hq_col_sum <- clair3_hq_col %>% group_by(sample, fusion_state) %>% 
  summarise(n_occurencsnv_per_sample = n()) #check how many SNV per sample 
# check of there are more SNV for the fused vs. the non-fused
wilcox.test(clair3_hq_col_sum[clair3_hq_col_sum$fusion_state == 'fused',]$n_occurencsnv_per_sample, clair3_hq_col_sum[clair3_hq_col_sum$fusion_state == 'non-fused',]$n_occurencsnv_per_sample, paired = F)

# check if the ones which are not in fell out because they have actually no SNV or bcause they were low cov
setdiff(clair3$sample, clair3_hq_col$sample)
# Exclude these timepoints because there was not enough read data sequenced
low_cov <- c('117Vc01D00','114Vc03D03','114Vc01D08','114Vc03D08','117Vc01D08','114Vc03D16') # see table S3 for reference
setdiff(setdiff(clair3$sample, clair3_hq_col$sample), low_cov)
# add this ("114Vc01D00") with number of SNV = 0
to_add <- data.frame(sample = '114Vc01D00', isolate = '114Vc01', fusion_state = 'non-fused', day='00', n_snv_per_sample = 0)
clair3_hq_col <- rbind(clair3_hq_col, to_add)

# remove the ones with an estimated read depth below 10 
clair3_hq_col <- clair3_hq_col[!clair3_hq_col$sample %in% low_cov,]
clair3_hq_col['participant'] <- gsub('(\\d{3})(Vc.*)', '\\1', clair3_hq_col$sample)

n_hq_snv <- ggplot(clair3_hq_col, aes(x=day, y=n_snv_per_sample, col=fusion_state)) +
  geom_point() +
  #geom_text(aes(label = isolate), hjust=0, vjust=0) + 
  theme_light() + facet_grid(participant~fusion_state) + scale_color_manual(values = c('lightblue', 'darkblue'), na.value = "grey91") + 
  #ggtitle('Clair 3 - hq\n(Estimated_read_depth >=10 &\nAllele_Frequency_ALT >=0.8 &\nQUAL >=20)') 
  ylab('Nr. SNV per sample') + xlab('')
n_hq_snv

pdf('output_figures//passaging_SNV.pdf', height = 3, width = 4.5)
n_hq_snv
dev.off()

# I had previously run clair3 with phasing (by mistake). For these SNV I had annotated them with SNPEff and checked the bakta annotation
# These are not not 100% consistent with the above (clair3 without phasing)
## import annotations (from snpEff)
#snpEff <- read.table('../../01_data/22_passaging/snpEff_sum.tab')
#colnames(snpEff) <- colnames(clair3) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','ANNOTATION', 'FORMAT','SAMPLE', 'sample')
#
#clair3_hq_col_annot <- merge(clair3_hq_col, snpEff, by = intersect(colnames(clair3_hq_col), colnames(snpEff)), all.x = T)
## seperate annotation column
#clair3_hq_col_annot_wide <- separate_wider_delim(clair3_hq_col_annot, cols = ANNOTATION, too_many = "merge", delim = "|",names = c('Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos_cDNA.length', 'CDS.pos_CDS.length', 'AA.pos_AA.length', 'Distance', 'ERRORS_WARNINGS_INFO'))
#
##check how many ncRNA
#ncRNA <- clair3_hq_col_annot_wide[clair3_hq_col_annot_wide$Gene_Name == 'V_AS5',]
#table(ncRNA$isolate, ncRNA$fusion_state)
## ncRNA 
#
## add gene info from bakta
#bakta_snv <- read.table('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/01_data/22_passaging/genes_with_snv.tab', sep='\t', quote = "", stringsAsFactors = FALSE)
#bakta_snv <- as.data.frame(bakta_snv)
#colnames(bakta_snv) <- c('SequenceId','Type','Start','Stop','Strand','LocusTag','Gene','Product','DbXrefs')
#bakta_snv['Ref_sample'] <- gsub('(.*\\:)(\\d{3}Vc\\d{2})(.)', '\\2', bakta_snv$SequenceId)
#
## remove the ncRNA ones, are all the same
#unique(bakta_snv[bakta_snv$Type == 'ncRNA',]$ProductDbXrefs)
#bakta_snv <- bakta_snv[bakta_snv$Type != 'ncRNA', ]
#
#clair3_hq_col_annot_wide_bakta <- merge(clair3_hq_col_annot_wide, bakta_snv, by.x='Gene_Name', by.y='LocusTag',)
## check if there are any missing other than ncRNA
#setdiff(clair3_hq_col_annot_wide$Gene_Name, clair3_hq_col_annot_wide_bakta$Gene_Name) # thats ok, the NA one is the line I added to plot for the sample with 0 SNV detected
#
## check if there are some which mutated in more than one sample
#clair3_hq_col_annot_wide_bakta <- clair3_hq_col_annot_wide_bakta %>% group_by(Product) %>% mutate(n_sample_per_product = n_distinct(Ref_sample))
## there are only two genes which are mutated in more than one sample: 
## --> I see (acfD, accessory colonization factor AcfD) mutated (missense mutations) in 3 experiments, one fused two unfused
## --> I see ompT mutated (synonymous SNV) in one fused and one non-fused sample
