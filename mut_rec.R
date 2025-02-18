# load packages
library(dplyr)
library(tidyr)
library(ggplot2)

# setwd
setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# import SNV called by medaka and annotated by snpEff
snpEff_medaka_out <- read.table('input_data/snv_in_mut_rec_1.tab')
colnames(snpEff_medaka_out) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT' ,'SAMPLE_ratio')
# add isolate name
snpEff_medaka_out['isolate'] <- gsub('(\\..*)', '', snpEff_medaka_out$CHROM)
length(unique(snpEff_medaka_out$isolate)) # found between 3-5 SNV in each sample (I guess these are mainly phylogenetic)

# seperate annotation column
snpEff_medaka_out_ <- separate_wider_delim(snpEff_medaka_out, cols = INFO, too_many = "merge", delim = "|",names = c('Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos_cDNA.length', 'CDS.pos_CDS.length', 'AA.pos_AA.length', 'Distance', 'ERRORS_WARNINGS_INFO'))
# check range
range(snpEff_medaka_out_$QUAL)
# only consider QUAL > 40
snpEff_medaka_out_$QUAL <- as.numeric(as.character(snpEff_medaka_out_$QUAL))
snpEff_medaka_out_ <- snpEff_medaka_out_[snpEff_medaka_out_$QUAL > 40,]

#check if all protein coding
table(snpEff_medaka_out_$Transcript_BioType)
# only consider nonsynonymous
table(snpEff_medaka_out_$Annotation)
snpEff_medaka_out_ <- snpEff_medaka_out_[snpEff_medaka_out_$Annotation !='synonymous_variant',]

# import info on fused vs. non-fused
QC <- read.csv2('input_data/QC_extended.csv')
# summarise the ones where two large contigs have been assembled as carrying two chromosomes (independent of whether these were asselmbled to circular chromosomes or not, as I do not have evidence for these to be fused)
QC$chr <- ifelse(QC$chr %in% c('Not all chromosomes circularised','Two circular chromosomes assembled'), 2, 1)

# add fusion state to snpEff
snpEff_medaka_out_ <- merge(snpEff_medaka_out_, QC[,c('Sample', 'chr')], by.x = 'isolate', by.y='Sample', all.x=T)

unique(snpEff_medaka_out_$Gene_Name) # SNV detect in 4 genes "NDFNKM_03440" "NDFNKM_02795" "NDFNKM_11405" "NDFNKM_01770"
# check one by one if difference between fused an non-fused 
# check how many different variants
SNV_in_NDFNKM_03440 <- snpEff_medaka_out_[snpEff_medaka_out_$Gene_Name == 'NDFNKM_03440',]
table(SNV_in_NDFNKM_03440$HGVS.c, SNV_in_NDFNKM_03440$chr) # Just one SNV and occurs in all non-fused and all but one fused strains so not unique to fused 

SNV_in_NDFNKM_02795 <- snpEff_medaka_out_[snpEff_medaka_out_$Gene_Name == 'NDFNKM_02795',]
table(SNV_in_NDFNKM_02795$HGVS.c, SNV_in_NDFNKM_02795$chr) # Just one SNV and occurs in all strains so not unique to fused 

SNV_in_NDFNKM_11405 <- snpEff_medaka_out_[snpEff_medaka_out_$Gene_Name == 'NDFNKM_11405',]
table(SNV_in_NDFNKM_11405$HGVS.c, SNV_in_NDFNKM_11405$chr) # two different SNV occure only in 90/409 non-fused

SNV_in_NDFNKM_01770 <- snpEff_medaka_out_[snpEff_medaka_out_$Gene_Name == 'NDFNKM_01770',]
table(SNV_in_NDFNKM_01770$HGVS.c, SNV_in_NDFNKM_01770$chr) # one SNV occur only in 58/409 non-fused

# no SNV in fused chr-specific SNV in mut & rec genes found




