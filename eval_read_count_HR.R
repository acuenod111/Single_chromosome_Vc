library('ggplot2')
library('tidyr')
library('dplyr')
library('seqinr')
library('data.table')
library('gggenes')
library('RColorBrewer')

setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# # # Evaluate blast from HS1 (12kb) queried against all genomes

# import blastout
overlap <- read.table('input_data/blastout_HR_all.tab')
colnames(overlap) <-  c('qseqid', 'sseqid', 'bitscore', 'pident', 'nident', 'mismatch', 'length', 'qcovs', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send')
# only consider hits when > 95% of the query was found
overlap <- overlap[overlap$length >= overlap$qlen*0.95,]
overlap['Sample'] <- gsub('(\\d{3}Vc\\d{2})(\\d{1})','\\1',overlap$sseqid)
# count the number of overlaps
overlap_sum <- overlap %>% group_by(Sample) %>% summarise(n_HR=n(), 
                                                          on_x_chr = length(unique(sseqid)))
table(overlap_sum$n_HR, overlap_sum$on_x_chr)
write.table(overlap_sum, 'input_data/blastout_HR_all_sum.tab', quote = F, row.names = F)

# export the HS1 location for each genome
for (i in overlap$Sample){
  subset <- NULL
  subset <- overlap[overlap$Sample == i, c("sseqid", "sstart", 'send')]
  subset['strand'] <- '+'
  colnames(subset) <- c("seqid","start","end","strand")
  write.table(subset, paste0('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/01_data/16_check_singular_contig/HS1_tab/', i, '.tab'), row.names = F, quote = F, sep = '\t')
}

# # # plot HS1 annotation (annotated by bakta)
# plot gene content on 12kb HR seq
HR_annotation <- read.table('input_data/annotation_overlap_134Vc041_134Vc042.tsv', comment.char="#", sep='\t', header = T)

HR_annotation['start_strand'] <- ifelse(HR_annotation$Strand == '-', HR_annotation$Stop, HR_annotation$Start)
HR_annotation['stop_strand'] <- ifelse(HR_annotation$Strand == '-', HR_annotation$Start, HR_annotation$Stop)

# count how many colors will be needed 
colourCount = length(unique(HR_annotation$Product))
# define getPalette function to get more colors than originally in Set3 pallette
getPalette = colorRampPalette(brewer.pal(9, "Set3"))

# change order of product (thereby legend) to fit the plot
HR_annotation$Product <- factor(HR_annotation$Product, levels = unique(HR_annotation[order(HR_annotation$Start),]$Product))

p <- ggplot(HR_annotation, aes(xmin = start_strand, xmax = stop_strand, y = Sequence.Id, fill = Product)) +
  geom_gene_arrow() +
  geom_gene_label(aes(label=Gene)) + scale_fill_manual(values = getPalette(colourCount)) + theme_classic() +
  theme(legend.position = 'bottom') 

pdf('output_figures/HR_gene_content.pdf', height=5, width=10)
p
dev.off()

# # # evaluate the count of reads which span the HS1
# import the output from when blasting the squencing adpaters against the reads which span the entire HS1
path <- 'input_data/adapterblast_spanning_reads/'
files <- list.files(path = path, pattern = '.tab')
temp <- lapply(paste0(path, files), fread, sep="\t")

blastout_q_full <- do.call(rbind, lapply(paste0(path, files), function(x) 
  transform(fread(x), query = gsub('.tab','',basename(x)))))
colnames(blastout_q_full) <-  c('qseqid', 'sseqid', 'bitscore', 'pident', 'nident', 'mismatch', 'length', 'qcovs', 'qlen', 'qstart', 'qend', 'qseq', 'slen', 'sstart', 'send', 'sseq', 'strand','query')
# decide for 95% qcovs and pident 
ggplot(blastout_q_full, aes(x=qcovs)) + geom_histogram(binwidth = 1) + geom_vline(xintercept = 95)
blastout_q <- blastout_q_full[blastout_q_full$qcovs > 95,]

ggplot(blastout_q, aes(x=pident)) + geom_histogram(binwidth = 1) + geom_vline(xintercept = 95)
blastout_q <- blastout_q[blastout_q$pident > 80,]

# F and R primer always found on +/- on same position, only consider F
# blastout_q <- blastout_q[blastout_q$query == 'adapter_F',] --> will anyways only export unique read ID
length(unique(blastout_q$sseqid))

# export in which reads an adapter was found. These were excluded from further analysis (these could potentially be hybrid reads (originating from two pieces of DNA which were accidentally fused while sequencing))
write.table(unique(blastout_q$sseqid), 'input_data/adapterblast_spanning_reads/read_id_adapter_found.txt', row.names = F, quote = F, col.names = F)

# import the read count over 12kb HR region. These are reads which span HS1 and 500pb either on a fused chromosome (chr1 on one side and chr2 on the other side) or on non-fused chromosomes (see below)
read_count_HR_junction <- read.table('input_data/read_count_spanning_HR_junction_adapter_excl.txt')
colnames(read_count_HR_junction) <- c('read_count','file')
read_count_HR_junction <- read_count_HR_junction[!read_count_HR_junction$file == 'total',]
read_count_HR_junction['Sample'] <- gsub('\\_.*', '', read_count_HR_junction$file)
read_count_HR_junction['region'] <- gsub('(\\d{3}Vc\\d{2}_reads_spanning_)(junction\\d{1})(.txt)', '\\2', read_count_HR_junction$file)
read_count_HR_junction$file <- NULL

# import read count non_fused
read_count_HR_non_fused <- read.table('input_data/read_count_spanning_HR_non-fused_adapter_excl.txt')
colnames(read_count_HR_non_fused) <- c('read_count','file')
read_count_HR_non_fused <- read_count_HR_non_fused[!read_count_HR_non_fused$file == 'total',]
read_count_HR_non_fused['Sample'] <- gsub('\\_.*', '', read_count_HR_non_fused$file)
read_count_HR_non_fused['region'] <- gsub('(\\d{3}Vc\\d{2}_reads_spanning_)(HR)(\\d{1})(.txt)', '\\2\\_chr\\3', read_count_HR_non_fused$file)
read_count_HR_non_fused$file <- NULL

read_counts <- merge(read_count_HR_junction, read_count_HR_non_fused, by=c('Sample', 'region', 'read_count'), all = T)

# import QC
QC <- read.csv('input_data/QC_extended.csv', sep = ';')
# check before merging
setdiff(QC$Sample, read_count_HR_junction$Sample)
read_counts <- merge(read_counts, QC[,c('Sample', 'Read_depth', 'chr')], by='Sample')

# import lineage assignment (medaka SNV calls)
lineage <- read.table('input_data/lineage_assignemnet_from_medaka.txt', header = T)
# summarise the 'Ind' lineages
lineage['lineage'] <- ifelse(lineage$Lineage_group %in% c('IND1.2', 'IND1.3'), 'IND1.2-3', lineage$Lineage_group)
# merge
read_counts <- merge(read_counts, lineage[,c('sample','lineage')], by.x = 'Sample', by.y = 'sample', all =T)
# normalise the read counts by the average genome coverage
read_counts['read_counts_norm'] <- read_counts$read_count / read_counts$Read_depth
# smmarise the assemblies for which two large contigs were assembled (circular or not) (there is no evidence that these are fused as opposed to the fused ones)
read_counts$chr <- ifelse(read_counts$chr %in% c('Not all chromosomes circularised', 'Two circular chromosomes assembled'), 'Two chromosomes assembled', read_counts$chr)
read_counts$chr <- factor(read_counts$chr, levels = c('Two chromosomes assembled', 'One chromosome assembled'))
read_counts$region <- ifelse(read_counts$region == 'HR_chr1', 'chr1-HR-chr1', 
                             ifelse(read_counts$region == 'HR_chr2', 'chr2-HR-chr2', 
                                    ifelse(read_counts$region == 'junction1', 'chr1-HR-chr2', 
                                           ifelse(read_counts$region == 'junction2', 'chr2-HR-chr1', NA))))
read_counts$region <- factor(read_counts$region, levels = c('chr1-HR-chr1', 'chr2-HR-chr2', 'chr1-HR-chr2', 'chr2-HR-chr1'))

# One 2xChr sample (112Vc05) has 5 and 4 reads spanning junction 1 and 2, respectively. All other 112Vc sample, one chr. only was assembled

p <- ggplot(read_counts, aes(x=region, y=read_counts_norm)) +
  #geom_boxplot(aes(fill=lineage)) + facet_grid(.~chr) + 
  geom_boxplot(position = position_dodge(preserve = "single"), aes(color=lineage, fill= lineage, alpha=0.9)) + scale_color_manual(values = c("#FDD58C", "#B2258F")) + scale_fill_manual(values = c("#FDD58C", "#B2258F")) +
  facet_grid(chr~.) + 
  ylab('Nr. reads spanning region / Mean read depth') + theme_light() + 
  theme(axis.text.x=element_text(angle=60,hjust=1), legend.position = 'bottom')


pdf('output_figures/read_count_HR_regions_excl_adapter.pdf', height=5, width=2.7)
p
dev.off()

