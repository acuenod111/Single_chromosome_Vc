library('data.table')
library('ggplot2')
library('dplyr')
library('tidyr')

# setwd
setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

## eval flanking region blast
# import the blast output from the virulence islands flanking sides
path <- 'input_data/island_flanking_blastout/'
files <- list.files(path = path, pattern = '.tab')
temp <- lapply(paste0(path, files), fread, sep="\t")

# merge to one file
blastout_q_full <- do.call(rbind, lapply(paste0(path, files), function(x) 
  transform(fread(x), query = gsub('.tab','',basename(x)))))
colnames(blastout_q_full) <-  c('qseqid', 'sseqid', 'bitscore', 'pident', 'nident', 'mismatch', 'length', 'qcovs', 'qlen', 'qstart', 'qend', 'qseq', 'slen', 'sstart', 'send', 'sseq', 'query')

blastout_q_full <- as.data.frame(blastout_q_full)
# remove seqs, are too heavy
blastout_q_full$qseq <-NULL
blastout_q_full$sseq <-NULL
# add strain and island
blastout_q_full$strain <- gsub('(\\d{3}Vc\\d{2})(\\d)', '\\1', blastout_q_full$sseqid)
blastout_q_full['island'] <- gsub('\\_.*', '', blastout_q_full$query)

table(blastout_q_full$query) 
# all queries where found in 467 genomes (in each genome once), except for VPI-1_A1552VC_RS03055, which was found twice per genome and VSP-II which was missing in one genome
ggplot(blastout_q_full[blastout_q_full$query == 'VPI-1_A1552VC_RS03055', ]) + geom_histogram(aes(x=pident))

# of these one is always much higher than the second hit, choose one hit per query and genome
blastout_q_VP <- blastout_q_full %>% group_by(query, strain) %>% slice_max(pident)
table(blastout_q_VP$query) 
range(blastout_q_VP$pident)
range(blastout_q_VP$qcovs) # VSP-II_A1552VC_RS01290 is only found with qcov if 66% in 115Vc09 (all others > 99%) (not at the end of a contig)
# find which genome is missing VSP-II
setdiff(blastout_q_VP[blastout_q_VP$query != 'VSP-II_A1552VC_RS01435',]$strain, blastout_q_VP[blastout_q_VP$query == 'VSP-II_A1552VC_RS01435',]$strain)
# Its 141Vc01. (has two contig and looks normal QC wise)

# convert to long to get range of genome position per island and genome
island_positions <- blastout_q_VP %>% 
  pivot_longer(cols=c(sstart,send), names_to = 's_which_location', values_to = 's_genome_location') %>%
  #get max and min position per genome and island
  group_by(strain, island) %>%
  mutate(s_min_loc = min(s_genome_location), 
         s_max_loc = max(s_genome_location)) %>%
  select(sseqid, strain, island, s_min_loc, s_max_loc) %>%
  distinct() # %>% group_by(strain, island) %>% mutate(n=n())

# Check length of supposed island
island_positions['length'] <- island_positions$s_max_loc - island_positions$s_min_loc
ggplot(island_positions) + geom_histogram(aes(x=length)) + facet_grid(island~., scales = 'free_x')
# check the ones which are > 500'000
island_positions_check <- island_positions[island_positions$length > 500000,]
island_positions_check$strain[duplicated(island_positions_check$strain)] # in 112Vc05 VPI-2 is split onto the two choromosomes (one flanking side one chr 1 and one on chr 2). I think this is actually also a 1xchr, as I found reads spanning the 12kb junction (similarly as all other 112Vc, which are 1x chr.)
# keep only one to count
island_positions_check <- island_positions_check[!island_positions_check$sseqid == '112Vc052',]
# this is very akward. VPI-2 us in 58 ca. 1,130,000 bp long. Check if these are (a) primarily on genomes with a fused (4m chromosome) or (b) on chromosomes which did not fully circularise 
blastout_q_full[blastout_q_full$strain %in% island_positions_check$strain & blastout_q_full$island == 'VPI-2',]
check <- blastout_q_full[blastout_q_full$strain %in% island_positions_check$strain & blastout_q_full$island == 'VPI-2',] 
# the long ones are all in fused chromosomes
# in 112Vc05 the two flanking regions where found on two different chromosomes
check <- check %>% group_by(strain) %>% slice_min(slen, with_ties = F)

# import QC
QC <- read.table("input_data/QC_extended.csv", sep=';', header = T)
QC['VPI-2_long'] <- ifelse(QC$Sample %in% check$strain, 'yes', 'no')
# 58 strains for which one fused chromosome was assembled have 'long' VP-2 and all straings with a 'long' VPI-2 have only one chromosome assembled (except for 112Vc05, where VPI-2 flanking regions are split between the two chromosomes)
table(QC$chr, QC$`VPI-2_long`)
# I figured out that this is because VPI-2 get split upon fusion (by the chr2 insertion) as HS1 is part of VPI-2
ggplot(island_positions) + geom_histogram(aes(x=length), bins = 20) + facet_grid(.~island, scales = 'free')

write.csv2(island_positions , 'input_data/island_positions_all.csv', quote = F, row.names = F)
