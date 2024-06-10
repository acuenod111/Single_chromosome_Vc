# public long read analysis
# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggforce)

setwd('/Users/alinecuenod/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alines_MacBook_Pro/Other/cholera/01_household_study/02_scripts/check/')

# import QC measures
QC <- read.table('input_data/public/long_reads/quality_4b.csv', sep=',', fill = TRUE, header = T, row.names = NULL)
QC <- QC[!grepl('analysed_samples', QC$Sample),]
samples <-  read.table('input_data/public/long_reads/samples_3.txt', sep=',', fill = TRUE, header = F, row.names = NULL)
setdiff(samples$V1, QC$Sample) # none

# merge so to have all 380 are in
QC <- merge(QC, samples, by.x = 'Sample', by.y = 'V1', all = T)
# some reads were sub-sampled and then assembled again, the sample list only includes the subsampled results, so does the QC file

# import species ID by GTDB-tk
gtdb <- read.table('input_data/public/long_reads/gtdbtk.bac120.summary_b.tsv', sep='\t', fill = TRUE, header = T, row.names = NULL)
gtdb_added <- read.table('input_data/public/long_reads/gtdbtk.bac120.summary_added_b.tsv', sep='\t', fill = TRUE, header = T, row.names = NULL)
gtdb <- merge(gtdb, gtdb_added, by = colnames(gtdb_added), all = T)
gtdb['species_classification'] <- gsub('.*s\\_{2}', '', gtdb$classification)
# Aliivibrio and Photobacterium are genera of the Vibrionales family
# check before merging
setdiff(gtdb$user_genome, QC$Sample) 
# merge
QC <- merge(QC, gtdb[,c('user_genome', 'species_classification')], by.x = 'Sample', by.y = 'user_genome')

# import kraken2/bracken outputs
bracken <-  read.table('input_data/public/long_reads/bracken_public_long_reads_b.tab', sep='\t', fill = TRUE, header = T, row.names = NULL)
bracken_added <-  read.table('input_data/public/long_reads/bracken_public_long_reads_added_b.tab', sep='\t', fill = TRUE, header = F, row.names = NULL)
colnames(bracken_added) <- colnames(bracken)
bracken <- merge(bracken, bracken_added, by = colnames(bracken), all = T)
bracken['Sample'] <- gsub('.bracken', '', bracken$filename)

# As some samples are classifies into multiple closely related species, I cannot just go for the bracken fraction to look for contamination. (I dont believe that these are contaminations, but rather differently assigned stretches of the same coherent genome)
# add info for which different genera were detected
bracken['Genus'] <- gsub('\\s+.*', '', bracken$name)
# add fraction per genus
bracken <- bracken %>% group_by(Sample, Genus) %>% 
  mutate(fraction_genus = sum(fraction_total_reads)) %>%
  # add random letter to each occurrence to make it unique before sumarising
  ungroup() %>% group_by(Sample) %>%
  mutate(fraction_genus_unique = paste0(fraction_genus, letters[dense_rank(Genus)]), 
         fraction_species_unique = paste0(fraction_total_reads, letters[dense_rank(name)])) %>%
  ungroup()

bracken_sum <- bracken %>% group_by(Sample) %>% 
  arrange(desc(fraction_total_reads)) %>% 
  summarise(genus_detected_bracken = paste0(unique(Genus), collapse = '_'), 
            n_genera_bracken = n_distinct(Genus),
            Most_abundant_genus_bracken = Genus[1], # is already sorted by abundance, so I can just take the first
            Most_abundant_genus_fraction_bracken = gsub('[[:alpha:]]','',fraction_genus_unique[1]),
            species_detected_bracken = paste0(unique(name), collapse = '_'),
            fraction_genus_bracken = gsub('[[:alpha:]]','', paste0(unique(fraction_genus_unique), collapse = '_')),
            species_fraction_bracken = gsub('[[:alpha:]]','', paste0(unique(fraction_species_unique), collapse = '_')))

# add to QC
setdiff(QC$Sample, bracken_sum$Sample)
QC <- merge(QC, bracken_sum, by = 'Sample')

# remove the following ones
# where main genus is not within order of Vibrionales or where > 10% or the contigs were not assigned as the most abundant genus (likely contaminations)
# remove all of the ones which start with SRR271* (these are all very short and most failed (but not all)), this is probably because they were not combined, but each run was submitted separately
# remove all with an average read depth < 5
# genome size smaller than 3.5 and larger than 7MB
# some strains were sequenced more than once, keep only repetition (ending with 'B'), remove first sequencing attempt
QC_1 <- QC[!(QC$Most_abundant_genus_bracken %in% c("Bacillus", "Enterococcus", "Staphylococcus","Klebsiella","Morganella", "Escherichia", "Acinetobacter", "Homo")| 
             as.numeric(QC$Most_abundant_genus_fraction < 0.9)| 
             grepl('^SRR271', QC$Sample) | 
             QC$Read_depth < 5 | 
             QC$Total_length_assembly < 3500000 | 
             QC$Total_length_assembly > 7000000),]
QC_1 <- QC_1[!is.na(QC_1$Sample),]
table(QC_1$species_classification)

remove <- QC[QC$Most_abundant_genus_bracken %in% c("Bacillus", "Enterococcus", "Staphylococcus","Klebsiella","Morganella", "Escherichia", "Acinetobacter", "Homo") | 
               as.numeric(QC$Most_abundant_genus_fraction < 0.9) | 
               grepl('^SRR271', QC$Sample) |
               QC$Read_depth < 5 | 
               QC$Total_length_assembly < 3500000 | 
               QC$Total_length_assembly > 7000000| 
               is.na(QC$Sample),]

# add whether a sample passed the QC or not and export the table for future reference
QC['quality_control'] <- ifelse(QC$Sample %in% QC_1$Sample, 'passed', 'failed')

write.table(QC[c('Sample', 'Median_read_length', 'N50_reads', 'Median_read_quality_nanostat', 'Read_depth', 'Contig_count', 'Total_length_assembly', 'GC_percent', 'species_classification', 'fraction_genus_bracken','quality_control')], 'input_data/public/long_reads/QC_evaluated_public_long_read.csv', quote = F, row.names = F, sep = ';')

# plot read quality and depth 
ggplot(QC_1, aes(y=Median_read_quality_nanostat, x=Read_depth)) + geom_point()
# the ones which have nanostat quality 0 were all pacbio sequenced (leave them in)
depth <- ggplot(QC_1, aes(x=Read_depth)) + geom_histogram(binwidth = 5)
table(QC_1$Contig_count) # some still have a lot of contigs, but they seem ok contamination wise, leave these in
n50_reads <- ggplot(QC_1, aes(x=N50_reads)) + geom_histogram()
assembly_length <- ggplot(QC_1, aes(x=Total_length_assembly)) + geom_histogram(aes(fill=species_classification))
# looks ok, all Vc around 4MB

# plot species overview
QC_1$species_classification <- ifelse(QC_1$species_classification == '', 'No GTDB classification', QC_1$species_classification)
QC_1$species_classification <- factor(QC_1$species_classification, levels = c(names(sort(table(QC_1$species_classification[!QC_1$species_classification == 'No GTDB classification']), decreasing =T)), 'No GTDB classification'))

species_breakdown_longread <- ggplot(QC_1, aes(x=species_classification)) + 
  geom_bar() + ylab('Nr. of long read assemblies') + xlab('GTDB classification') + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1))
species_breakdown_longread

pdf('output_figures/species_breakdown_longread.pdf', width = 5, height =3)
species_breakdown_longread 
dev.off()

# read assembly stats 
assembly_info <- read.table('input_data/public/long_reads/assembly_info_all_b.txt')
assembly_info_added <- read.table('input_data/public/long_reads/assembly_info_repeat_b.txt')
assembly_info <- merge(assembly_info, assembly_info_added, by = colnames(assembly_info_added), all=T)
colnames(assembly_info) <- c('seq_name','length','cov.','circ.','repeat','mult.','alt_group','graph_path','readfile_name')
assembly_info['sample'] <- gsub('((_1)|(A|B))*.fastq', '', assembly_info$readfile_name)

# check for which sample multiple input files were assembled
assembly_info_sample <- assembly_info[!duplicated(assembly_info$readfile_name),c('readfile_name', 'sample')]

# add how assembly is called
assembly_info['assembly_name'] <- gsub('.fastq(.gz)*', '', assembly_info$readfile_name)
assembly_info$assembly_name <- gsub('raw', '', assembly_info$assembly_name)
assembly_info$assembly_name <- gsub('_1$', '', assembly_info$assembly_name)
assembly_info$assembly_name <- gsub('_subreads$', '', assembly_info$assembly_name)
# '_1' and 'subreads' was not consistently taken out, add back 
assembly_info$assembly_name <- gsub('SRR17631716', 'SRR17631716_subreads', assembly_info$assembly_name)
assembly_info$assembly_name <- gsub('SRR18554162', 'SRR18554162_subreads', assembly_info$assembly_name)
assembly_info$assembly_name <- gsub('SRR23692898', 'SRR23692898_subreads', assembly_info$assembly_name)
assembly_info$assembly_name <- gsub('SRR12699728', 'SRR12699728_1', assembly_info$assembly_name)

# only keep the ones which passed QC
setdiff(QC_1$Sample,assembly_info$assembly_name)
assembly_info <- assembly_info[assembly_info$assembly_name %in% QC_1$Sample,]

# plot contig lengths
ggplot(assembly_info, aes(x=length)) +
  geom_histogram(binwidth = 100000, aes(fill = circ.))

ggplot(assembly_info[assembly_info$length > 800000,], aes(x=length)) +
  geom_histogram(binwidth = 100000, aes(fill = circ.))

# add information of how many large contigs (i.e. chromosomes) where assembled per strain
assembly_info <- assembly_info %>% group_by(sample) %>% mutate(Nr.of.large_contigs = sum(length>800000)) %>% ungroup()

# add information of how many large contigs were circularly (ie chromosomes) where assembled per strain
assembly_info <- assembly_info %>% group_by(sample) %>% mutate(Nr.of.large_contigs.circ = sum(length>800000 & circ.== 'Y')) %>% ungroup()
# remove the ones for which no large contig could be assembled
assembly_info_sel_large_contigs_circ <- assembly_info[assembly_info$Nr.of.large_contigs.circ > 0,]
assembly_info_sel_large_contigs_circ$Nr.of.large_contigs.circ <- factor(assembly_info_sel_large_contigs_circ$Nr.of.large_contigs.circ)

# check the ones for which more than 2 have been assembled (should not be the case) # --> none
assembly_info_sel_large_contigs_circ[as.numeric(as.character(assembly_info_sel_large_contigs_circ$Nr.of.large_contigs.circ)) >2,]

# add the information on total assembly length and length of cirucularly assembled large contigs
assembly_info <- assembly_info %>% group_by(sample) %>% mutate(total_assembly_length = sum(length), 
                                                             length_in_circ_larger = sum(length[length> 800000 & circ. == 'Y']), 
                                                             length_in_larger = sum(length[length> 800000])) %>% ungroup()

# plot only the ones where > 90% of assembly length is in long circular contigs (only here the number on chromosomes can be counted)
assembly_info_sel <- assembly_info
assembly_sel_plot <- assembly_info_sel[assembly_info_sel$length > 800000 & assembly_info_sel$total_assembly_length*0.9 <=assembly_info_sel$length_in_circ_larger,]

assembly_sel_plot$Nr.of.large_contigs.circ <- factor(assembly_sel_plot$Nr.of.large_contigs.circ, levels = c('2', '1'))
assembly_sel_plot['length_MB'] <- as.numeric(as.character(assembly_sel_plot$length))/1000000
sum_contig_length <- ggplot(assembly_sel_plot, aes(x=length_MB)) + theme_light() + theme(legend.position = 'bottom') +
  geom_histogram(binwidth = 0.1, aes(fill = factor(Nr.of.large_contigs.circ))) + scale_fill_manual(values = c('darkblue', 'lightblue')) +
  ylab('Nr. of long read assemblies') + xlab('Length of circular contig [MB]')
sum_contig_length

pdf('output_figures/sum_contig_length_public_long_reads.pdf', width = 3, height =3)
sum_contig_length 
dev.off()

length(unique(assembly_sel_plot$sample))
assembly_sel_plot$sample[assembly_sel_plot$Nr.of.large_contigs == '1']
# ERR9364039: Vc sample from Ghana, submitter Bernhard Nocht Institute for Tropical Medicine, 2022 (nothing stated about single chromosome)
  # https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC9257098&blobtype=pdf
  # strain-name: Iso02508
# ERR9364050: Vc sample from Ghana, submitter Bernhard Nocht Institute for Tropical Medicine, 2022 (nothing stated about single chromosome)
  # https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC9257098&blobtype=pdf
  # strain-name: Iso02520
# --> both strains mentioned as 1-contig assembled in suplementary table 1, but nothing mentioned in main text
# SRR10208201: 
#   - Vibrio cholerae El Tor strains isolated during cholera complications in Siberia and at the Far East of Russia in 90s
#   - Published 2020 here: https://www.sciencedirect.com/science/article/pii/S1567134819303223?via%3Dihub
#   - Collected in 1994 in Novosibirsk, strainnaame: 	I-1181 
#   - is described as being highly drug resistant and to have the entire VSP-II and its likely origin to be Bangladesh
# SRR24502230: V. natriensis strain, submitter Max Planck Institute for Terrestrial Microbiology, sequences published 04/2023, 
# check publication: https://www.biorxiv.org/content/10.1101/2023.05.26.541695v1.full.pdf $
# --> is a synthetically fused strain! no growth differences observed
# SRR24138850: This is the 5th large contig (of an assembly has two large circular contigs) is of a new unclassified (new?) acquatic vibrio species. I did not find a published study to this. 
# the published ion torrent / ONT hybrid assembly has three contigs of 4,055,171 / 2,244,890	and 37,555 bp

# check what lineage they are (only check the fused ones)
# I called SNV (and filtered) using medaka to assign the lineages
nr_snv <- read.table('input_data/public/long_reads/fused_medaka_hq_snv_counts.txt', sep = ':')
colnames(nr_snv) <- c('comparison', 'nr.snv')
# extract sample and ref from comparison column
nr_snv['sample'] <- gsub('(.*)(\\/.*)', '\\1', nr_snv$comparison)
nr_snv['ref'] <- gsub('(.*\\_)([[:alnum:]]*)(\\.fna.*)', '\\2', nr_snv$comparison)
# import which is which
sel <- read.table('input_data/sel_lineage_assignment_refs.csv', sep=';', header = T)
# check before merge
setdiff(nr_snv$ref, sel$Accession_merge)

nr_snv <- merge(nr_snv, sel[,c('Accession_merge','Lineage_group')], by.x='ref', by.y = 'Accession_merge', all.x=T)
nr_snv$nr.snv <- as.numeric(as.character(nr_snv$nr.snv))
nr_snv <- nr_snv %>% group_by(sample) %>% 
  mutate(rank_lineage = dense_rank(nr.snv)) # dense rank gives the same rank when there are ties
nr_snv[nr_snv$rank_lineage == '1',] # --> the ones isolated in Ghana are closest lineages IND1.1 whereas the ones isolated in Russia is closest to lineage IND.2

# subset to chromosomes 1 and 2 such that these can be blasted against each other
# subset assemblies in which 2 chromosomes where assembled and were 90% of the sequence length is in 
assembly_info_sel_2large_contigs <- assembly_info_sel[assembly_info_sel$Nr.of.large_contigs == '2' & assembly_info_sel$length > 800000 & assembly_info_sel$total_assembly_length*0.9 <=assembly_info_sel$length_in_larger,]
assembly_info_sel_2large_contigs <- assembly_info_sel_2large_contigs %>% group_by(sample) %>% mutate(chr = paste0('chr.',dense_rank(desc(length))))

# check chr length
quantile(assembly_info_sel_2large_contigs[assembly_info_sel_2large_contigs$chr == 'chr.1',]$length)
quantile(assembly_info_sel_2large_contigs[assembly_info_sel_2large_contigs$chr == 'chr.2',]$length)

# add the samplename to the contigname, to make them unique between the samples
assembly_info_sel_2large_contigs['contig_names'] <- paste0(assembly_info_sel_2large_contigs$assembly_name, '_',assembly_info_sel_2large_contigs$seq_name)

#write.table(assembly_info_sel_2large_contigs[assembly_info_sel_2large_contigs$chr == 'chr.1',]$contig_names, row.names = F, quote = F, 'input_data/public/long_reads//chr1_names.txt')
#write.table(assembly_info_sel_2large_contigs[assembly_info_sel_2large_contigs$chr == 'chr.2',]$contig_names, row.names = F, quote = F, 'input_data/public/long_reads//chr2_names.txt')
#write.table(unique(assembly_info_sel_2large_contigs$assembly_name), row.names = F, quote = F, 'input_data/public/long_reads/assemblies_2chr.txt')
# these were exported and the two chromosomes of each strains were compared against each pther using blastn


# evaluate the blast output
blastout_q_full <- read.table('input_data/public/long_reads/all_chr_blastout.tab')
# add columnnames
colnames(blastout_q_full) <-  c('qseqid', 'sseqid', 'bitscore', 'pident', 'nident', 'mismatch', 'length', 'qcovs', 'qlen', 'qstart', 'qend',  'slen', 'sstart', 'send', 'qseq', 'sseq')
# only include the ones which passed QC
blastout_q_full <- as.data.frame(blastout_q_full)
blastout_q_full['Sample'] <- gsub('_contig.*','',blastout_q_full$qseqid)
# check the ones which were not in blast although they passed QC
check <- assembly_info[assembly_info$assembly_name %in% setdiff(QC_1$Sample, c(blastout_q_full$Sample, remove$Sample)),]
check <- check[!duplicated(check$assembly_name),]
# the ones with 0 or 1 assembled large circular contig (should not be in)

# save as numeric
blastout_q_full$pident <- as.numeric(as.character(blastout_q_full$pident))
blastout_q_full$length<- as.numeric(as.character(blastout_q_full$length))
# only consider samples which passed the QC
blastout_q_full <- blastout_q_full[blastout_q_full$Sample %in% QC_1$Sample,]
# only consider sequences with 90% identity or more
blastout_q <- blastout_q_full[blastout_q_full$pident >= 90,]

# plot the number of blasthits their length found
blastout_q <- blastout_q %>% group_by(qseqid) %>% mutate(n_matches = n())
ggplot(blastout_q, aes(x=n_matches, y=length, fill=qseqid)) +
  geom_boxplot() + theme(legend.position = 'none') +
  facet_zoom(xlim = c(0, 300))

# add species information
blastout_q <- merge(blastout_q, QC_1[,c("Sample", "species_classification")], by = 'Sample')
# only keep species which for which more than 2 genomes were available
blastout_q <- blastout_q[blastout_q$species_classification %in% names(table(QC_1$species_classification)[table(QC_1$species_classification) >2]),]

# summarise the number of blasthits per sample, their mean and total length (is total overlap between chr1 and chr2)
blastout_q_unique <- blastout_q %>% group_by(qseqid, species_classification) %>% summarise(n_matches = n(), 
                                                                 mean_length = mean(length), 
                                                                 total_length = sum(length))
# plot the number of matches
p_n_matches <- ggplot(blastout_q_unique, aes(x=n_matches)) + 
  geom_histogram(binwidth = 30) +
  facet_zoom(xlim = c(0, 1200)) + 
  xlab('Nr. of blasthits per sample')
# convert to numeric
blastout_q_unique$mean_length <- as.numeric(as.character(blastout_q_unique$mean_length))
blastout_q_unique$total_length <- as.numeric(as.character(blastout_q_unique$total_length))
# plot mean and total lengths
p_mean_length <- ggplot(blastout_q_unique, aes(x=n_matches, y=mean_length)) + # color ba species when kraken/bracken done
  geom_point() +
  facet_zoom(xlim = c(0, 1200)) + 
  xlab('Nr. of blasthits per sample') 
p_total_length <- ggplot(blastout_q_unique, aes(x=n_matches, y=total_length)) + # color ba species when kraken/bracken done
  geom_point(aes(colour  = species_classification)) +
  facet_zoom(xlim = c(0, 650), ylim = c(2500,200000)) + 
  xlab('Nr. of blasthits per sample') 
# --> not sure how much it actually makes sense to summarise all matches per sample, plot each match individually and color by species

# as there are many more short than long hits, I will split these into two different plots (to get two different Y axis)
p_length_per_hit_short <- ggplot(blastout_q, aes(x=length)) + 
  geom_histogram(bins = 200, aes(fill  = species_classification)) +
  #facet_zoom(xlim = c(5000, 30000), ylim = c(1,30)) + 
  facet_zoom(xlim = c(0, 5000)) + 
  xlab('length [bp]')
p_length_per_hit_large <- ggplot(blastout_q, aes(x=length)) + 
  geom_histogram(bins = 50, aes(fill  = species_classification)) +
  facet_zoom(xlim = c(5000, 20000), ylim = c(1,30)) + 
  #facet_zoom(xlim = c(0, 5000)) + 
  xlab('length [bp]')
p_length_per_hit_short / p_length_per_hit_large

# try splittng up per species
# shorten genus names for plot
blastout_q$species_classification <- gsub('Vibrio', 'V.', blastout_q$species_classification)
blastout_q$species_classification <- gsub('Photobacterium', 'P.', blastout_q$species_classification)

blastout_q_sum <- blastout_q %>% mutate(length_bin=cut(length, breaks=seq(0,15000, by=100))) %>% # bin the blasthits into the same defined bin for all (histogram same bin for all species)
  group_by(length_bin) %>% mutate(avg_bin_position = mean(length)) %>% # add mean per length 
  ungroup() %>%
  group_by(species_classification) %>% mutate(n_samples_per_spp = n_distinct(Sample)) %>%
  group_by(length_bin, species_classification, n_samples_per_spp, avg_bin_position) %>% 
  summarise(n_hits_per_bin = n()) %>% mutate(n_hits_per_bin_norm = n_hits_per_bin / n_samples_per_spp) # normalise by the number of samples (as there are many more Vc samples than samples of other species). 
# order species, such that Vc is plotted first
blastout_q_sum$species_classification <-factor(blastout_q_sum$species_classification, levels = c('V. cholerae', 'V. parahaemolyticus', 'V. anguillarum', 'P. damselae', 'V. campbellii', 'V. natriegens'))

# split into short (smaller than 3000bp) and long (longer than 3000bp) hots
blasthits_norm_short <- ggplot(blastout_q_sum[blastout_q_sum$avg_bin_position < 3000,], aes(x=avg_bin_position, y= n_hits_per_bin_norm)) + geom_bar(stat="identity") +
  facet_grid(species_classification~., scales = 'free_y') + 
  xlab('Length of BLASTn hit [bp] (binsize=100)') + ylab('Nr. of BLASTn hits between the two chromosomes /\nNr. of assemblies') +
  theme_light()
blasthits_norm_long <- ggplot(blastout_q_sum[blastout_q_sum$avg_bin_position >= 3000,], aes(x=avg_bin_position, y= n_hits_per_bin_norm)) + geom_bar(stat="identity") +
  facet_grid(species_classification~., scales = 'free_y', drop = FALSE) + 
  xlab('Length of BLASTn hit [bp] (binsize=100)') + ylab('Nr. of BLASTn hits between the two chromosomes /\nNr. of assemblies') +
  theme_light()
blasthits_norm_short | blasthits_norm_long

# export
pdf('output_figures/chr_blast_public_long_reads.pdf', width =6, height =7)
blasthits_norm_short | blasthits_norm_long
dev.off()


# Next, I will look more closely at the hits for Vc
# subset for Vc hits
blastout_q_Vc <- blastout_q[blastout_q$species_classification == 'V. cholerae',]
# Compare all the Vc hits > 10000 using fastANI (to check if they are all the same sequence)
write_out <- blastout_q_Vc[blastout_q_Vc$length > 10000,]

# export
#for (i in 1:nrow(write_out)){
#  print(i)
#  write.fasta(write_out[i,'sseq'], write_out[i,'sseqid'], paste0('./Vc_over_10000_chr_overlap_',i,'.fasta'), open = "w", nbchar = 60, as.string = FALSE)
#}

# subset for Vc sample
blastout_q_unique_Vc <- blastout_q_unique[blastout_q_unique$species_classification == 'Vibrio cholerae',]

# import the fastANI outputs
# I ran fastANI on all the Vc hits > 10'000 as well as a reference for VSP-I and for HS1 (as I by comparing the locations on the genome I expect these to match)
ANI_out <- read.table('input_data/public/long_reads/fastANI_out_incl_ref.txt')
colnames(ANI_out) <- c('query_genome', 'reference_genome', 'ANI_value', 'count_bidirectional_fragment_mappings', 'total_query_fragments')
ANI_out$query_genome <- gsub('.*\\/', '', ANI_out$query_genome)
ANI_out$reference_genome <- gsub('.*\\/', '', ANI_out$reference_genome)

# plot heatmap to visualise cluster
# shorten name, only use sample labels as names
ANI_out$query_genome <- gsub('Vc_over_10000_chr_overlap_', 'V.cholerae chr. overlap ', ANI_out$query_genome)
ANI_out$query_genome <- gsub('.fasta', '', ANI_out$query_genome)
ANI_out$reference_genome <- gsub('Vc_over_10000_chr_overlap_', 'V.cholerae chr. overlap ', ANI_out$reference_genome)
ANI_out$reference_genome <- gsub('.fasta', '', ANI_out$reference_genome)
# use label HS1
ANI_out$reference_genome <- gsub('overlap_134Vc041_134Vc042', '12kb HS1', ANI_out$reference_genome)
ANI_out$query_genome <- gsub('overlap_134Vc041_134Vc042', '12kb HS1', ANI_out$query_genome)

#order according to values, to better visualise clusters
ANI_out_1 <- ANI_out[ANI_out$reference_genome == 'V.cholerae chr. overlap 1',]
ANI_out_2 <- ANI_out[ANI_out$reference_genome == 'V.cholerae chr. overlap 9',]
order_plot <- c(ANI_out_1[order(ANI_out_1$ANI_value),'query_genome'], ANI_out_2[order(ANI_out_2$ANI_value),'query_genome'])

# convert to factors
ANI_out$reference_genome <- factor(ANI_out$reference_genome, levels = order_plot)
ANI_out$query_genome <- factor(ANI_out$query_genome, levels = order_plot)

# plot heatmap
ANI_overlaps_10000_Vc <- ggplot(ANI_out, aes(query_genome, reference_genome, fill= ANI_value)) + 
  geom_tile() + theme_light() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_fill_gradient(low = 'darkblue', high = 'orange') +xlab('') + ylab('') + theme(legend.position = 'bottom')
# two clusters are formed, take one of each (Vc_over_10000_chr_overlap_1.fasta and Vc_over_10000_chr_overlap_9.fasta) and annotate these with bakta
ANI_overlaps_10000_Vc
# no ANI values are reported if the are below 70% (as they are not regarded reliable), this is why these values are missing

pdf('output_figures/ANI_overlaps_10000_Vc.pdf', width =7.5, height =7)
ANI_overlaps_10000_Vc
dev.off()


# I check how many copies of HS1 I find in each of the Vc long read genomes using blastn
blastout_12kbHR <- read.table('input_data/public/long_reads/blastout_HR_Vc_longread_public.tab')
colnames(blastout_12kbHR) <-  c('qseqid', 'sseqid', 'bitscore', 'pident', 'nident', 'mismatch', 'length', 'qcovs', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send')

ggplot(blastout_12kbHR, aes(x=pident)) +
  geom_histogram()
ggplot(blastout_12kbHR, aes(x=qcovs)) +
  geom_histogram()
# use 95 and 95 as threshold for pident and qcovs
blastout_12kbHR_Sel <- blastout_12kbHR[blastout_12kbHR$pident >=95 & (blastout_12kbHR$length >= blastout_12kbHR$qlen*0.90),]
# count per assembly
blastout_12kbHR_Sel['sample'] <- gsub('_contig_.*', '', blastout_12kbHR_Sel$sseqid)
blastout_12kbHR_Sel_sum <- blastout_12kbHR_Sel %>% group_by(sample) %>% summarise(n_hits = n(), n_sseqs = n_distinct(sseqid), sseq_lengths = paste0(unique(slen), collapse = '_'))

# check how many HS1 are present in the public assemblies which assembled to one chromosome 
blastout_12kbHR_Sel_sum[blastout_12kbHR_Sel_sum$sample %in% c("ERR9364039","ERR9364050","SRR10208201"),] # --> they all carry 2

# check if orientation of HS1 is always the same
blastout_12kbHR_Sel['HS1_orientation'] <- ifelse(blastout_12kbHR_Sel$sstart < blastout_12kbHR_Sel$send, 'F', 'R')
blastout_12kbHR_Sel[blastout_12kbHR_Sel$sample %in% c("ERR9364039","ERR9364050","SRR10208201"),]
