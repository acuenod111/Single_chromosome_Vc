library('dplyr')
library('tidyr')
library('ggplot2')
library('ggforce')

# import list of references (chosen form Monir et al. NatCom 2023, supplementary tables 1 and 3)
sel <- read.csv2('input_data/sel_lineage_assignment_refs.csv',  sep=';')

# I previously tested if I could use ANI, but this gave be very similar and ambiguous values to multiple references, so I decided to include high quality SNVs only. I called SNV using medaka and filtered for a min quality score of 40 
nr_snv <- read.table('input_data/medaka_hq_snv_counts.txt', sep = ':')
colnames(nr_snv) <- c('comparison', 'nr.snv')
# extract sample and ref from comparison column
nr_snv['sample'] <- gsub('(\\d{3}Vc\\d{2})(\\/.*)', '\\1', nr_snv$comparison)
nr_snv['ref'] <- gsub('(.*\\_)([[:alnum:]]*)(\\.fna.*)', '\\2', nr_snv$comparison)
# check before merge
setdiff(nr_snv$ref, sel$Accession_merge)

# merge to reference strain list 
nr_snv <- merge(nr_snv, sel[,c('Accession_merge','Lineage_group')], by.x='ref', by.y = 'Accession_merge', all.x=T)

# add rank (give rank 1 to the reference where the fewest SNV were called)
nr_snv <- nr_snv %>% group_by(sample) %>% 
  mutate(rank_lineage = dense_rank(nr.snv)) # dense rank gives the same rank when there are ties

# convert to binary (assign the lineage with the lowest number of SNV detected)
nr_snv['rank_lineage_binary'] <- ifelse(nr_snv$rank_lineage == '1', '1', 'lower_rank')

# some lineages with a lower rank have fewer SNV than lineages with rank 1 (can be true for closely related lineages, still something to keep in mind)
check_lineage <- ggplot(nr_snv, aes(x=rank_lineage_binary, y=nr.snv)) +
  geom_jitter(width = 0.2) +
  facet_zoom(ylim = c(0, 400))

# not 100% conclusive
table(nr_snv[nr_snv$rank_lineage_binary == '1',]$Lineage_group) # I would have expected B2 and B1.2, its unexpected (but not impossible) that it gives be IND.3

# check the ones which were assigned an IND sublineage and how many more SNV were called to the BD1.2 reference
check_ind <- nr_snv[nr_snv$sample %in% nr_snv[nr_snv$rank_lineage_binary == '1' & nr_snv$Lineage_group %in% c('IND1.2','IND1.3'),]$sample,]
check_ind_w <- check_ind %>% 
  tidyr::pivot_wider(., id_cols = sample, names_from = Lineage_group, values_from = nr.snv)
check_ind_w['diff_BD1.2_min'] <- check_ind_w$BD1.2 - do.call(pmin, check_ind_w[,c("IND1.1","IND1.3","IND1.2","BD2","BD1.1","BD1.2","IND2","BD1","IND1","AS1","T12","AS2","T13","LAT3")])
check_ind_w['diff_BD1.1_min'] <- check_ind_w$BD1.1 - do.call(pmin, check_ind_w[,c("IND1.1","IND1.3","IND1.2","BD2","BD1.1","BD1.2","IND2","BD1","IND1","AS1","T12","AS2","T13","LAT3")])

quantile(check_ind_w$diff_BD1.2_min)
quantile(check_ind_w$diff_BD1.1_min)
# --> For the BD1 lineages less SNV more high quality SNV are called, so I think these are actually closer to IND3 / IND2

# check if there are any ties with rank 1 assignments
check <- nr_snv[nr_snv$rank_lineage == '1',]
check[duplicated(check$sample),] # none of the rank 1 are duplicated

# check how much the difference between first and second rank assignment
nr_snv_check <- nr_snv %>% 
  # if rank 2 are ties arndomply choose one
  group_by(sample, rank_lineage) %>%
  slice(1) %>%
  ungroup() %>% group_by(sample) %>%
  summarise(diff_rank1_2 = nr.snv[rank_lineage == '2']  - nr.snv[rank_lineage == '1'])

nr_snv_check <- merge(nr_snv_check, nr_snv[nr_snv$rank_lineage == '1',c('sample', 'Lineage_group')], all.x=T, all.y=F)
ggplot(nr_snv_check) + geom_histogram(aes(x=diff_rank1_2, fill = Lineage_group), binwidth = 5) + facet_zoom(xlim = c(0,250))
# I thought the BD1 references might be incorrect, I doublechecked and it says that this should be BD1. 
# I could not find a tree with BD1 strains and IND strains, it just states that the BD1 strains are closer to the IND strains than the BD2 strains, but were collected from Bangladesh.
# Our samples were collected before BD1.2 was defined as such (outbreak 2022))
# I think it makes sense to call the lineage to whiich the smallest number of SNV were called. 
# Export the lineage with rank 1
write.table(nr_snv[nr_snv$rank_lineage_binary == '1',c('sample','ref','nr.snv', 'Lineage_group', 'rank_lineage')], 'input_data/lineage_assignemnet_from_medaka.txt', quote = F, row.names = F)

