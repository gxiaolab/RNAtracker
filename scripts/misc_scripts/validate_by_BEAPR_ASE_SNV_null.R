rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(ComplexUpset)

args = commandArgs(trailingOnly=TRUE)
CELL_LINE = args[1]
WGS_LABEL = args[2]
EFFECT_SIZE_CUTOFF = 0.1

metadata = fread('../prepare_metadata', header = FALSE)
classified_genes = fread('../../../HIER_MODEL_V2/snv_types/ASRS_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed')
all_genes = fread('../../../HIER_MODEL_V2/bed_files/all_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed')


cell_line_spec_ALL = filter(all_genes, V8 == CELL_LINE) %>% unique()
cell_line_spec_ALL = cell_line_spec_ALL %>% mutate(id = paste0(V1, ':', V3, ':', V6, ':', V7))

read_dat = function(WGS_LABEL){
  in_fn = paste0('../REAL_DATA_RECOVERED_SNVS_10_2/', WGS_LABEL, '.ALL_TIMEPOINTS.RECOVER_SNVs_INCLUDE_ZERO.RData')
  load(in_fn)
  return(all_timepoints_complete)
}

ACTD_DAT = read_dat(WGS_LABEL)

filter(ACTD_DAT, gene_name %in% cell_line_spec_ALL$V4) %>% select(id) %>% unique()
filter(ACTD_DAT, gene_name %in% cell_line_spec_ALL$V4) %>% select(gene_name) %>% unique()

## compare with Bru direction
bru_res = fread('../../../HIER_MODEL_V2/results/full_ALL_GENES_counts_data.txt')
spec_bru_dat =  bru_res %>% filter(cell_line %in% CELL_LINE) #%>% filter(gene_state %in% c('state1', 'state2')) 

nonASE_0h_SNV = filter(ACTD_DAT, p_adj_0h > 0.05) #%>% mutate(AR_0h = MODEL_COUNT_0h/total_counts_0h) %>% mutate(delta_AR_0h = abs(AR_0h - 0.5))
# for ASE timepoints, the allelic imbalance should be the same across replicates
nonASE_0h_SNV %>% filter(gene_name %in% cell_line_spec_ALL$V4) %>% select(id) %>% unique()
nonASE_0h_SNV %>% filter(gene_name %in% cell_line_spec_ALL$V4) %>% select(gene_name) %>% unique()

ASE_0h_SNV = filter(ACTD_DAT, p_adj_0h < 0.05) %>% filter(gene_name %in% cell_line_spec_ALL$V4) %>% select(id) %>% unique()
filter(ACTD_DAT, p_adj_0h < 0.05) %>% filter(gene_name %in% cell_line_spec_ALL$V4) %>% select(gene_name) %>% unique()


check_replicate_imbalance_direction = function(timepoint){
  spec_dat = nonASE_0h_SNV %>% select(id, gene_name, region, replicate, contains(timepoint))
  names(spec_dat) = gsub(paste0('_', timepoint), '', names(spec_dat))
  spec_dat$AR = spec_dat$MODEL_COUNT/spec_dat$total_counts
  spec_dat$IMBALANCE_direction = sign(spec_dat$AR - 0.5)
  spec_dat$delta_AR = abs(spec_dat$AR - 0.5)

  ## get testable snvs and testable genes
  testable_dat = filter(spec_dat, id %in% cell_line_spec_ALL$id)
  
  consistent_replicates = testable_dat %>% filter(p_adj < 0.05) %>% filter(gene_name %in% cell_line_spec_ALL$V4) %>% filter(!is.na(IMBALANCE_direction)) %>% select(id, IMBALANCE_direction) %>% unique() %>% count(id) %>% arrange(desc(n)) %>% filter(n==1)
  validated_snvs = filter(testable_dat, id %in% consistent_replicates$id) %>% filter(gene_name %in% cell_line_spec_ALL$V4)
  # now apply effect size cutoff: at least one replicate should pass
  pass_snvs = validated_snvs %>% filter(delta_AR >= EFFECT_SIZE_CUTOFF)
  validated_snvs = filter(validated_snvs, id %in% pass_snvs$id)

  ## annotate snv with the direction of allelic imbalance
  compare_res = inner_join(spec_bru_dat %>% select(id, rep1_AR_0h, rep1_AR_2h, rep1_AR_6h, rep2_AR_0h, rep2_AR_2h, rep2_AR_6h), validated_snvs %>% select(id, IMBALANCE_direction) %>% unique(), by = 'id')
  stable_ref_2h = compare_res %>% filter(rep1_AR_2h > 0.5 & rep2_AR_2h > 0.5 & IMBALANCE_direction == 1)
  stable_ref_6h = compare_res %>% filter(rep1_AR_6h > 0.5 & rep2_AR_6h > 0.5 & IMBALANCE_direction == 1)
  stable_alt_2h = compare_res %>% filter(rep1_AR_2h < 0.5 & rep2_AR_2h < 0.5 & IMBALANCE_direction == -1)
  stable_alt_6h = compare_res %>% filter(rep1_AR_6h < 0.5 & rep2_AR_6h < 0.5 & IMBALANCE_direction == -1)
  consistent_direction_snvs = unique(c(stable_ref_2h$id, stable_ref_6h$id, stable_alt_2h$id, stable_alt_6h$id))
  testable_dat$snv_validation = testable_dat$id %in% consistent_direction_snvs

  total_snv_counts = left_join(testable_dat, cell_line_spec_ALL %>% select(id, V5), by = 'id', relationship = "many-to-many") %>% group_by(V5) %>% count(V5) %>% rename(total_testable_snvs = n) # many-to-many relationship because some snvs are in more than one gene 
  testable_snv_counts = left_join(testable_dat, cell_line_spec_ALL %>% select(id, V5), by = 'id', relationship = "many-to-many") %>% group_by(V5) %>% count(snv_validation) %>% rename(total_status_snvs = n)
  summary_dat = left_join(total_snv_counts, testable_snv_counts, by = 'V5') %>% mutate(prop_snvs = total_status_snvs/total_testable_snvs)
  out = list(consistent_direction_snvs, summary_dat)
  return(out)
}

validate_2h = check_replicate_imbalance_direction('2h')
validate_8h = check_replicate_imbalance_direction('8h')
validate_24h = check_replicate_imbalance_direction('24h')


## get full information 
full_info_2h = nonASE_0h_SNV %>% filter(id %in% validate_2h[[1]]) %>% select(id, gene_name, region, cell_line, replicate, contains('0h'), contains('2h')) %>% select(-contains('timepoint')) %>% mutate(timepoint = '2h')
names(full_info_2h) = gsub('2h', 'later_timepoint', names(full_info_2h))
full_info_2h$spec_label = 'nonASE_0h_2h'

full_info_8h = nonASE_0h_SNV %>% filter(id %in% validate_8h[[1]]) %>% select(id, gene_name, region, cell_line, replicate, contains('0h'), contains('8h')) %>% select(-contains('timepoint')) %>% mutate(timepoint = '8h')
names(full_info_8h) = gsub('8h', 'later_timepoint', names(full_info_8h))
full_info_8h$spec_label = 'nonASE_0h_8h'

full_info_24h = nonASE_0h_SNV %>% filter(id %in% validate_24h[[1]]) %>% select(id, gene_name, region, cell_line, replicate, contains('0h'), contains('24h')) %>% select(-contains('timepoint')) %>% mutate(timepoint = '24h')
names(full_info_24h) = gsub('24h', 'later_timepoint', names(full_info_24h))
full_info_24h$spec_label = 'nonASE_0h_24h'

all_method1_res = rbind(full_info_2h, full_info_8h, full_info_24h)



##### NOW TEST THE SNVs that are ASE at 0h
ASE_SNV_0h = filter(ACTD_DAT, p_adj_0h < 0.05) %>% filter(gene_name %in% cell_line_spec_ALL$V4)
ASE_SNV_0h$AR_0h = ASE_SNV_0h$MODEL_COUNT_0h/ASE_SNV_0h$total_counts_0h
ASE_SNV_0h$delta_AR_0h = ASE_SNV_0h$AR_0h - 0.5
check_ASE_0h_later_timepoints = function(timepoint){
  spec_dat = ASE_SNV_0h %>% select(id, gene_name, region, replicate, AR_0h, delta_AR_0h, contains(timepoint))
  names(spec_dat) = gsub(paste0('_', timepoint), '', names(spec_dat))
  spec_dat$AR = spec_dat$MODEL_COUNT/spec_dat$total_counts
  spec_dat$delta_AR = spec_dat$AR - 0.5

  ## get testable snvs and testable genes
  testable_dat = filter(spec_dat, id %in% cell_line_spec_ALL$id)

  spec_dat = filter(spec_dat, p_adj < 0.05) %>% filter(abs(delta_AR) > abs(delta_AR_0h)) %>% filter(abs(delta_AR) >= EFFECT_SIZE_CUTOFF)


  spec_dat$IMBALANCE_direction = sign(spec_dat$delta_AR)
  
  consistent_replicates = spec_dat %>% filter(gene_name %in% cell_line_spec_ALL$V4) %>% filter(!is.na(IMBALANCE_direction)) %>% select(id, IMBALANCE_direction) %>% unique() %>% count(id) %>% arrange(desc(n)) %>% filter(n==1)
  validated_snvs = filter(spec_dat, id %in% consistent_replicates$id) %>% filter(gene_name %in% cell_line_spec_ALL$V4)

  # now apply effect size cutoff: at least one replicate should pass
  pass_snvs = validated_snvs %>% filter(abs(delta_AR) >= EFFECT_SIZE_CUTOFF)
  validated_snvs = filter(validated_snvs, id %in% pass_snvs$id)

  ## annotate snv with the direction of allelic imbalance
  compare_res = inner_join(spec_bru_dat %>% select(id, rep1_AR_0h, rep1_AR_2h, rep1_AR_6h, rep2_AR_0h, rep2_AR_2h, rep2_AR_6h), validated_snvs %>% select(id, IMBALANCE_direction) %>% unique(), by = 'id')
  stable_ref_2h = compare_res %>% filter(rep1_AR_2h > 0.5 & rep2_AR_2h > 0.5 & IMBALANCE_direction == 1)
  stable_ref_6h = compare_res %>% filter(rep1_AR_6h > 0.5 & rep2_AR_6h > 0.5 & IMBALANCE_direction == 1)
  stable_alt_2h = compare_res %>% filter(rep1_AR_2h < 0.5 & rep2_AR_2h < 0.5 & IMBALANCE_direction == -1)
  stable_alt_6h = compare_res %>% filter(rep1_AR_6h < 0.5 & rep2_AR_6h < 0.5 & IMBALANCE_direction == -1)
  consistent_direction_snvs = unique(c(stable_ref_2h$id, stable_ref_6h$id, stable_alt_2h$id, stable_alt_6h$id))

  testable_dat$snv_validation = testable_dat$id %in% consistent_direction_snvs

  total_snv_counts = left_join(testable_dat, cell_line_spec_ALL %>% select(id, V5), by = 'id', relationship = "many-to-many") %>% group_by(V5) %>% count(V5) %>% rename(total_testable_snvs = n) # many-to-many relationship because some snvs are in more than one gene 
  testable_snv_counts = left_join(testable_dat, cell_line_spec_ALL %>% select(id, V5), by = 'id', relationship = "many-to-many") %>% group_by(V5) %>% count(snv_validation) %>% rename(total_status_snvs = n)
  summary_dat = left_join(total_snv_counts, testable_snv_counts, by = 'V5') %>% mutate(prop_snvs = total_status_snvs/total_testable_snvs)
  out = list(consistent_direction_snvs, summary_dat)
  return(out)
}

validate_2h_ASE = check_ASE_0h_later_timepoints('2h')
validate_8h_ASE = check_ASE_0h_later_timepoints('8h')
validate_24h_ASE = check_ASE_0h_later_timepoints('24h')


## get full information 
full_ASE_info_2h = ASE_SNV_0h %>% filter(id %in% validate_2h_ASE[[1]]) %>% select(id, gene_name, region, cell_line, replicate, contains('0h'), contains('2h')) %>% select(-contains('timepoint')) %>% mutate(timepoint = '2h')
names(full_ASE_info_2h) = gsub('2h', 'later_timepoint', names(full_ASE_info_2h))
full_ASE_info_2h$spec_label = 'ASE_0h_2h'

full_ASE_info_8h = ASE_SNV_0h %>% filter(id %in% validate_8h_ASE[[1]]) %>% select(id, gene_name, region, cell_line, replicate, contains('0h'), contains('8h')) %>% select(-contains('timepoint')) %>% mutate(timepoint = '8h')
names(full_ASE_info_8h) = gsub('8h', 'later_timepoint', names(full_ASE_info_8h))
full_ASE_info_8h$spec_label = 'ASE_0h_8h'

full_ASE_info_24h = ASE_SNV_0h %>% filter(id %in% validate_24h_ASE[[1]]) %>% select(id, gene_name, region, cell_line, replicate, contains('0h'), contains('24h')) %>% select(-contains('timepoint')) %>% mutate(timepoint = '24h')
names(full_ASE_info_24h) = gsub('24h', 'later_timepoint', names(full_ASE_info_24h))
full_ASE_info_24h$spec_label = 'ASE_0h_24h'

all_method2_res = rbind(full_ASE_info_2h, full_ASE_info_8h, full_ASE_info_24h)


## OUTPUT FULL VALIDATED SNV INFO FOR SUPPLEMENTARY TABLE
supp_out = rbind(all_method1_res , all_method2_res %>% select(-c(AR_0h, delta_AR_0h)))
supp_out = left_join(supp_out, spec_bru_dat %>% select(gene_name, gene_state, id) %>% unique(), by = c('id', 'gene_name')) %>% filter(gene_state %in% c('state1', 'state2'))

out_supp_fn =  paste0('../SUPPLEMENTARY_TABLES/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV_ALL_STATE_GENES.txt')
supp_out$cell_line = CELL_LINE
write.table(supp_out, file = out_supp_fn, sep = '\t', row.names = FALSE, quote = F)



## now make a data frame annotating all the validated snvs
df1 = data.frame(id = validate_2h[[1]], spec_label = 'nonASE_0h_validate_2h', general_label = 'nonASE_0h')
df2 = data.frame(id = validate_8h[[1]], spec_label = 'nonASE_0h_validate_8h', general_label = 'nonASE_0h')
df3 = data.frame(id = validate_24h[[1]], spec_label = 'nonASE_0h_validate_24h', general_label = 'nonASE_0h')
df4 = data.frame(id = validate_2h_ASE[[1]], spec_label = 'ASE_0h_validate_2h', general_label = 'ASE_0h')
df5 = data.frame(id = validate_8h_ASE[[1]], spec_label = 'ASE_0h_validate_8h', general_label = 'ASE_0h')
df6 = data.frame(id = validate_24h_ASE[[1]], spec_label = 'ASE_0h_validate_24h', general_label = 'ASE_0h')
all_df = rbind(df1, df2, df3, df4, df5, df6)

## annotate with gene
all_dat = left_join(all_df, spec_bru_dat %>% select(gene_name, gene_state, id) %>% unique(), by = 'id') # multiple matches are due to same snv being in multiple snvs



## get ALL testable snvs and genes 
ALL_TESTABLE_DAT = filter(ACTD_DAT, id %in% cell_line_spec_ALL$id) %>% filter(gene_name %in% cell_line_spec_ALL$V4)
ALL_TESTABLE_DAT = left_join(ALL_TESTABLE_DAT, spec_bru_dat %>% select(gene_name, gene_state, id) %>% unique(), by = c('gene_name', 'id')) 

all_results = left_join(ALL_TESTABLE_DAT %>% select(id, gene_name, gene_state) %>% unique(), all_dat, by = c('id', 'gene_name', 'gene_state'), relationship = "many-to-many")
all_results$spec_label[is.na(all_results$spec_label)] = 'not validated'
all_results$general_label[is.na(all_results$general_label)] = 'not validated'
all_results$validation_status = !(all_results$general_label == 'not validated')

out_data =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV_ALL_STATE_GENES.RData')
save(all_results, file = out_data)


