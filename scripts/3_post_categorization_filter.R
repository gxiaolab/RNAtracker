rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)

#' This function is from
#' [StackOverflow post 5577221](https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file)
loadRData <- function(file_name){
  load(file_name)
  get(ls()[ls() != "file_name"])
}



args = commandArgs(trailingOnly=TRUE)
CELL_LINE = args[1]
INPUT_DIR = args[2] # ../data
RESULTS_DIR = args[3] # ../results


MAX_PROB_CUTOFF =0.95
STATE_QUANT_CUTOFF = 0.333333333333333
PASS_PROP_THRESHOLD = 0.50 # at ASE timepoints (i.e. 2h and 6h for state1 asRS genes, all three timepoints for asRT genes), what proportion of SNVs must have consistent allelic imbalance direction across replicates?



get_real_data = function(CELL_LINE, WHICH_MODEL, MAX_PROB_CUTOFF, STATE_QUANT_CUTOFF){
  output_fn = paste0(RESULTS_DIR, '/', CELL_LINE,  '_', WHICH_MODEL, '_model.est.Rdata')
  load(output_fn)
  all_df = final_step_post
  cutoff_df = all_df %>% group_by(predicted_state) %>% summarise(quantile = scales::percent(c(STATE_QUANT_CUTOFF)),quant_cutoff = quantile(max_probability, c(STATE_QUANT_CUTOFF)))
  processed_all_df = left_join(all_df, cutoff_df, by = 'predicted_state')
  passing_genes = processed_all_df %>% filter(max_probability >= MAX_PROB_CUTOFF | max_probability >= quant_cutoff)
  processed_all_df$FILTER_PASS = processed_all_df$gene_id %in% passing_genes$gene_id
  return(processed_all_df)
}


get_snv_count = function(CELL_LINE){
  timepoints = c('zero_hr', 'two_hr', 'six_hr')
   input_fn = paste0(INPUT_DIR, '/', CELL_LINE, '.snv_info.txt.gz')
   all_data = fread(input_fn)
  # now load the classified gene dat
  ASE_0h_out = get_real_data(CELL_LINE, 'ASE_0h', MAX_PROB_CUTOFF, STATE_QUANT_CUTOFF)
  nonASE_0h_out = get_real_data(CELL_LINE, 'nonASE_0h', MAX_PROB_CUTOFF, STATE_QUANT_CUTOFF)
  all_results_dat = rbind(ASE_0h_out %>% select(-contains('comp')), nonASE_0h_out %>% select(-contains('comp')))
  all_results_dat$predicted_state = paste0('state', all_results_dat $predicted_state)
  
  combine_dat = left_join(all_results_dat %>% rename(gene_name = gene_id), all_data, by = c('gene_name', 'cell_line'))

  testable_zero_hr = filter(all_data, !is.na(p_adj_0h))
  testable_two_hr = filter(all_data, !is.na(p_adj_2h))
  testable_six_hr = filter(all_data, !is.na(p_adj_6h))
  testable_later_timepoints = unique(c(testable_two_hr$id, testable_six_hr$id))
  
  zero_hr = combine_dat %>% select(simple_id, id, gene_name, region, MODEL_COUNT_0h, total_counts_0h, predicted_state, FILTER_PASS, replicate) %>% rename(MODEL_COUNT = MODEL_COUNT_0h, TOTAL_COUNT = total_counts_0h) %>% mutate(timepoint = 'zero_hr') %>% filter(id %in% testable_zero_hr$id)
  two_hr = combine_dat %>%  select(simple_id, id, gene_name, region, MODEL_COUNT_2h, total_counts_2h, predicted_state, FILTER_PASS, replicate) %>% rename(MODEL_COUNT = MODEL_COUNT_2h, TOTAL_COUNT = total_counts_2h) %>% mutate(timepoint = 'two_hr') %>% filter(id %in% testable_later_timepoints)
  six_hr = combine_dat %>% select(simple_id, id, gene_name, region, MODEL_COUNT_6h, total_counts_6h, predicted_state, FILTER_PASS, replicate) %>% rename(MODEL_COUNT = MODEL_COUNT_6h, TOTAL_COUNT = total_counts_6h) %>% mutate(timepoint = 'six_hr') %>% filter(id %in% testable_later_timepoints)
  
  
  all_snv_dat = rbind(zero_hr, two_hr, six_hr)
  all_snv_dat = all_snv_dat %>% mutate(delta_AR = (MODEL_COUNT/TOTAL_COUNT) - 0.5)
  
  all_snv_dat$timepoint_numeric = all_snv_dat$timepoint
  all_snv_dat$timepoint_numeric[all_snv_dat$timepoint == 'zero_hr'] = '0h'
  all_snv_dat$timepoint_numeric[all_snv_dat$timepoint == 'two_hr'] = '2h'
  all_snv_dat$timepoint_numeric[all_snv_dat$timepoint == 'six_hr'] = '6h'
  all_snv_dat$cell_line = CELL_LINE 
  return(all_snv_dat)
  
}




get_filtered_count = function(CELL_LINE, PASS_PROP_THRESHOLD, return_gene_list = FALSE){
  all_snv_dat = get_snv_count(CELL_LINE)
  rep1_dat = all_snv_dat %>% filter(replicate == 'rep1') %>% rename(MODEL_COUNT_rep1 = MODEL_COUNT, TOTAL_COUNT_rep1 = TOTAL_COUNT, delta_AR_rep1 = delta_AR)
  rep2_dat = all_snv_dat %>% filter(replicate == 'rep2') %>% rename(MODEL_COUNT_rep2 = MODEL_COUNT, TOTAL_COUNT_rep2 = TOTAL_COUNT, delta_AR_rep2 = delta_AR)
  all_rep_dat = left_join(rep1_dat, rep2_dat, by = c('id', 'gene_name', 'predicted_state', 'FILTER_PASS', 'timepoint', 'timepoint_numeric'))
  
  
  all_rep_dat = all_rep_dat %>% mutate(consistency_check = (sign(delta_AR_rep1) == sign(delta_AR_rep2)))
  all_rep_dat$consistency_check[all_rep_dat$delta_AR_rep2 == 0] = TRUE
  all_rep_dat$consistency_check[all_rep_dat$delta_AR_rep1 == 0] = TRUE
  

  total_snvs = all_rep_dat %>% group_by(predicted_state, timepoint, FILTER_PASS)  %>% count(gene_name) %>% rename(snvs_per_gene = n)
  consistent_snvs = all_rep_dat %>% group_by(predicted_state, timepoint, gene_name, FILTER_PASS) %>% filter(consistency_check) %>% count(consistency_check) %>% rename(consistent_snvs_per_gene = n)
  
  
  plot_dat = left_join(total_snvs, consistent_snvs, by = c('predicted_state', 'timepoint', 'gene_name', 'FILTER_PASS'))
  plot_dat$consistent_snvs_per_gene[is.na(plot_dat$consistent_snvs_per_gene)] = 0
  plot_dat$prop_cons = plot_dat$consistent_snvs_per_gene/plot_dat$snvs_per_gene
  plot_dat$ASE_state = "nonASE"
  plot_dat$ASE_state[plot_dat$predicted_state == 'state3'] = "ASE"
  plot_dat$ASE_state[plot_dat$predicted_state == 'state1' & plot_dat$timepoint %in% c('two_hr', 'six_hr')] = "ASE"
  plot_dat$ASE_state[plot_dat$predicted_state == 'state2' & plot_dat$timepoint %in% c('six_hr')] = "ASE"
  plot_dat$ASE_state[plot_dat$predicted_state == 'state4' & plot_dat$timepoint %in% c('zero_hr', 'six_hr')] = "ASE"
  plot_dat$ASE_state[plot_dat$predicted_state == 'state5' & plot_dat$timepoint %in% c('zero_hr')] = "ASE"
  plot_dat$ASE_state[plot_dat$predicted_state == 'state6' & plot_dat$timepoint %in% c('zero_hr', 'two_hr')] = "ASE"
  
  
  plot_dat$timepoint[plot_dat$timepoint == 'zero_hr'] = '0h'
  plot_dat$timepoint[plot_dat$timepoint == 'two_hr'] = '2h'
  plot_dat$timepoint[plot_dat$timepoint == 'six_hr'] = '6h'
  
  plot_dat$FILTER_PASS = paste0('FILTER_PASS = ', plot_dat$FILTER_PASS)
  all_genes = plot_dat %>% group_by(gene_name) %>% count() %>% filter(n==3)
  pass_genes = plot_dat %>% filter(prop_cons >= PASS_PROP_THRESHOLD | ASE_state == 'nonASE')  %>% group_by(gene_name) %>% count() %>% filter(n==3)
  
  pass_genes_count = plot_dat %>% filter(gene_name %in% pass_genes$gene_name) %>% ungroup() %>% select(gene_name, FILTER_PASS, predicted_state) %>% group_by(FILTER_PASS) %>% unique() %>% count(predicted_state)
  total_genes_count = plot_dat  %>% ungroup() %>% select(gene_name, FILTER_PASS, predicted_state) %>% group_by(FILTER_PASS) %>% unique() %>% count(predicted_state) %>% rename(total_genes = n)
  pass_genes_count$PASS_PROP_THRESHOLD = PASS_PROP_THRESHOLD 
  pass_genes_count$cell_line = CELL_LINE
  pass_genes_count = left_join(pass_genes_count, total_genes_count, by = c('predicted_state', 'FILTER_PASS'))
  if(return_gene_list){
    return(pass_genes)
  } else {
  return(pass_genes_count)
  }
}


## output file of snvs in passing genes 
output_res = function(CELL_LINE, PASS_PROP_THRESHOLD){
  all_pass = get_filtered_count(CELL_LINE, PASS_PROP_THRESHOLD, return_gene_list = TRUE)
  all_snv_dat = get_snv_count(CELL_LINE)
  out_res = all_snv_dat %>% filter(gene_name %in% all_pass$gene_name) %>% filter(FILTER_PASS) %>% mutate(cell_line = CELL_LINE) %>% tidyr::separate(id, into = c('chrom', 'pos1', 'strand', 'type_snp'), convert = TRUE, sep = ':') %>% mutate(pos0 = pos1 - 1) %>% select(chrom, pos0, pos1, gene_name, predicted_state, strand, type_snp, cell_line, region) 
  out_res = unique(out_res)
  out_fn = paste0(RESULTS_DIR, '/', CELL_LINE, '_PASS_PROP_THRESHOLD_', PASS_PROP_THRESHOLD, '_gene_categorizations.txt')
  write.table(out_res, file = out_fn, sep = '\t', row.names = F, quote = F)
  return(out_res)
}

out = output_res(CELL_LINE, PASS_PROP_THRESHOLD)
