rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)

setwd('~/Desktop/project-gxxiao4.2/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')

ALL_COUNT_DAT = fread('../results/full_ALL_GENES_counts_data.txt')

MAX_PROB_CUTOFF = 0.95
STATE_QUANT_CUTOFF = 0.333333333333333
TOP2DIFF_CUTOFF = 1
SNP_COUNT_FILTER = 2


process_res = function(CELL_LINE, process_dat_fn){

load(process_dat_fn)
all_df = final_step_post
total_genes = nrow(final_step_post)
all_df$predicted_state = paste0('state', all_df$predicted_state)
all_df$cell_line = CELL_LINE
all_df$SNP_COUNT_FILTER = SNP_COUNT_FILTER
# filter by requiring max probability to be at least some amount greater than the second max probability
all_df$second_max_prob = apply(all_df %>% select(contains('comp')), 1, function(x) max( x[x!=max(x)] ))
all_df$top_prob_diff = all_df$max_probability - all_df$second_max_prob
all_df$second_predicted_state = apply(all_df %>% select(contains('comp')), 1, function (x) which(x == max(x[x!=max(x)]))[1]) # some rows have ties, so just take the first item
# subtract 1 for indexing purposes
all_df$second_predicted_state = paste0('state', as.numeric(all_df$second_predicted_state) - 1)


# scale each predicted state's max probability
predicted_states = unique(all_df$predicted_state)
tmp_list = list()
for (spec_state in predicted_states){
tmp_df = all_df %>% filter(predicted_state == spec_state)
tmp_df$scaled_prob = scale(tmp_df$max_probability)
tmp_list[[spec_state]] = tmp_df
}
all_df_scaled = do.call("rbind", tmp_list)

cutoff_df = all_df_scaled %>% group_by(predicted_state) %>% summarise(quantile = scales::percent(c(STATE_QUANT_CUTOFF )),quant_cutoff = quantile(max_probability, c(STATE_QUANT_CUTOFF )))
all_df_scaled = left_join(all_df_scaled, cutoff_df, by = 'predicted_state')

zero_hr_assignments = all_df_scaled %>% filter(max_probability >= MAX_PROB_CUTOFF | top_prob_diff >= TOP2DIFF_CUTOFF | max_probability >= quant_cutoff)
return(zero_hr_assignments)

}

output_results = function(CELL_LINE, requested_output){
  old_result_fn = paste0("../REAL_DATA_balanced_g1/10_2/zero_hr_only/", CELL_LINE, "_GATK_VAR_merge_10_2_zero_hr_only.est.Rdata")
  
  new_result_fn = paste0("../REAL_DATA_balanced_g1_ALT_SPLICED_REMOVED/10_2/zero_hr_only/", CELL_LINE, "_GATK_VAR_merge_10_2_zero_hr_only.est.Rdata")
  
  old_ASE_result_fn = paste0("../REAL_DATA_balanced_g1/10_2REQUIRE_TESTABLE_2_TIMEPOINTS/ASE_0h_method_original/UNIFORM_CUTOFFS_MAX_PROB_CUTOFF_0.95_TOP2DIFF_CUTOFF_1_STATE_QUANT_CUTOFF_0.333333333333333/", CELL_LINE, "_GATK_VAR_merge_10_2_ASE_0h_model.est.Rdata")
  
  new_ASE_result_fn = paste0("../REAL_DATA_balanced_g1_ALT_SPLICED_REMOVED/10_2REQUIRE_TESTABLE_2_TIMEPOINTS/ASE_0h_method_original/UNIFORM_CUTOFFS_MAX_PROB_CUTOFF_0.95_TOP2DIFF_CUTOFF_1_STATE_QUANT_CUTOFF_0.333333333333333/", CELL_LINE, "_GATK_VAR_merge_10_2_ASE_0h_model.est.Rdata")
  
  old_nonASE_result_fn = paste0("../REAL_DATA_balanced_g1/10_2REQUIRE_TESTABLE_2_TIMEPOINTS/ASE_0h_method_original/UNIFORM_CUTOFFS_MAX_PROB_CUTOFF_0.95_TOP2DIFF_CUTOFF_1_STATE_QUANT_CUTOFF_0.333333333333333/", CELL_LINE, "_GATK_VAR_merge_10_2_nonASE_0h_model.est.Rdata")
  
  new_nonASE_result_fn = paste0("../REAL_DATA_balanced_g1_ALT_SPLICED_REMOVED/10_2REQUIRE_TESTABLE_2_TIMEPOINTS/ASE_0h_method_original/UNIFORM_CUTOFFS_MAX_PROB_CUTOFF_0.95_TOP2DIFF_CUTOFF_1_STATE_QUANT_CUTOFF_0.333333333333333/", CELL_LINE, "_GATK_VAR_merge_10_2_nonASE_0h_model.est.Rdata")
  
  
process_res_old = process_res(CELL_LINE, old_result_fn)
process_res_new = process_res(CELL_LINE, new_result_fn)


comparison_dat = inner_join(process_res_old %>% select(gene_id, predicted_state) %>% rename(old_state = predicted_state), process_res_new %>% select(gene_id, predicted_state) %>% rename(new_state = predicted_state), by = 'gene_id')

table(comparison_dat$old_state == comparison_dat$new_state) 


## now eval stage 2
process_res_old_ASE = process_res(CELL_LINE, old_ASE_result_fn)
process_res_new_ASE = process_res(CELL_LINE, new_ASE_result_fn)


process_res_old_nonASE = process_res(CELL_LINE, old_nonASE_result_fn)
process_res_new_nonASE = process_res(CELL_LINE, new_nonASE_result_fn)


all_old = rbind(process_res_old_ASE %>% select(gene_id, predicted_state, max_probability), process_res_old_nonASE %>% select(gene_id, predicted_state, max_probability))
all_new = rbind(process_res_new_ASE %>% select(gene_id, predicted_state, max_probability), process_res_new_nonASE %>% select(gene_id, predicted_state, max_probability))
comparison_dat = full_join(all_old %>% select(gene_id, predicted_state) %>% rename(old_state = predicted_state), all_new %>% select(gene_id, predicted_state) %>% rename(new_state = predicted_state), by = 'gene_id')
table(comparison_dat$old_state == comparison_dat$new_state) 

testable_in_old_only = filter(comparison_dat, is.na(new_state))
testable_in_new_only = filter(comparison_dat, is.na(old_state))
label_is_same = filter(comparison_dat, old_state == new_state)
label_changes = filter(comparison_dat, old_state != new_state)
summary_dat = data.frame(cell_line = CELL_LINE, testable_in_old_only = nrow(testable_in_old_only), 
                         testable_in_new_only = nrow(testable_in_new_only), label_is_same = nrow(label_is_same), label_changes = nrow(label_changes))

## assess genes with label changes
alt_info_fn = paste0('../PSI_calculator/', CELL_LINE, '_alt_spliced.testable_snvs')
alt_info = fread(alt_info_fn)
alt_info = mutate(alt_info, id = paste0(V1, ':', V3, ':', V6))

m_dat = filter(ALL_COUNT_DAT, cell_line == CELL_LINE, gene_name %in% label_changes$gene_id) %>% select(gene_name, id, contains('AR')) %>% melt(id.vars = c('gene_name', 'id'))
m_dat = m_dat %>% tidyr::separate(variable, into = c('rep', 'type', 'timepoint'), sep = '_')
m_dat = left_join(m_dat, label_changes %>% mutate(label = paste0(old_state, '|', new_state)) %>% rename(gene_name = gene_id), by = 'gene_name')
m_dat = mutate(m_dat, simple_id = gsub(':[A-Z]>[A-Z]', '', id))
m_dat$alt_spliced = m_dat$simple_id %in% alt_info$id

m_dat = m_dat %>% mutate(label = paste0(gene_name, "\n", label))
p = ggplot(m_dat, aes(x = timepoint, y = value, group = id, color = alt_spliced)) + geom_point() + geom_line() + facet_wrap(.~label) + theme_bw() + ylab('allelic ratio') + geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'red') + ggtitle(CELL_LINE)
#out = list(summary_dat, p)
if(requested_output == "summary_dat"){
  out = summary_dat
} else{
  out = p
}
return(out)
}


cell_lines = fread('../cell_lines_11_final', header = FALSE)
all_res = lapply(cell_lines$V1, output_results, "summary_dat")
all_res_df = do.call("rbind", all_res)
all_res_df$prop_label_change = all_res_df$label_changes/(all_res_df$label_changes + all_res_df$label_is_same)

all_figs = lapply(cell_lines$V1, output_results, "fig")

pdf('../figures/genes_with_state_changes_after_removing_alt_spliced.pdf', height = 10, width = 10)
all_figs
dev.off()


## make barplot for counts of genes affected by alt spliced regions
m_dat = all_res_df %>% select(cell_line, label_is_same, label_changes) %>% melt(id.vars = 'cell_line')

names(m_dat) = c('cell_line', 'gene_state_change', 'number_genes')
m_dat$gene_state_change = as.character(m_dat$gene_state_change)
m_dat$gene_state_change[m_dat$gene_state_change == 'label_is_same'] = 'NO'
m_dat$gene_state_change[m_dat$gene_state_change == 'label_changes'] = 'YES'
m_dat = m_dat %>% arrange(cell_line)
m_dat$cell_line = factor(m_dat$cell_line, levels = rev(unique(m_dat$cell_line)))

cell_line_order = m_dat %>% arrange(desc(number_genes)) %>% filter(gene_state_change == 'NO') %>% select(cell_line)
m_dat$cell_line = factor(m_dat$cell_line, levels = rev(cell_line_order$cell_line))
p = ggplot(m_dat, aes(cell_line, y = number_genes, fill = gene_state_change)) + geom_bar(position = 'stack', stat = 'identity', color = 'black') + coord_flip() + theme_classic() + scale_fill_brewer(palette = 'Blues') + theme(legend.position = 'bottom')



 write.table(m_dat, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig2c.txt', sep = '\t', row.names = F, quote = F)


pdf('../figures/barplot_genes_before_after_remove_alt_spliced.pdf', width = 9, height = 6.5)
p + ylab('Number of Genes') + xlab('Cell line')
dev.off()


pdf('../MANUSCRIPT_FIGURES/extended_fig2/extended_fig2c.pdf', width = 11, height = 4.5)
p + ylab('Number of Genes') + xlab('Cell line')
dev.off()



