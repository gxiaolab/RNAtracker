rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)

COVERAGE_CUTOFF = 10
MISMATCH_CUTOFF = 2
ASE_0h_METHOD = 'original'
ASE_0h_METHOD = 'original'
output_label = 'GATK_VAR_merge'
output_label = paste0(output_label, "_", COVERAGE_CUTOFF, "_", MISMATCH_CUTOFF)
filter_label = 'MAX_PROB_CUTOFF_0.95_TOP2DIFF_CUTOFF_1_STATE_QUANT_CUTOFF_0.333333333333333'


CELL_LINE = 'GM12878'
WHICH_MODEL = 'ASE_0h'

results_dir = paste0('../REAL_DATA_balanced_g1/', COVERAGE_CUTOFF, '_', MISMATCH_CUTOFF, 'REQUIRE_TESTABLE_2_TIMEPOINTS/ASE_0h_method_', ASE_0h_METHOD, '/UNIFORM_CUTOFFS_', filter_label,'/')



get_PI_est = function(CELL_LINE, WHICH_MODEL){
output_fn = paste0(CELL_LINE, '_', output_label,  '_', WHICH_MODEL, '_model/')
in_fn = paste0(results_dir, CELL_LINE, '_', output_label, '_', WHICH_MODEL, '_model.est.Rdata')
load(in_fn)
out_res = data.frame(cell_line = CELL_LINE, PI_EST = final_step_param$pi_est, number_genes_tested = length(final_step_post$gene_id))
out_res$state = rownames(out_res)
return(out_res)
}

get_number_tested_genes_0h = function(CELL_LINE){
  zero_hr_results = paste0('../REAL_DATA_balanced_g1/', COVERAGE_CUTOFF, '_', MISMATCH_CUTOFF, '/zero_hr_only/', CELL_LINE, '_GATK_VAR_merge_10_2_zero_hr_only.est.Rdata')
  load(zero_hr_results)
  out_res = data.frame(cell_line = CELL_LINE, number_genes_tested = length(final_step_post$gene_id))
}




CELL_LINES = fread('../cell_lines_14_final', header = F)

out_ASE = do.call("rbind", lapply(CELL_LINES$V1, get_PI_est,'ASE_0h'))
out_nonASE = do.call("rbind", lapply(CELL_LINES$V1, get_PI_est,'nonASE_0h'))
all_res = rbind(out_ASE, out_nonASE)
total_genes_tested = all_res %>% dplyr::select(cell_line, number_genes_tested) %>% unique() %>% group_by(cell_line) %>% summarize(total_genes_tested = sum(number_genes_tested))

all_res_0h = do.call("rbind", lapply(CELL_LINES$V1, get_number_tested_genes_0h))
all_res_0h %>% arrange(number_genes_tested) # less clean of a cut-off to remove Caco-2 and PC-9 because these are over 100
# add K562 and Panc1 as 0


total_genes_tested = rbind(total_genes_tested, data.frame(cell_line = c('K562', 'Panc1'), total_genes_tested = c(0,0)))
total_genes_tested = total_genes_tested %>% arrange(total_genes_tested)
total_genes_tested$cell_line = factor(total_genes_tested$cell_line, levels = total_genes_tested$cell_line)

p0 = ggplot(total_genes_tested, aes(x = cell_line, y = total_genes_tested, label = total_genes_tested)) + 
  geom_bar(stat = 'identity', color = 'black', fill = 'steelblue', alpha = 0.67) + geom_label() + ylab('Total tested genes') + xlab('Cell line') + theme_classic()

pdf('../figures/total_genes_tested_stage2.pdf', width = 10, height = 5)
p0 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('../MANUSCRIPT_FIGURES/extended_fig1/extended_fig1b.pdf', width = 5, height = 6.5)
p0 + coord_flip() # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


## also plot # of genes that are ultimately classified after all the filtering
all_snvs = fread('../bed_files/all_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.25.bed')
total_genes_classified = all_snvs %>% select(V4, V5, V8) %>% unique() %>% group_by(V8) %>% count()
names(total_genes_classified) = c('cell_line', 'total_genes_classified')
p1 = ggplot(total_genes_classified, aes(x = cell_line, y = total_genes_classified, label = total_genes_classified)) + geom_bar(stat = 'identity') + geom_label()





all_res = left_join(all_res, total_genes_tested, by = 'cell_line') %>% mutate(scaling_factor = number_genes_tested/total_genes_tested)

all_res$adjusted_PI_est = all_res$PI_EST * all_res$scaling_factor

all_res$gene_state = all_res$state
all_res$gene_state[all_res$gene_state == 'comp0.post'] = 'nonASE'
all_res$gene_state[all_res$gene_state %in% c('comp1.post', 'comp2.post')] = 'asRS'
all_res$gene_state[all_res$gene_state %in% c('comp3.post')] = 'asRT'
all_res$gene_state[all_res$gene_state %in% c('comp4.post', 'comp5.post', 'comp6.post')] = 'mixed'
sum_all_res = all_res %>% group_by(cell_line, gene_state) %>% summarize(adjusted_PI_est = sum(adjusted_PI_est))
pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'mixed' = '#cbc0d3')
p = ggplot(sum_all_res %>% filter(gene_state!='nonASE') %>% filter(!cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))) + geom_linerange(aes(x = cell_line, ymin = 0, ymax = adjusted_PI_est, group = gene_state),  position = position_dodge(width = 0.75))+
  geom_point( pch = 21, color = 'black', size = 6, aes(x = cell_line, y = adjusted_PI_est, fill = gene_state, group = gene_state), position = position_dodge(width = 0.75)) + scale_fill_manual(values = pal)

## plot original prev plot just for comparison
load('../results/bru_gene_state_prev.RData')
p2 = ggplot(bru_plot_dat %>% filter(`gene type`!='nonASE') %>% filter(!cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))) + geom_linerange(aes(x = cell_line, ymin = 0, ymax = prevalence, group = `gene type`),  position = position_dodge(width = 0.75))+
  geom_point( pch = 21, color = 'black', size = 6, aes(x = cell_line, y = prevalence, fill = `gene type`, group = `gene type`), position = position_dodge(width = 0.75)) + scale_fill_manual(values = pal)



# recalculate prev using total genes tested
combine_dat = left_join(bru_plot_dat, all_res %>% select(cell_line, total_genes_tested) %>% unique(), by = 'cell_line' ) %>% mutate(adj_prev = total_state_genes/total_genes_tested)
p3 = ggplot(combine_dat %>% filter(`gene type`!='nonASE') %>% filter(!cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))) + geom_linerange(aes(x = cell_line, ymin = 0, ymax = adj_prev, group = `gene type`),  position = position_dodge(width = 0.75))+
  geom_point( pch = 21, color = 'black', size = 6, aes(x = cell_line, y = adj_prev, fill = `gene type`, group = `gene type`), position = position_dodge(width = 0.75)) + scale_fill_manual(values = pal)




pdf('../figures/prev_from_pi_est_lollipop_plot.pdf', width = 10, height = 10)
plot_grid(p2 + ylim(c(0,0.5)),p3 + ylim(c(0, 0.5)), p + ylim(c(0, 0.5)), nrow=3)
dev.off()



