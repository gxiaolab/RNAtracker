rm(list=ls())
library(data.table)
library(dplyr)
library(gtools)
source("../../../../softwares/loadRData.R")
library(doParallel)
library(foreach)
library(ggplot2)
library(cowplot)


COVERAGE_CUTOFF = 10
MISMATCH_CUTOFF = 2
SNP_COUNT_FILTER = 2


MAX_PROB_CUTOFF = 0.95
TOP2DIFF_CUTOFF = 1
STATE_QUANT_CUTOFF = 1

options(scipen = 99)

make_plot = function(MU0, RHO0, MU1, RHO1){
  OUT_LABEL = paste0("MU0_", MU0, "_RHO0_", RHO0, "_MU1_", MU1, "_RHO1")
  
  in_fn = paste0('../results/MBASED_PURE_SIM_RESULTS_', OUT_LABEL, '.txt')
  combined_res = fread(in_fn)
  combined_res$gene_state = combined_res$state
  
  combined_res = combined_res %>% select(-quant_cutoff) %>% unique()
  
  combined_res$label = OUT_LABEL
  
  pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'mixed' = '#cbc0d3', 'nonASE' = 'grey')
  
  p1 = ggplot(combined_res, aes(x = state, y = Precision, fill = gene_state)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = pal)
  
  p2 = ggplot(combined_res, aes(x = state, y = Recall, fill = gene_state)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = pal)
  
  p3 = ggplot(combined_res, aes(x = state, y = Recall_classified, fill = gene_state)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = pal)
  print(combined_res %>% group_by(gene_state) %>% summarize(mean_Precision = mean(Precision), mean_Recall = mean(Recall)))
  return(combined_res)
}
res1 = make_plot(0.5, 0.004, 0.7, 0.004)
res2 = make_plot(0.5, 0.004, 0.8, 0.004)
res3 = make_plot(0.5, 0.004, 0.9, 0.004)

all_res = rbind(res1, res2, res3)

pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'mixed' = '#cbc0d3', 'nonASE' = 'grey')
p1 = ggplot(all_res, aes(x = label, y = Precision, fill = gene_state)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = pal) + xlab('Simulation hyperparameters') + ylim(0.5, 1)
p2 = ggplot(all_res, aes(x = label, y = Recall, fill = gene_state)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = pal) + xlab('Simulation hyperparameters') + ylim(0.5, 1)
p3 = ggplot(all_res, aes(x = label, y = Recall_classified, fill = gene_state)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = pal) + xlab('Simulation hyperparameters') + ylim(0.5, 1)


out_fn = paste0('../figures/pure_sim_res_precision_recall_boxplot_MBASED_HYPERPARAM_SETUPS_STRINGENT_CUTOFFS.pdf')


pdf(out_fn, width = 10)
plot_grid(p1,p2,p3, nrow = 3)
dev.off()


pdf('../MANUSCRIPT_FIGURES/extended_fig3/extended_fig3a.pdf', height = 5, width = 5)
p1 + theme(legend.position = 'none')
dev.off()

pdf('../MANUSCRIPT_FIGURES/extended_fig3/extended_fig3b.pdf', height = 5, width = 5)
p2 + theme(legend.position = 'none')
dev.off()


pdf('../MANUSCRIPT_FIGURES/extended_fig3/extended_fig3c.pdf', height = 5, width = 5)
p3 + ylim(0.9, 1) + theme(legend.position = 'none')
dev.off()

 write.table(all_res, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig3.txt', sep = '\t', row.names = F, quote = F)


