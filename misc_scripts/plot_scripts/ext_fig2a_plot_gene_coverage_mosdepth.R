library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd('~/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')
cov_dat = fread('../mosdepth_gene_coverage/ASRS_ASRT_MIXED_coverage_comparison_genes_labels.txt')
# color pal
pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'mixed' ='#cbc0d3')

cov_dat = cov_dat %>% select(-ensembl_transcript) %>% unique()
mean_cov_dat = cov_dat %>% group_by(gene_name, gene_state, cell_line, timepoint, replicate) %>% summarize(mean_log2 = mean(log2_count))


## remove polyploid cell lines
mean_cov_dat = mean_cov_dat %>% filter(!cell_line %in% c('PC-9', 'PC-3', 'Caco-2'))

cov_dat %>% group_by(gene_state) %>% summarize(mean_log2 = mean(log2_count))

ggplot(mean_cov_dat, aes(x = timepoint, y = mean_log2, fill = gene_state)) + geom_boxplot() + facet_wrap(.~cell_line) + theme_bw() + scale_fill_manual(values = pal)

p1 = ggplot(mean_cov_dat, aes(x = timepoint, y = mean_log2, fill = gene_state)) + geom_boxplot(aes(x = timepoint, y = mean_log2)) + facet_wrap(.~cell_line) + theme_bw() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1)')
p2 = ggplot(mean_cov_dat %>% filter(gene_state !='mixed'), aes(x = timepoint, y = mean_log2, fill = gene_state)) + geom_boxplot() + facet_wrap(.~cell_line) + theme_bw() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1)')
p3 = ggplot(mean_cov_dat %>% filter(gene_state !='mixed'), aes(x = gene_state, y = mean_log2, fill = gene_state)) + geom_boxplot()  + theme_bw() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1)') + stat_compare_means()

#ggplot(mean_cov_dat %>% filter(gene_state !='mixed'), aes(x = gene_state, y = mean_log2)) + geom_boxplot()  + facet_wrap(.~timepoint) + theme_bw() + scale_fill_manual(values = pal) + ylab('mean(log2(count)') + stat_compare_means()

mean_cov_dat %>% group_by(gene_state) %>% summarize(mean_log2 = mean(mean_log2))


p1 = ggplot(mean_cov_dat %>% filter(timepoint == 'zero_hr'), aes (x = gene_state, y = mean_log2, fill = gene_state)) + geom_boxplot(aes(x = gene_state, y = mean_log2)) + facet_wrap(.~cell_line) + 
  theme_classic() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1)') + ggtitle('0h') #+ stat_compare_means(label = 'p.signif', comparisons = list(c('asRS', 'asRT'), c('asRS', 'mixed'), c('asRT', 'mixed')), hide.ns = TRUE, symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) + ggtitle('0h')

p2 = ggplot(mean_cov_dat %>% filter(timepoint == 'two_hr'), aes (x = gene_state, y = mean_log2, fill = gene_state)) + geom_boxplot(aes(x = gene_state, y = mean_log2)) + facet_wrap(.~cell_line) + 
  theme_classic() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1)') + ggtitle('2h') #+ stat_compare_means(label = 'p.signif', comparisons = list(c('asRS', 'asRT'), c('asRS', 'mixed'), c('asRT', 'mixed')), hide.ns = TRUE, symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) + ggtitle('2h')

p3 = ggplot(mean_cov_dat %>% filter(timepoint == 'six_hr'), aes (x = gene_state, y = mean_log2, fill = gene_state)) + geom_boxplot(aes(x = gene_state, y = mean_log2)) + facet_wrap(.~cell_line) + 
  theme_classic() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1)') + ggtitle('6h') #+ stat_compare_means(label = 'p.signif', comparisons = list(c('asRS', 'asRT'), c('asRS', 'mixed'), c('asRT', 'mixed')), hide.ns = TRUE, symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) + ggtitle('6h')


mean_across_timepoints = mean_cov_dat %>% group_by(gene_name, gene_state, cell_line, replicate) %>% summarize(mean_log2_across_timepoints = mean(mean_log2))

p4 = ggplot(mean_across_timepoints, aes (x = gene_state, y = mean_log2_across_timepoints, fill = gene_state)) + geom_boxplot(aes(x = gene_state, y = mean_log2_across_timepoints)) + facet_wrap(.~cell_line, nrow = 2) + 
  theme_classic() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1) across timepoints') #+ stat_compare_means(label = 'p.signif', comparisons = list(c('asRS', 'asRT'), c('asRS', 'mixed'), c('asRT', 'mixed')), hide.ns = TRUE, symnum.args <- list(cutpoints = c(0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns"))) 


manuscript_fig = ggplot(mean_cov_dat, aes (x = cell_line, y = mean_log2, fill = gene_state)) + geom_boxplot(aes(x = cell_line, y = mean_log2)) +
  theme_classic() + scale_fill_manual(values = pal) + ylab('Mean log2(read count + 1)') + xlab('Cell line') # #+ stat_compare_means(label = 'p.signif', comparisons = list(c('asRS', 'asRT'), c('asRS', 'mixed'), c('asRT', 'mixed')), hide.ns = TRUE, symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) + ggtitle('6h')




 write.table(mean_cov_dat, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig2a.txt', sep = '\t', row.names = F, quote = F)


pdf('../MANUSCRIPT_FIGURES/extended_fig2/extended_fig2a.pdf', width = 9, height = 4)
manuscript_fig + theme(legend.position = 'bottom')
dev.off()




library(cowplot)
plot_grid(p1,p2,p3, nrow = 1)
pdf('../figures/asrs_asrt_mosdepth_gene_cov.pdf', height = 8.5, width = 10)
p1
p2
p3
dev.off()


pdf('../figures/asrs_asrt_mosdepth_gene_cov_plot_grid.pdf', height = 10, width = 25)
plot_grid(p1 + theme(legend.position = 'none'),p2 + theme(legend.position = 'none'),p3 + theme(legend.position = 'none'), nrow = 1)
dev.off()

pdf('../figures/asrs_asrt_mosdepth_gene_cov_plot_grid_average_across_timepoints.pdf', height = 4, width = 10)
p4 + theme(legend.position ='none')
dev.off()



