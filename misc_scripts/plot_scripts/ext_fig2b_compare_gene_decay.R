rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)


load('../decay_results/ALL_CELL_LINES_asrs_asrt_mixed_default_model_half_life_estimate.RData')
all_snv_fn = '../bed_files/all_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed'
all_snvs = fread(all_snv_fn)
names(all_snvs) = c('chr', 'pos0', 'pos1', 'gene_name', 'gene_state', 'strand', 'alleles', 'cell_line', 'region')

gene_status = all_snvs
gene_status$gene_state[gene_status$gene_state == 'state3'] = 'asRT'
gene_status$gene_state[gene_status$gene_state %in% c('state1', 'state2')] = 'asRS'
gene_status$gene_state[gene_status$gene_state %in% c('state4', 'state5', 'state6')] = 'mixed'
gene_status = filter(gene_status, gene_state %in% c('asRT', 'asRS', 'mixed'))
#gene_status = filter(gene_status, gene_state %in% c('asRT', 'asRS'))
pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'mixed' ='#cbc0d3')

plot_dat = left_join(gene_status %>% select(gene_name, gene_state, cell_line) %>% unique(), all_res_df, by = c('cell_line', 'gene_name'))

## remove polyploid cell lines
plot_dat = plot_dat %>% filter(!cell_line %in% c('PC-9', 'PC-3', 'Caco-2'))


plot_dat$half_life = as.numeric(plot_dat$half_life)
p1 = ggplot(plot_dat, aes(x = gene_state, y = half_life, fill = gene_state)) + geom_boxplot() + theme_bw()  + 
  facet_wrap(.~cell_line, nrow = 2) + scale_fill_manual(values = pal) + theme_classic()
p2 = ggplot(plot_dat, aes(x =  gene_state, y = half_life, fill = gene_state)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = pal)



manuscript_fig = ggplot(plot_dat, aes(x =  cell_line, y = half_life, fill = gene_state)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = pal) + ylab('Estimated half-life') + xlab('Cell_line') # #+ stat_compare_means(label = 'p.signif', comparisons = list(c('asRS', 'asRT'), c('asRS', 'mixed'), c('asRT', 'mixed')), hide.ns = TRUE, symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) + ggtitle('6h')



 write.table(plot_dat, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig2b.txt', sep = '\t', row.names = F, quote = F)

pdf('../MANUSCRIPT_FIGURES/extended_fig2/extended_fig2b.pdf', width = 9, height = 4)
manuscript_fig 
dev.off()



pdf('../figures/asrs_asrt_gene_decay.pdf', height = 4, width = 10)
p1 + theme(legend.position = 'none') + xlab('Gene state') + ylab('Half-life')
p2 + theme(legend.position = 'none')
dev.off()

