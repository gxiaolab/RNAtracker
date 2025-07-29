plot_gene_state_prev_intronic_ASE_comparison.Rrm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(gtools)
library(cowplot)

#setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')


load('../results/bru_gene_state_prev.RData') # bru_plot_dat, bru_stars_dat 

intronic_dat = fread('../results/include_intron_new_ASE_gene_count_0h.txt')

## need to apply 0.5 PROP threshold
all_intronic_res = fread('../bed_files/all_snvs_intronic_0h_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed')


# remove polyploid cell lines
intronic_dat = filter(intronic_dat, !CELL_LINE %in% c('PC-3', 'PC-9', 'Caco-2'))
bru_plot_dat = filter(bru_plot_dat, !cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))
bru_stars_dat = filter(bru_stars_dat, !cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))



bru_asrt = filter(bru_plot_dat, `gene type` == 'asRT')
intronic_dat$total_intron_classified = intronic_dat$intron_new_ASE + intronic_dat$intron_new_nonASE
intronic_dat = intronic_dat %>% rename(cell_line = CELL_LINE) %>% mutate(`gene type` = 'intronic asRT') %>% rename(total_state_genes = intron_new_ASE) %>% rename(total_classified_genes = total_intron_classified) %>% mutate(prevalence = total_state_genes/total_classified_genes)





# plot "adjusted" asRS prevalence
# we need to use the mixed row counts twice for both adjusted asRS and adjusted asRT
tmp_dat = rbind(bru_plot_dat, intronic_dat %>% select(names(bru_plot_dat)))
asrs_dat = filter(tmp_dat, `gene type` %in% c('asRS', 'mixed'))
#asrs_dat_summarized = asrs_dat %>% group_by(cell_line) %>% summarize(adjusted_prev = sum(prevalence), total_state_genes = sum(total_state_genes), total_classified_genes = total_classified_genes) %>% unique() 
asrs_dat_summarized = asrs_dat %>% group_by(cell_line) %>% summarize(total_state_genes = sum(total_state_genes), total_classified_genes = total_classified_genes) %>% mutate(adjusted_prev = total_state_genes/total_classified_genes) %>% unique()

# plot "adjusted" asRT prevalence WITHOUT intronic asRT 
orig_asrt_dat = filter(tmp_dat, `gene type` %in% c('asRT', 'mixed'))
orig_asrt_dat_summarized = orig_asrt_dat %>% group_by(cell_line) %>% summarize(total_state_genes = sum(total_state_genes), total_classified_genes = total_classified_genes) %>% mutate(adjusted_prev = total_state_genes/total_classified_genes) %>% unique()

# plot "adjusted" asRT prevalence WITH intronic asRT 
asrt_dat = filter(tmp_dat, `gene type` %in% c('asRT', 'mixed', 'intronic asRT'))
#asrt_dat_summarized = asrt_dat %>% group_by(cell_line) %>% summarize(adjusted_prev = sum(total_state_genes)/sum(unique(total_classified_genes)), total_state_genes = sum(total_state_genes), total_classified_genes = sum(unique(total_classified_genes))) 
asrt_dat_summarized = asrt_dat %>% group_by(cell_line) %>% summarize(total_state_genes = sum(total_state_genes), total_classified_genes = sum(unique(total_classified_genes))) %>% mutate(adjusted_prev = total_state_genes/total_classified_genes) %>% unique()





## actually plot mixed prev separately like in the original fig
# we need to use the mixed row counts twice for both adjusted asRS and adjusted asRT
#tmp_dat = rbind(bru_plot_dat, intronic_dat %>% select(names(bru_plot_dat)))
#asrs_dat = filter(tmp_dat, `gene type` %in% c('asRS'))
#asrs_dat_summarized = asrs_dat %>% group_by(cell_line) %>% summarize(adjusted_prev = sum(prevalence), total_state_genes = sum(total_state_genes), total_classified_genes = total_classified_genes) %>% unique() 
#asrt_dat = filter(tmp_dat, `gene type` %in% c('asRT', 'intronic asRT'))
#asrt_dat_summarized = asrt_dat %>% group_by(cell_line) %>% summarize(adjusted_prev = sum(total_state_genes)/sum(unique(total_classified_genes)), total_state_genes = sum(total_state_genes), total_classified_genes = sum(unique(total_classified_genes))) 




combine_dat = rbind(asrt_dat_summarized %>% mutate(gene_type = 'adjusted asRT'), orig_asrt_dat_summarized %>% mutate(gene_type = 'adjusted asRT (exclude intron-based asRT)'), asrs_dat_summarized %>% mutate(gene_type = 'adjusted asRS'))




#### PLOTTING
theme_set(theme_classic())
#pal = c('asRS' = '#875f9a' , 'adjusted asRT' = '#48d1cc', 'mixed' = '#cbc0d3')
#pal = c('adjusted asRS' = '#875f9a' , 'adjusted asRT' = '#48d1cc')

pal = c('asRS + mixed' = '#875f9a' , 'asRT + intron-based asRT + mixed' = '#48d1cc', 'asRT + mixed' = '#90e0ef')
combine_dat$gene_type[combine_dat$gene_type == 'adjusted asRS'] = 'asRS + mixed'
combine_dat$gene_type[combine_dat$gene_type == 'adjusted asRT'] = 'asRT + intron-based asRT + mixed'
combine_dat$gene_type[combine_dat$gene_type == 'adjusted asRT (exclude intron-based asRT)'] = 'asRT + mixed'



 write.table(combine_dat, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig4c.txt', sep = '\t', row.names = F, quote = F)

cell_line_order = combine_dat %>% filter(gene_type == 'asRS + mixed') %>% arrange(desc(adjusted_prev))
combine_dat$cell_line = factor(combine_dat$cell_line, levels = cell_line_order$cell_line)

p3 = ggplot(combine_dat) +   geom_linerange(aes(x = cell_line, ymin = 0, ymax = adjusted_prev, group = gene_type), position = position_dodge(width = 0.75))+
  geom_point(pch = 21, color = 'black', size = 6, aes(x = cell_line, y = adjusted_prev, fill = gene_type, group = gene_type), position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pal, breaks = c('asRS + mixed', 'asRT + intron-based asRT + mixed', 'asRT + mixed'), name = 'Gene Type')  + ylab('Prevalence')
combine_dat %>% group_by(cell_line) %>% summarize(combine_prev = sum(adjusted_prev)) %>% arrange(combine_prev)


pdf('../figures/asrs_asrt_lollipop_plot_prev_comparison_with_and_without_intronic_ASE_0h.pdf', width = 8, height = 5)
p3 + theme(legend.position = 'bottom') + xlab('Cell line')
dev.off()


