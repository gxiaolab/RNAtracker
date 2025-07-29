rm(list=ls())
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

#combine_dat = left_join(bru_asrt %>% rename(CELL_LINE = cell_line), intronic_dat, by = 'CELL_LINE') %>% as.data.frame() %>% mutate(adjusted_prop = (total_state_genes + intron_new_ASE) / (total_classified_genes + total_intron_classified)) 
#combine_dat$`gene type` = 'adjusted asRT'
#tmp1 = combine_dat %>% select(CELL_LINE, intron_new_ASE, total_state_genes) 



## compare new ASE genes (from 0h intronic) to previously identified ASE genes
#m_dat = tmp1 %>% melt(id.vars = 'CELL_LINE')
#m_dat$variable = as.character(m_dat$variable)
#m_dat$variable[m_dat$variable == 'intron_new_ASE'] = 'intron ASE 0h'
#m_dat$variable[m_dat$variable == 'total_state_genes'] = 'original asRT'

plot_dat1 = rbind(bru_asrt, intronic_dat %>% select(names(bru_asrt)))
p1 = ggplot(plot_dat1, aes(x = cell_line, y = total_state_genes)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black', aes(fill = `gene type`)) + 
  ylab('number of genes') + scale_fill_manual(values = c('#48d1cc', '#48d1cc40'))

p2 = ggplot(plot_dat1, aes(x = cell_line, y = prevalence)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black', aes(fill = `gene type`)) + 
  ylab('Prevalence') + scale_fill_manual(values = c('#48d1cc', '#48d1cc40'))


## or do one line plot and one bar plot
plot_dat2 = plot_dat1 %>% filter(`gene type` == 'intronic asRT')
plot_dat2 = plot_dat2 %>% arrange(desc(total_state_genes))
plot_dat2$cell_line = factor(plot_dat2$cell_line, levels = plot_dat2$cell_line)

part1 = ggplot(plot_dat2, aes(cell_line, y = total_state_genes )) + 
  geom_bar(stat = 'identity', fill = '#90e0ef', color = 'black', alpha = 0.67)  + 
  theme_classic() + ylab('Number of intronic asRT genes') + xlab('Cell line') 

part2 = ggplot(plot_dat2, aes(cell_line, y = prevalence, group = 1)) + 
  geom_line(stat = 'identity', color = 'red')  + 
  geom_point(fill = 'black', color = 'black', size = 2, pch = 21) + 
  theme_classic() + ylab('Proportion') + xlab('Cell line') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p4 = plot_grid(part2, part1, nrow = 2, rel_heights = c(0.2, 0.8))



write.table(plot_dat2, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig4b.txt', sep = '\t', row.names = F, quote = F)



pdf('../MANUSCRIPT_FIGURES/extended_fig4/extended_fig4b.pdf')
p4
dev.off()



# plot "adjusted" asRT prevalence
# we need to use the mixed row counts twice for both adjusted asRS and adjusted asRT
tmp_dat = rbind(bru_plot_dat, intronic_dat %>% select(names(bru_plot_dat)))
asrs_dat = filter(tmp_dat, `gene type` %in% c('asRS', 'mixed'))
#asrs_dat_summarized = asrs_dat %>% group_by(cell_line) %>% summarize(adjusted_prev = sum(prevalence), total_state_genes = sum(total_state_genes), total_classified_genes = total_classified_genes) %>% unique() 
asrs_dat_summarized = asrs_dat %>% group_by(cell_line) %>% summarize(total_state_genes = sum(total_state_genes), total_classified_genes = total_classified_genes) %>% mutate(adjusted_prev = total_state_genes/total_classified_genes) %>% unique()

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




combine_dat = rbind(asrt_dat_summarized %>% mutate(gene_type = 'adjusted asRT'), asrs_dat_summarized %>% mutate(gene_type = 'adjusted asRS'))

## now compare the prevalence using fisher's exact test
CELL_LINES = unique(combine_dat$cell_line)
OR_vec = rep(0, length(CELL_LINES))
p_vec = rep(0, length(CELL_LINES))
for (i in 1:length(CELL_LINES)){
  spec_cell_line = CELL_LINES[i]
  print(spec_cell_line)
  spec_row_asrs = combine_dat %>% filter(cell_line == spec_cell_line) %>% filter(gene_type == 'adjusted asRS')
  spec_row_asrt = combine_dat %>% filter(cell_line == spec_cell_line) %>% filter(gene_type == 'adjusted asRT')
  fisher_dat = rbind(c(spec_row_asrs$total_state_genes, spec_row_asrt$total_state_genes), c(spec_row_asrs$total_classified_genes - spec_row_asrs$total_state_genes, spec_row_asrt$total_classified_genes - spec_row_asrt$total_state_genes))
  fisher_results = fisher.test(fisher_dat)
  OR_vec[i] = fisher_results$estimate
  p_vec[i] = fisher_results$p.value
  
}

stats_df = data.frame(cell_line = CELL_LINES, OR = OR_vec, p_value = p_vec)
stats_df$stars.pvalue = stars.pval(stats_df$p_value)
stats_df = left_join(stats_df, combine_dat %>% group_by(cell_line) %>% summarize(max_prev = max(adjusted_prev)), by = 'cell_line') # for y-coordinate plotting 

#plot_dat2 = rbind(bru_plot_dat, intronic_dat %>% select(names(bru_plot_dat)))

#######
#### PLOTTING
theme_set(theme_classic())
#pal = c('asRS' = '#875f9a' , 'adjusted asRT' = '#48d1cc', 'mixed' = '#cbc0d3')
#pal = c('adjusted asRS' = '#875f9a' , 'adjusted asRT' = '#48d1cc')

pal = c('stability' = '#875f9a' , 'transcriptional reg.' = '#48d1cc')
combine_dat$gene_type[combine_dat$gene_type == 'adjusted asRS'] = 'stability'
combine_dat$gene_type[combine_dat$gene_type == 'adjusted asRT'] = 'transcriptional reg.'



cell_line_order = combine_dat %>% filter(gene_type == 'stability') %>% arrange(desc(adjusted_prev))
combine_dat$cell_line = factor(combine_dat$cell_line, levels = cell_line_order$cell_line)

p3 = ggplot(combine_dat) +   geom_linerange(aes(x = cell_line, ymin = 0, ymax = adjusted_prev, group = gene_type), position = position_dodge(width = 0.75))+
  geom_point(pch = 21, color = 'black', size = 6, aes(x = cell_line, y = adjusted_prev, fill = gene_type, group = gene_type), position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pal, breaks = c('stability', 'transcriptional reg.'), name = 'Gene Type') + geom_text(data = stats_df, aes(label = stars.pvalue, x = cell_line, y = max_prev), size = 3, nudge_y = 0.04) + ylab('Prevalence')
combine_dat %>% group_by(cell_line) %>% summarize(combine_prev = sum(adjusted_prev)) %>% arrange(combine_prev)

out = left_join(combine_dat, stats_df, by = 'cell_line') %>% as.data.frame() %>% select(-max_prev) %>% rename(prevalence=adjusted_prev)
write.table(out, file = '../../../SOURCE_DATA/MAIN_FIGS/fig1d.txt', sep = '\t', row.names = F, quote = F)

pdf('../figures/asrs_asrt_lollipop_plot_prev_comparison_add_intronic_ASE_0h.pdf', width = 8, height = 5)
p1 + theme(axis.text.x = element_text(angle = 30, hjust=1)) + theme(legend.position = 'bottom')
p2 + theme_classic(base_size = 13) + theme(axis.text.x = element_text(angle = 30, hjust=1)) + theme(legend.position = 'bottom')
p3 + theme_classic(base_size = 13) + theme(axis.text.x = element_text(angle = 30, hjust=1)) + theme(legend.position = 'bottom')
p4 + theme_classic(base_size = 13) + theme(axis.text.x = element_text(angle = 30, hjust=1)) #+ theme(legend.position = 'bottom')
dev.off()


pdf('../MANUSCRIPT_FIGURES/fig1/fig1d_barplot.pdf', width = 6, height = 3)
p3  + theme(axis.text=element_text(size=6),
            axis.title=element_text(size=8)) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) + theme(legend.position = 'top') + xlab('Cell line')

dev.off()


