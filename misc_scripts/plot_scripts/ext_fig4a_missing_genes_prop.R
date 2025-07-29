rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)

prop_dat = fread('../mosdepth_gene_coverage/ALL_CELL_LINES_coverage_genes_labels.txt')

prop_dat = filter(prop_dat, !cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))
### note that missed_genes_intronic_rdd + missed_genes_from_no_rdd != missed_genes_from_no_rdd because we take the intersection across all cell lines (so some genes may be in one category in one sample but in another category at another timepoint)


prop_dat$prop_missing = prop_dat$missed_genes/prop_dat$coverage_pass_genes
prop_dat$prop_missing_protein_coding = prop_dat$missing_genes_protein_coding/prop_dat$coverage_pass_protein_coding_genes


m_dat_all = prop_dat %>% dplyr::mutate(genes_with_genic_het_snvs = coverage_pass_genes - missed_genes) %>% dplyr::rename(genes_without_genic_het_snvs = missed_genes) %>% dplyr::select(cell_line, genes_with_genic_het_snvs, missed_genes_intronic_rdd, missed_genes_from_no_rdd) %>% melt()
#m_dat_all = prop_dat %>% mutate(genes_with_genic_het_snvs = coverage_pass_genes - missed_genes) %>% rename(genes_without_genic_het_snvs = missed_genes) %>% select(cell_line, genes_with_genic_het_snvs, genes_without_genic_het_snvs) %>% melt()
#m_dat_pc = prop_dat  %>% select(cell_line, coverage_pass_protein_coding_genes, missing_genes_protein_coding) %>% melt()
names(m_dat_all) = c('cell_line', 'gene_type', 'number_genes')


m_dat_all$gene_type = as.character(m_dat_all$gene_type)
m_dat_all$gene_type[m_dat_all$gene_type == 'genes_with_genic_het_snvs'] = '>1 het mRNA SNV'
m_dat_all$gene_type[m_dat_all$gene_type == 'missed_genes_intronic_rdd'] = '>1 intronic het SNV'
m_dat_all$gene_type[m_dat_all$gene_type == 'missed_genes_from_no_rdd'] = '0 het SNVs'

cell_line_order = m_dat_all %>% group_by(cell_line) %>% summarize(number_genes = sum(number_genes)) %>% arrange(desc(number_genes))
m_dat_all$cell_line = factor(m_dat_all$cell_line, levels = cell_line_order$cell_line)

p = ggplot(m_dat_all, aes(cell_line, y = number_genes, fill = gene_type)) + geom_bar(position = 'stack', stat = 'identity', color = 'black') + coord_flip() + theme_classic() + scale_fill_brewer(palette = 'Greens')

## now get prop of missing genes with intronic rdds across all timepoints
prop_dat$prop_missed_genes_intronic_rdd = round(prop_dat$missed_genes_intronic_rdd/prop_dat$missed_genes,2)
p2 = ggplot(prop_dat, aes(cell_line, y = missed_genes_intronic_rdd, label = prop_missed_genes_intronic_rdd )) + 
  geom_bar(stat = 'identity', fill = 'steelblue', color = 'black', alpha = 0.67) + coord_flip() + 
  theme_bw() + geom_label()  + ylab('Number of missed genes with intronic het. snvs') + xlab('Cell line')


prop_dat = prop_dat %>% arrange(desc(missed_genes_intronic_rdd))
prop_dat$cell_line = factor(prop_dat$cell_line, levels = prop_dat$cell_line)

## or do one line plot and one bar plot
part1 = ggplot(prop_dat, aes(cell_line, y = missed_genes_intronic_rdd )) + 
  geom_bar(stat = 'identity', fill = 'steelblue', color = 'black', alpha = 0.67)  + 
  theme_classic() + ylab('Number of missed genes with intronic het. snvs') + xlab('Cell line') 

part2 = ggplot(prop_dat, aes(x = cell_line, y = prop_missed_genes_intronic_rdd, group = 1 )) + 
  geom_line(stat = 'identity', color = 'red')  + 
  geom_point(fill = 'black', color = 'black', size = 2, pch = 21) + 
  theme_classic() + ylab('Proportion') + xlab('Cell line') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p3 = plot_grid(part2, part1, nrow = 2, rel_heights = c(0.2, 0.8))




write.table(m_dat_all, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig4a.txt', sep = '\t', row.names = F, quote = F)

pdf('../figures/prop_of_genes_missing_genic_het_snvs.pdf')
p + theme(legend.position = 'bottom') + xlab('Cell line') + ylab('Number of genes')
p2
p3
dev.off()


pdf('../MANUSCRIPT_FIGURES/extended_fig4/extended_fig4a.pdf', height = 5, width = 5)
p + theme(legend.position = 'bottom') + xlab('Cell line') + ylab('Number of genes')
dev.off()



