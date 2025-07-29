library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

load('../results/miRNA_SNP_expressed_miRNAs_prop.RData')

prop_dat = all_res %>% select(cell_line, asrs_prop, ctrl_prop)

all_res$not_asrs = all_res$total_asrs - all_res$total_asrs_hits
all_res$not_ctrl = all_res$total_ctrl - all_res$total_ctrl_hits

OR_vec = rep(0, nrow(all_res))
p_vec = rep(0, nrow(all_res))
for(i in 1:nrow(all_res)){
	row_info = all_res[i,]
	fisher_dat = rbind(c(row_info$total_asrs_hits, row_info$total_ctrl_hits), c(row_info$not_asrs, row_info$not_ctrl))	
  fisher_dat = fisher_dat + 1 # add pseudocount
  res = fisher.test(fisher_dat)
  OR_vec[i] = res$estimate
  p_vec[i] = res$p.value
}

all_res$OR = OR_vec
all_res$p_value = p_vec


# remove the polyploid cell lines
prop_dat = filter(prop_dat, !cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))
write.table(all_res %>% filter(!cell_line %in% c('PC-3', 'PC-9', 'Caco-2')), file = '../../../SOURCE_DATA/MAIN_FIGS/fig3d.txt', sep = '\t', row.names = F, quote = F)


names(prop_dat) = c('celltypes', 'asRS', 'Control')
m_dat = melt(prop_dat)
names(m_dat) = c('celltypes', 'SNV_type', 'prop')

my_comparisons = list(c('asRS', 'Control'))
p0 = ggplot(m_dat, aes(x = SNV_type, y = prop)) + geom_boxplot(aes(fill = SNV_type), width = 0.5) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) + 
  scale_fill_manual(name = 'Variant Type', values = c("#875f9a", "#a5a5a5")) +
  xlab('Variant Type') + ylab('Proportion of SNVs with miRNA target gain/loss effect') + stat_compare_means(comparisons = my_comparisons, paired = TRUE) #+ ggtitle('RNAtracker asRS SNVs')

prop_dat2 = all_res %>% select(cell_line, prop_of_mirna_asrs, prop_of_mirna_ctrl)
names(prop_dat2) = c('celltypes', 'asRS', 'Control')
m_dat2 = melt(prop_dat2)
names(m_dat2) = c('celltypes', 'SNV_type', 'prop')

p1 = ggplot(m_dat2, aes(x = SNV_type, y = prop)) + geom_boxplot(aes(fill = SNV_type), width = 0.5) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) + theme_classic(base_size = 12) + 
  scale_fill_manual(name = 'Variant Type', values = c("#875f9a", "#a5a5a5")) +
  xlab('Variant Type') + ylab('Proportion of miRNA target gain/loss SNPs\nthat overlap asRS or Control variant') + stat_compare_means(comparisons = my_comparisons, paired = TRUE, size = 5) #+ ggtitle('RNAtracker asRS SNVs')


pdf('../figures/mirna_snp_boxplots.pdf')
p0 
p1
dev.off()

write.table(m_dat, file = '../../../SOURCE_DATA/MAIN_FIGS/fig3d.txt', sep = '\t', row.names = F, quote = F)
write.table(m_dat2, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig8b.txt', sep = '\t', row.names = F, quote = F)



pdf('../MANUSCRIPT_FIGURES/fig3/fig3d.pdf', width = 2, height = 2.5)
p0 + theme_classic()  + guides(fill=guide_legend(title="Variant Type"))  + theme(legend.position="none") + theme(axis.text=element_text(size=8),axis.title=element_text(size=8))
dev.off()

pdf('../MANUSCRIPT_FIGURES/extended_fig8/extended_fig8b.pdf', width = 2, height = 2.5)
p1 + theme(legend.position="none") + theme(axis.text=element_text(size=8),axis.title=element_text(size=8))
dev.off()




