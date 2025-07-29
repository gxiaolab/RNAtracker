library(data.table)
library(dplyr)
library(ggplot2)
library(geomtextpath)
#library(ggbrace)
library(ggpubr)



ASB_fn = '../results/ASB_overlapping_all_ASRS_and_ctrl_prop'
ASB_prop = fread(ASB_fn)

ASB_prop = filter(ASB_prop, !celltypes %in% c('Caco-2', 'PC-3', 'PC-9'))


write.table(ASB_prop, file = '../../../SOURCE_DATA/MAIN_FIGS/fig3b.txt', quote = F, sep = '\t', row.names = F)
ASB_prop = ASB_prop %>% select(celltypes, asrs_prop, nonase_prop)
names(ASB_prop) = c('cell_line', 'asRS', 'nonASE')
m_out = melt(ASB_prop, id.vars = "cell_line", value.name = 'prop', variable.name = 'snv_type')
my_comparisons = list(c('asRS', 'nonASE'))
p0 = ggplot(m_out, aes(x = snv_type, y = prop)) +
  geom_boxplot(aes(fill = snv_type), width = 0.5) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.6) + 
  scale_fill_manual(values = c("#875f9a", "#a5a5a5")) + xlab('Variant Type') + ylab('Proportion of variants in ASB site')+ stat_compare_means(comparisons = my_comparisons, paired = TRUE, size = 2.5)
p = p0 + theme(legend.position = "none")




pdf('../figures/ASB_boxplot.pdf', height = 4, width = 4)
p + theme_bw()
dev.off()

pdf('../MANUSCRIPT_FIGURES/fig3/fig3b.pdf',width = 2, height = 2.5)
p + theme_classic() + theme(axis.text=element_text(size=8),
                            axis.title=element_text(size=8)) + guides(fill=guide_legend(title="Variant Type"))  + theme(legend.position="none")
dev.off()
