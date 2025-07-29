library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gtools)
#setwd('~/Desktop/project-gxxiao4/ASE_ASRS/gene_strand_fix/scripts')

out = fread('../results/asrs_variant_celltype_sharing_info.combos.txt')

# remove polyploidy cell lines
out = filter(out, !celltype1_vec %in% c('PC-9', 'PC-3', 'Caco-2')) %>% filter(!celltype2_vec %in% c('PC-9', 'PC-3', 'Caco-2'))

## plot some of the numbers for vis
cutoff = 400

plot_bar = function(cutoff) {
  sig_binom = filter(out, common_test_vec > cutoff)
  sig_binom$label = paste0(sig_binom$celltype1_vec, '_', sig_binom$celltype2_vec)
  plot_dat = sig_binom %>% select(label, expected_prop, prop_overlap)
  m_dat = melt(plot_dat, id.vars = "label")
  names(m_dat) = c('label', 'prop.type', 'proportion')
  m_dat$variable[m_dat$variable == 'prop_overlap'] <- 'actual.overlap.prop'
  m_dat$variable[m_dat$variable == 'expected_prop'] <- 'expected.overlap.prop'
  
  ordering = m_dat %>% filter(prop.type == 'prop_overlap') %>% arrange(desc(proportion))
  m_dat$label = factor(m_dat$label, levels = ordering$label)
  #m_dat$label = gsub('_', '\n', m_dat$label)
  
  p = ggplot(m_dat, aes(fill = prop.type, y = proportion, x = label)) +
    geom_bar(position="dodge", stat="identity", width = 0.7) + coord_flip() + theme_classic() + ggtitle(paste('at least ', cutoff, ' common testable genes'))
  return(p)
}



## plot only binom pval significant results
sig_binom = out %>% arrange(desc(common_test_vec)) %>% filter(binom_pval < 0.05) %>% filter(prop_overlap > expected_prop)
sig_binom$pstar = stars.pval(sig_binom$binom_p)
#sig_binom$label = paste0(sig_binom$celltype1_vec, '/', sig_binom$celltype2_vec, "\n", sig_binom$pstar, " (", sig_binom$common_test_vec, ")")
sig_binom$label =  paste0(sig_binom$celltype1_vec, '/', sig_binom$celltype2_vec)


plot_dat = sig_binom %>% select(label,  expected_prop, prop_overlap)
m_dat = melt(plot_dat, id.vars = c("label"))
stars_dat = melt(sig_binom %>% select(label, pstar))
names(m_dat) = c('label',  'prop.type', 'Proportion')
m_dat$prop.type = as.character(m_dat$prop.type)

m_dat = rename(m_dat, `Proportion Type` = prop.type)
m_dat$`Proportion Type`[m_dat$`Proportion Type` == 'expected_prop'] <- "Expected"
m_dat$`Proportion Type`[m_dat$`Proportion Type` == 'prop_overlap'] <- "Actual"

ordering = m_dat %>% filter(`Proportion Type` == 'Actual') %>% arrange(desc(Proportion))
m_dat$label = factor(m_dat$label, levels = ordering$label)
  #m_dat$label = gsub('_', '\n', m_dat$label)
#m_dat$pstar = stars.pval(m_dat$binom_p)

all_dat = left_join(m_dat, stars_dat, by = 'label')

ordering = all_dat %>% filter(`Proportion Type` == 'Actual') %>% arrange(desc(Proportion))
all_dat$label = factor(all_dat$label, levels = ordering$label)

p = ggplot(all_dat, aes(fill = `Proportion Type`, y = Proportion, x = label)) + geom_text(data = all_dat %>% filter(`Proportion Type` == 'Actual'), aes(label = pstar, x = label, y = Proportion), nudge_y = 0.0015, size = 2.5) + 
    geom_bar(position="dodge", stat="identity", width = 0.7, col = 'black', size = 0.25) + theme_classic() + scale_fill_manual(values = c("#875f9a", "#d9d9d9" )) + xlab("") + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))#+ coord_flip() 

#p = ggplot(all_dat, aes(fill = `Proportion Type`, y = Proportion, x = label)) + geom_text(data = all_dat %>% filter(`Proportion Type` == 'Actual'), aes(label = pstar, x = label, y = Proportion), nudge_y = 0.002, size = 3.5, angle = 90) + 
#    geom_bar(position="dodge", stat="identity", width = 0.7, col = 'black', size = 0.25) + theme_classic() + scale_fill_manual(values = c("#875f9a", "#d9d9d9" )) + xlab("") #+ coord_flip() 


write.table(sig_binom, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig5f.txt', sep = '\t', row.names = F, quote = F)

pdf('../figures/asrs_variant_celltype_sharing_binom_sig.pdf')
p
dev.off()

pdf('../MANUSCRIPT_FIGURES/extended_fig5/extended_fig5f.pdf', height = 4, width = 4)
p + theme(legend.position = 'bottom')
dev.off()











