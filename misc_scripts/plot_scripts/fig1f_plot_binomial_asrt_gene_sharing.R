library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gtools)

out = fread('../results/asrt_plus_adjusted_asrt_celltype_sharing_info.combos.txt')

out = filter(out, !celltype1_vec %in% c('PC-9', 'PC-3', 'Caco-2')) %>% filter(!celltype2_vec %in% c('PC-9', 'PC-3', 'Caco-2'))

## plot only binom pval significant results
sig_binom = out %>% arrange(desc(common_test_vec)) %>% filter(binom_pval < 0.05) #%>% filter(prop_overlap > expected_prop)

sig_binom$pstar = stars.pval(sig_binom$binom_p)
#sig_binom$label = paste0(sig_binom$celltype1_vec, '/', sig_binom$celltype2_vec, "\n", sig_binom$pstar, " (", sig_binom$common_test_vec, ")")
sig_binom$label =  paste0(sig_binom$celltype1_vec, '/', sig_binom$celltype2_vec)

plot_dat = sig_binom %>% select(label,  expected_prop, prop_overlap)
m_dat = melt(plot_dat, id.vars = c("label"))
stars_dat = melt(sig_binom %>% select(label, pstar))
names(m_dat) = c('label',  'prop.type', 'Proportion')
m_dat$prop.type = as.character(m_dat$prop.type)

m_dat = rename(m_dat, `Proportion Type` = prop.type)
m_dat$`Proportion Type`[m_dat$`Proportion Type` == 'expected_prop'] <- "Background Expectation"
m_dat$`Proportion Type`[m_dat$`Proportion Type` == 'prop_overlap'] <- "Actual"

ordering = m_dat %>% filter(`Proportion Type` == 'Actual') %>% arrange(desc(Proportion))
m_dat$label = factor(m_dat$label, levels = ordering$label)
  #m_dat$label = gsub('_', '\n', m_dat$label)
#m_dat$pstar = stars.pval(m_dat$binom_p)

all_dat = left_join(m_dat, stars_dat, by = 'label')

ordering = all_dat %>% filter(`Proportion Type` == 'Actual') %>% arrange(desc(Proportion))
all_dat$label = factor(all_dat$label, levels = ordering$label)


print(nrow(all_dat))

write.table(sig_binom %>% select(-hyper_pval_vec), file = '../../../SOURCE_DATA/MAIN_FIGS/fig1f.txt', sep = '\t', row.names = F, quote = F)

p = ggplot(all_dat, aes(fill = `Proportion Type`, y = Proportion, x = label)) + geom_text(data = all_dat %>% filter(`Proportion Type` == 'Actual'), aes(label = pstar, x = label, y = Proportion), nudge_y = 0.0005, size = 3) + 
    geom_bar(position="dodge", stat="identity", width = 0.7, col = 'black', size = 0.25) + theme_classic() + scale_fill_manual(values = c("#48D1CC", "#d9d9d9" )) + xlab("") + theme(axis.text.x = element_text(angle = 90 , vjust = 1, hjust=1))#+ coord_flip() 

pdf('../figures/asrt_celltype_sharing_binom_sig.pdf')
p
dev.off()


pdf('../MANUSCRIPT_FIGURES/fig1/fig1f.pdf', width = 5, height = 4)
p  + theme(axis.text=element_text(size=6),
           axis.title=element_text(size=8)) +
  theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6), legend.position = 'top')  #change legend text font size 
dev.off()

