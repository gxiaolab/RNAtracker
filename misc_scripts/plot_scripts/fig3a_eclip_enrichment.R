library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(cowplot)


load('../results/REMOVE_5_CELL_LINES_eclip.bedtools.strand.fixed.plot.dat.control.snp.subset.RData')
plot_scatter = function(dat) {
  dat$logp = -log10(dat$p_vec)
  RBP_interest =  dat %>% filter(p_vec < 0.05) %>% filter(OR_vec > 1)
  
  RBP_func_fn = '../van_nostrand_supp1.txt'
  RBP_func = fread(RBP_func_fn)
  RBP_func = rename(RBP_func, RBP_name = V1) 
  RBP_func %>% filter(RBP_name %in% RBP_interest$RBP)
  
  stab_RBP = RBP_func %>% filter(RBP_name %in% RBP_interest$RBP) %>% filter(`RNA stability & decay` == 1)
  #stab_RBP = RBP_func %>% filter(`RNA stability & decay` == 1)
  
  dat$`RBP type` = 'not sig'
  dat$`RBP type`[dat$p_vec < 0.05 & dat$OR > 1] = 'Other'
  dat$`RBP type`[dat$RBP %in% stab_RBP$RBP_name] = "Stability-related"
  color_pal = c("not sig" = "gray", "Other" = '#80B1D3', "Stability-related" = "#875F9A")
  dat$`RBP type` = factor(dat$`RBP type`, levels = c('Stability-related', 'Other', 'not sig'))
  dat = rename(dat, OR = OR_vec)
 
    p = ggplot(dat, aes(y = logp, x = OR)) + geom_hline(yintercept = -log10(0.05), linetype = 'dashed') + geom_vline(xintercept = 1, linetype = 'dashed') +
      geom_point(aes(fill = `RBP type`, shape = `RBP type`, size = `RBP type`),color = "black") + 
      geom_text_repel(size = 5, dat = dat %>% filter(p_vec < 0.05) %>% filter(OR >1) %>% arrange(p_vec), aes(label = RBP, color = `RBP type`), min.segment.length = unit(0, 'lines'), max.overlaps = 100, seed = 123, nudge_y = 0.1, nudge_x = 0.1)  + 
      theme_classic(base_size = 12.5) + ylab('-log10(p)') + scale_color_manual(values = color_pal) + scale_fill_manual(values = color_pal)  + xlab('Odds Ratio') + scale_shape_manual(values = c(21,24,21)) + scale_size_manual(values = c(3,3,2))
    
    
    return(p)
}

p1 = plot_scatter(results)


pdf('../figures/eclip.bedtools.strand.fix.enrichment.results.pdf', width = 8, height = 5)
p1
dev.off()


pdf('../MANUSCRIPT_FIGURES/fig3/fig3a.pdf', width = 7, height = 6)
p1 + theme(legend.position = 'bottom')#+ theme(axis.text=element_text(size=7),
   #        axis.title=element_text(size=7), axis.text.x = element_text(size = 7)) + theme(legend.position = 'bottom')
dev.off()


