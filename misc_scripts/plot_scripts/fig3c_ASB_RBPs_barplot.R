rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(geomtextpath)
library(ggbrace)



#setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts/')

all_snv_fn = '../bed_files/all_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed'
all_snvs = fread(all_snv_fn)
names(all_snvs) = c('chr', 'pos0', 'pos1', 'gene_name', 'gene_state', 'strand', 'alleles', 'cell_line', 'region')

## remove the polyploid cell lines
all_snvs = filter(all_snvs, !cell_line %in% c('PC-3', 'PC-9', 'Caco-2'))


all_asrs_snv_dat = all_snvs %>% filter(gene_state %in% c('state1', 'state2'))
all_transcrip_reg_snv_dat = all_snvs %>% filter(gene_state == 'state3')


#ASB_fn = '../results/match_ASB_and_ASRS_transcripreg_SNV.txt' # out by ASB_snv_enrichment.R
ASB_fn = '../results/ASB_overlapping_all_ASRS_and_ctrl_prop'
ASB = fread(ASB_fn)

###
ASB_interest_region_fn = '../../../ASB/beapr_K562_and_HepG2_combined_ASB_hg38.CUP_removed.txt'
ASB_interest_region_df = fread(ASB_interest_region_fn)


### combine ASB and asRS SNV info
plot_dat = inner_join(all_asrs_snv_dat %>% mutate(id = paste0(chr, ':', pos1, ':', strand)), ASB_interest_region_df %>% mutate(id = paste0(hg38_chr, ':', hg38_coordinate, ':', strand)), by = c('id', 'strand')) %>% tidyr::separate(alleles, into = c('ref', 'alt'), sep = '>')


plot_dat2 =  left_join(plot_dat %>% select(id, ref, alt, RBP) %>% unique() %>% count(RBP) %>% rename(asrs_number_RBP = n), ASB_interest_region_df %>% count(RBP) %>% rename(total_number_RBP = n), by = 'RBP')
plot_dat2$prop = plot_dat2$asrs_number_RBP/plot_dat2$total_number_RBP


RBP_func_fn = '../van_nostrand_supp1.txt'
RBP_func = fread(RBP_func_fn)
RBP_func = rename(RBP_func, RBP = V1)

stab_RBP = RBP_func %>% filter(RBP %in% plot_dat$RBP) %>% filter(`RNA stability & decay` == 1)
plot_dat$color = 'white'
plot_dat$color[plot_dat$RBP %in% stab_RBP$RBP] = "purple"



## remove intronic ASB sites
library(tidyverse)


plot_dat2$prop = signif(plot_dat2$asrs_number_RBP/plot_dat2$total_number_RBP,2)
#plot_dat2 = left_join(plot_dat2, plot_dat %>% select(RBP, color), by = 'RBP')
plot_dat2$color = 'white'
plot_dat2$color[plot_dat2$RBP %in% stab_RBP$RBP] = "purple"
#plot_dat2$RBP = reorder(plot_dat2$RBP, plot_dat2$RBP, length)
#plot_dat2$prop.overlap = factor(plot_dat2$prop.overlap, levels =  plot_dat2$prop.overlap[order(plot_dat2$prop.overlap, decreasing = TRUE)])
#plot_dat2 = as.data.frame(plot_dat2)
plot_dat2$perc = plot_dat2$prop * 100

plot_dat2 = plot_dat2 %>% rename(`RBP type` = color)
plot_dat2$`RBP type`[plot_dat2$`RBP type` == 'purple']  = 'Stability-related'
plot_dat2$`RBP type`[plot_dat2$`RBP type` == 'white']  = 'Other'
plot_dat2$`RBP type` = factor(plot_dat2$`RBP type`, levels = c('Stability-related', 'Other'))



write.table(plot_dat2, file = '../../../SOURCE_DATA/MAIN_FIGS/fig3c.txt', sep = '\t', row.names = F, quote = F)

write.table(plot_dat2, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig8a.txt', sep = '\t', row.names = F, quote = F)

p = ggplot(data = plot_dat2, aes(x = reorder(RBP, perc),  y = perc, fill = `RBP type`, label = asrs_number_RBP)) + geom_bar(stat = 'identity', color = 'white') + theme_classic() + coord_flip() + scale_fill_manual(values = c("#875f9a","#80B1D3")) + theme(legend.position="none") + ylab('% of ASB sites') + xlab('RBP') + geom_text(nudge_y = 0.1, size = 2)



plot_dat2$RBP = paste0(plot_dat2$RBP, "\n", plot_dat2$perc, '%, ', plot_dat2$asrs_number_RBP)


brk <- ggplot_build(p0)$layout$panel_params[[1]]$y$breaks
#brk <- brk[-c(1, length(brk))]
brk = c(0:15)
#p0 + 
#  annotate('text', x = 13.2, y = brk, label = as.character(brk)) + coord_polar()




p0 = ggplot(plot_dat2 %>% arrange(desc(perc)) %>% head(n = 15 )) + geom_hline(aes(yintercept = y), data.frame(y = c(0:13)), color = 'lightgrey') +
  geom_col(aes(x = reorder(RBP, perc), y = perc, fill = `RBP type`), position = 'dodge2', alpha = 0.45, color = 'black') +
  geom_point(aes(x = reorder(RBP, perc), y = asrs_number_RBP, fill = `RBP type`,  pch = `RBP type`), size = 3, color = 'black') +
  geom_segment(aes(x = reorder(RBP, perc), y = 0, xend = reorder(RBP, perc), yend = 13), linetype = 'dashed', color = 'gray12') + 
  #coord_curvedpolar()  +
  coord_polar() + 
  scale_fill_manual(values = c("#875f9a","#80B1D3")) +   scale_color_manual(values = c("#875f9a","#80B1D3")) + theme_bw() + scale_shape_manual(values = c(21, 24)) #+ #geom_brace(aes(c(4,7), c(4, 4.5)), inherit.data=F, rotate = 90, labelrotate=90)



p1 = p0 + xlab("") + ylab("")   + scale_y_continuous(
  limits = c(-5, 13),
  expand = c(0, 0),
  breaks = NULL
) 

# + theme(legend.position = "bottom")


out_fn = '../figures/ASB_RBP_barplot.pdf'
out_fn2 = '../figures/top_10_perc_ASB_circle_plot.pdf'

pdf(out_fn, height = 14)
##p0
p
dev.off()

pdf(out_fn2)
p1
dev.off()


pdf('../MANUSCRIPT_FIGURES/fig3/fig3c.pdf', width = 7
    , height = 6.5)
p1 + theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10)) + theme(legend.position = 'bottom')
dev.off()


pdf('../MANUSCRIPT_FIGURES/extended_fig8/extended_fig8a.pdf', width = 3.3, height = 7)
p + theme(axis.text=element_text(size=6),
        axis.title=element_text(size=8)) 
dev.off()

