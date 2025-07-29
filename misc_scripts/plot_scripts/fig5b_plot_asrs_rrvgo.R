rm(list=ls())
library(data.table)
library(dplyr)
library(rrvgo)
library(org.Hs.eg.db)
library(unikn)
library(sjmisc)
library(forcats)
library(geomtextpath)
library(ggrepel)
library(RColorBrewer)


MIN_OCCUR = 5
SIM_THRESHOLD = 0.7

out_fn = paste0('../GO_RESULTS/REMOVE_5_CELL_LINES_rrvgo_results_min_occur', MIN_OCCUR, '_SIM_THRESHOLD_', SIM_THRESHOLD, '_GENE_LENGTH_RPKM_MATCHED_BKG.Rdata')


load(out_fn)

color_values = usecol(pal_unikn_pref)
color_values = c(color_values, "#FB9A99", "#8DA0CB") # "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99","#386CB0","#F0027F",  "#BF5B17", "#1B9E77", "#7570B3")

n <- 100
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color_values = c(color_values, setdiff(col_vector, color_values))



## SET VARIABLES TO TEST
asrs_reducedTerms = left_join(asrs_reducedTerms_size %>% dplyr::select(-score), asrs_reducedTerms %>% dplyr::select(term, score), by = 'term')
MIN_TERMS_PER_CLUSTER = 2
reducedTerms = asrs_reducedTerms #%>% filter(parentSimScore >= 0.55)
simMatrix = asrs_simMatrix



# only plot clusters with at least n members
bigClusters = table(reducedTerms$parentTerm) %>% as.data.frame() %>% filter(Freq >= MIN_TERMS_PER_CLUSTER)
bigClusters %>% arrange(desc(Freq))


reducedTerms = filter(reducedTerms, parentTerm %in% bigClusters$Var1)
simMatrix = simMatrix[reducedTerms$go,reducedTerms$go]
reducedTerms$parentTerm = word_wrap(reducedTerms$parentTerm, wrap = 15)

## get p-value scores

write.table(reducedTerms %>% filter(score<=25) %>% mutate(parentTerm = gsub("\n", " ", parentTerm)), file = '../../../SOURCE_DATA/MAIN_FIGS/fig5b.txt', sep = '\t', row.names = F, quote = F)


# plot basics
p0 = ggplot(data = reducedTerms %>% filter(score<=25)) + geom_linerange(aes(x = parentTerm, ymin = 0, ymax = score, group = term), 
                                                                        position = position_dodge(width = 0.5)) +
  geom_point(pch = 21, color = 'black', size = 4.5, aes(x = parentTerm,y = score, group = term, fill = parentTerm),
             position = position_dodge(width = 0.5)) 


label_pos = 1.5

p1 = p0 + scale_fill_manual(values = color_values) + theme_bw() + ylab("-log10(p)") + theme(legend.position = 'none')#+ 


p2 = p1 + scale_y_continuous(limits = c(0, 11))
p3 = p1 +scale_y_continuous(limits = c(20, 25), breaks = c(20,22.5, 25)) + xlab("") +  theme(axis.text.x = element_blank())


library(cowplot)


final_out = plot_grid(p3 , p2 , nrow = 2, rel_heights = c(5/11,1))

terms_to_plot = c('innate immune response', 'adaptive immune response', 'T cell activation', 'cellular response to interferon-gamma', 'positive regulation of viral genome replication')

p0 + geom_label_repel(data = reducedTerms %>% filter(term %in% terms_to_plot), aes(x = parentTerm, y = score, label = term, group = term),  position = position_dodge(width = 0.5), min.segment.length = 0)





pdf('../figures/LINEAR_custom_asrs_rrvgo_circleplot.pdf', height = 8.5, width = 11)
final_out
dev.off()

pdf('../MANUSCRIPT_FIGURES/fig5/fig5b.pdf', height = 8.5, width = 11)
final_out
dev.off()



### NOW DO THE SAME FOR ASRT


## SET VARIABLES TO TEST
asrt_reducedTerms = left_join(asrt_reducedTerms_size %>% dplyr::select(-score), asrt_reducedTerms %>% dplyr::select(term, score), by = 'term')
MIN_TERMS_PER_CLUSTER = 2
reducedTerms = asrt_reducedTerms #%>% filter(parentSimScore >= 0.55)
simMatrix = asrt_simMatrix

#lostTerms = asrs_reducedTerms_size %>% filter(parentSimScore < 0.55)


# only plot clusters with at least n members
bigClusters = table(reducedTerms$parentTerm) %>% as.data.frame() %>% filter(Freq >= MIN_TERMS_PER_CLUSTER)
bigClusters %>% arrange(desc(Freq))


reducedTerms = filter(reducedTerms, parentTerm %in% bigClusters$Var1)
simMatrix = simMatrix[reducedTerms$go,reducedTerms$go]

reducedTerms$parentTerm = word_wrap(reducedTerms$parentTerm, wrap = 15)





# plot basics
p0 = ggplot(data = reducedTerms) + geom_linerange(aes(x = parentTerm, ymin = 0, ymax = score, group = term), 
                                                  position = position_dodge(width = 0.5)) +
  geom_point(pch = 21, color = 'black', size = 4.5, aes(x = parentTerm,y = score, group = term, fill = parentTerm),
             position = position_dodge(width = 0.5)) 


label_pos = 1.5

p1 = p0 + scale_fill_manual(values = color_values) + theme_bw() + 
  geom_textvline(xintercept = label_pos, label = "-log10(p)", linetype = 'dashed', hjust = 1, vjust = 1.1) + coord_curvedpolar()



MAX_LIMIT = 25
pdf('../figures/custom_asrt_rrvgo_circleplot.pdf', height = 8.5, width = 8.5)
p1 +  scale_y_continuous(
  limits = c(-2,   MAX_LIMIT),
  expand = c(0, 0),
  breaks = seq(0,MAX_LIMIT,2)) + ylab("-log10(p)")
dev.off()

