rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexUpset)
library(ComplexHeatmap)
library(RColorBrewer)


setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')

number_upset_intersections = 30

all_snvs_fn = '../bed_files/all_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed'

all_snvs_dat = fread(all_snvs_fn)

names(all_snvs_dat) = c('chrom', 'pos0', 'pos1', 'gene_name', 'gene_state', 'gene_strand', 'type_snp', 'cell_line', 'genomic_region')
all_snvs_dat = all_snvs_dat %>% filter(!cell_line %in% c('PC-9', 'PC-3', 'Caco-2'))

# are there any genes that are testable across all cell lines?
all_snvs_dat %>% select(cell_line, gene_name) %>% unique() %>% count(gene_name) %>% arrange(desc(n)) # the most is that a gene is testable in 8 cell lines


all_intronic_res = fread('../bed_files/all_snvs_intronic_0h_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed')
names(all_intronic_res) = c('chrom', 'pos0', 'pos1', 'gene_name', 'gene_state', 'gene_strand', 'type_snp', 'cell_line', 'genomic_region')
all_intronic_res$gene_state[all_intronic_res$gene_state == 1] = 'intronic_asRT'

all_snvs_dat = rbind(all_snvs_dat, all_intronic_res)



all_snvs_dat = all_snvs_dat %>% mutate(id = paste0(chrom, ':', pos1, ':', gene_strand, ':', type_snp))

# 
all_snvs_dat = all_snvs_dat %>% filter(!cell_line %in% c('PC-9', 'PC-3', 'Caco-2'))

states_of_interest = c('state1', 'state2')
PLOT_PREFIX_TITLE = 'asRS'
PLOT_COLOR = '#875f9a'

make_upset = function(states_of_interest, PLOT_PREFIX_TITLE, PLOT_COLOR){
  all_gene_list = list()
  all_snv_list = list()
  
  cell_lines = unique(all_snvs_dat$cell_line)
  for (i in 1:length(cell_lines)){
    interest_dat = all_snvs_dat %>% filter(cell_line == cell_lines[i]) %>% filter(gene_state %in% states_of_interest)
    all_gene_list[[i]] = unique(interest_dat$gene_name)
    all_snv_list[[i]] = unique(interest_dat$id)
  }
  names(all_gene_list) = cell_lines
  names(all_snv_list) = cell_lines
  
  gene_mat = list_to_matrix(all_gene_list)
  gene_tbl = as_tibble(gene_mat)
  
  snv_mat = list_to_matrix(all_snv_list)
  snv_tbl = as_tibble(snv_mat)
  p1 = upset(gene_tbl, names(gene_tbl), n_intersections = number_upset_intersections, name = paste0(PLOT_PREFIX_TITLE, " genes"), 
             base_annotations = list('Intersection size' = intersection_size(counts = FALSE, mapping = aes(fill = 'bars_color'), width = 0.6) +
                                       scale_fill_manual(values = c('bars_color' = PLOT_COLOR), guide = 'none')),
             width_ratio = 0.2, height_ratio = 1, themes=upset_default_themes(text=element_text(size=7)), matrix=(intersection_matrix(geom=geom_point(size=1.2))),   set_sizes=(upset_set_size()
                                                                                                                                                              + theme(axis.text.x=element_text(angle=90))))
  
  p2 = upset(snv_tbl, names(snv_tbl), n_intersections = number_upset_intersections, name = paste0(PLOT_PREFIX_TITLE, " SNVs"), 
             base_annotations = list('Intersection size' = intersection_size(counts = FALSE, mapping = aes(fill = 'bars_color'), width = 0.6) +
                                       scale_fill_manual(values = c('bars_color' = PLOT_COLOR), guide = 'none')),
             width_ratio = 0.2, height_ratio = 1, themes=upset_default_themes(text=element_text(size=7)), matrix=(intersection_matrix(geom=geom_point(size=1.2))),   set_sizes=(upset_set_size()
                                                                                                                                                                      + theme(axis.text.x=element_text(angle=90))))
  out = list(p1, p2, all_gene_list, all_snv_list)
  return(out)
  #print(p1)
  #print(p2)
}


#number_upset_intersections = 50

pdf("../figures/asrs_asrt_gene_snv_upset_plot.pdf", width = 10)
make_upset(c("state1", "state2"), "asRS", '#875f9a')
make_upset(c("intronic_asRT", "asRT"), "asRT" , '#48D1CC')
make_upset(c("state0", "state1", "state2", "state3", "state4", "state5", "state6"), "all classified", 'darkblue')
dev.off()



final_asrs_gene_upset_fn = '../MANUSCRIPT_FIGURES/extended_fig5/extended_fig5a.pdf'
final_het_tested_snv_upset_fn = '../MANUSCRIPT_FIGURES/extended_fig5/extended_fig5b.pdf'
final_het_tested_gene_fn = '../MANUSCRIPT_FIGURES/extended_fig5/extended_fig5d.pdf'
final_asrs_snv_upset_fn = '../MANUSCRIPT_FIGURES/extended_fig5/extended_fig5e.pdf'
final_asrt_gene_upset_fn = '../MANUSCRIPT_FIGURES/extended_fig6/extended_fig6a.pdf'
final_asrt_snv_upset_fn = '../MANUSCRIPT_FIGURES/extended_fig6/extended_fig6b.pdf'

pdf(final_asrs_gene_upset_fn, width = 3.5, height = 2.5)
make_upset(c("state1", "state2"), "asRS", '#875f9a')[[1]]
dev.off()

pdf(final_asrs_snv_upset_fn, width = 3.5, height = 2.5)
make_upset(c("state1", "state2"), "asRS", '#875f9a')[[2]]
dev.off()

pdf(final_het_tested_gene_fn, width = 3.5, height = 2.5)
make_upset(c("state0", "state1", "state2", "state3", "state4", "state5", "state6"), "all classified", 'darkblue')[[1]]
dev.off()

pdf(final_het_tested_snv_upset_fn, width = 3.5, height = 2.5)
make_upset(c("state0", "state1", "state2", "state3", "state4", "state5", "state6"), "all classified", 'darkblue')[[2]]
dev.off()

pdf(final_asrt_gene_upset_fn, width = 3.5, height = 2.5)
make_upset(c("intronic_asRT", "asRT"), "asRT" , '#48D1CC')[[1]]
dev.off()

pdf(final_asrt_snv_upset_fn, width = 3.5, height = 2.5)
make_upset(c("intronic_asRT", "asRT"), "asRT" , '#48D1CC')[[2]]
dev.off()





het_gene = make_upset(c("state0", "state1", "state2", "state3", "state4", "state5", "state6"), "all classified", 'darkblue')[[3]]
het_snv = make_upset(c("state0", "state1", "state2", "state3", "state4", "state5", "state6"), "all classified", 'darkblue')[[4]]


capture.output(het_gene, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig5d.txt')
capture.output(het_snv, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig5b.txt')


