# script to filter parent/ancestor GO terms from inputted list of GO terms enriched
# to retain the more specific child terms

#source("http://bioconductor.org/biocLite.R")
#A set of annotation maps describing the entire Gene Ontology
#biocLite("GO.db")
rm(list=ls())
library(GO.db)
library(data.table)
library(dplyr)
library('unikn')

args = commandArgs(trailingOnly=TRUE)
#plot_title = ""
#input_GO_file = '../results/GO_enrichment_of_differentially_edited_genes.txt'

#setwd('~/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')

#trait = 'all_asrt'
#trait = 'all_asrs_asrt_bg_ensembl_genes'
#input_GO_file = paste0("GO_queries/", trait, ".txt_GO_results")
#input_GO_file = '../GO_RESULTS/ALL_state1_state2_GO_enrichment_res.txt'
#input_GO_file = '../GO_RESULTS/ALL_state3_bg_state0_GO_enrichment_res.txt'
#input_GO_file = '../GO_RESULTS/ALL_state1_state2_bg_state3_GO_enrichment_res.txt'


min_occurences = 5
FDR_cutoff = 0.05
term_cutoff = 60 # number of rows to include in output figure
plot_col = "#8DA0CB"



make_GO_plot = function(input_GO_file, trait){
  
  plot_title = paste0(trait, ' ', 'min occur = ', min_occurences)
  output_pdf = paste0("../figures/", trait, "_min_occur_", min_occurences, '_top_', term_cutoff, "_terms.pdf")
  #plot_title = "BrainGVEX and CMC WGCNA disease module genes"
  #input_GO_file="../results/PE_CMC_module_all_cohort_overlapping_ENSIDs_111_GO_results.txt"
  #output_pdf=paste0("../plots/PE_CMC_module_all_cohort_overlapping_ENSIDs_111_GO_results_min",min_occurences,"_occur_FDR_",FDR_cutoff,".filtered_parent.pdf")
  #plot_title = paste("cluster_1")
  #input_GO_file = "../mfuzz_results/data.s_cluster_comparison/cluster_size_6/GO_enrichment_of_cluster_1"
  #output_pdf = paste0("../figures/differential_editing_GO_plot.filtered_parent_plot.min.occur", min_occurences, ".pdf")
  
  GO_table = read.table(input_GO_file,sep="\t",stringsAsFactors=FALSE,header=TRUE,quote="\"",check.names=FALSE, fill=T)

  #keep only FDR significant GO categories
  plotting_table = GO_table[GO_table$FDR <= FDR_cutoff,]
  #keep only GO categories with a min number of genes in the category
  plotting_table = plotting_table[plotting_table$ocurrences_in_query >= min_occurences,]
  #plotting_table = plotting_table %>% filter(namespace_1003 == 'biological_process')
  #plotting_table = plotting_table %>% arrange(FDR)
  
  plotting_table = head(plotting_table, n = term_cutoff)
  
  #plotting_table = head(plotting_table %>% filter(namespace_1003 == 'biological_process'), n = term_cutoff)
  print(nrow(plotting_table))
  plotting_table = head(plotting_table, n = term_cutoff)
  plotting_table$name_1006 = factor(plotting_table$name_1006,levels=rev(plotting_table$name_1006))
  
  #get the - value of the logp value
  
  #get the - value of the logp value
  
  zero_values = which(plotting_table$enrichment_p_values == 0)
  if (length(zero_values)>=1){
    plotting_table[zero_values,"enrichment_p_values"] = 1e-1000
    plotting_table[zero_values,"FDR"] = 1e-1000
  }
  plotting_table$enrichment_p_values = -1 * log10(plotting_table$enrichment_p_values)
  plotting_table$neg_log_FDR = -log10(plotting_table$FDR)
  
  #plot_title = 'GO enrichment on asRS genes (min. occur = 10)'
  library(ggplot2)
  #library(viridis)
  #theme_set(theme_bw())
  #theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #theme_update(plot.title = element_text(size=6), axis.title.y = element_text(size=5), axis.title.x = element_text(size=5), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4))
  #modify some legend attributes
  #theme_update(legend.key.size = unit(3,'mm'), legend.text = element_text(size=3), legend.title = element_text(size=4), legend.key = element_rect(color=NA), legend.margin = unit(0,'mm'))
  sig_value = -log10(0.05)
  pdf(output_pdf,width=15,height=13)
  
  #main = ggplot(plotting_table %>% filter(namespace_1003 != 'cellular_component'), aes(x=name_1006,y=(enrichment_p_values), fill = label)) + geom_bar(stat='identity',aes(fill = label), color = 'black', alpha = 0.67) + scale_fill_manual(values = color_pal) + coord_flip()
  main = ggplot(plotting_table, aes(x=name_1006,y=(enrichment_p_values), fill = label)) + geom_bar(stat='identity',aes(fill = label), color = 'black', alpha = 0.67) #+ scale_fill_manual(values = color_pal) + coord_flip()
  
  main = ggplot(plotting_table, aes(x=name_1006,y=(enrichment_p_values))) + geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip()
  
  #ylabel = ylab("-log10 FDR")
  ylabel = ylab("-log10(adj p)")
  title = ggtitle("Gene Ontology Enrichment")
  vline = geom_hline(yintercept = sig_value)
  xlabel = xlab("")
  #p = main+geom+options+coordinates+FDR_cutoff_annotation+options+ylabel+xlabel+title
  p = main+vline+ylabel+xlabel+title+theme_classic(base_size = 18)
  print(p)
  
  dev.off()
  return(plotting_table)
  
}


#asrt = make_GO_plot('../GO_RESULTS/ALL_state3_bg_state0_GO_enrichment_res.txt', 'asrt')
#asrs = make_GO_plot('../GO_RESULTS/ALL_state1_state2_bg_state0_GO_enrichment_res.txt', 'asrs')
#asrs_asrt_bg = make_GO_plot('../GO_RESULTS/ALL_state1_state2_bg_state3_GO_enrichment_res.txt', 'asrs_asrt_bg')

#asrt_yes = make_GO_plot('../GO_RESULTS/REMOVE_5_CELL_LINES_ALL_state3_bg_GENE_LENGTH_RPKM_CONTROLLED_THROW_OUT_GENES_yesGO_enrichment_res_gene_length_RPKM_controlled.txt', 'asrt_yes')

asrt_yes = make_GO_plot('../GO_RESULTS/REMOVE_5_CELL_LINES_ALL_ASRT_AND_INTRONIC_ASRT_bg_GENE_LENGTH_RPKM_CONTROLLED_THROW_OUT_GENES_yesGO_enrichment_res_gene_length_RPKM_controlled.txt', 'asrt_add_intronic_yes')
asrs_yes = make_GO_plot('../GO_RESULTS/REMOVE_5_CELL_LINES_ALL_state1_state2_bg_GENE_LENGTH_RPKM_CONTROLLED_THROW_OUT_GENES_yesGO_enrichment_res_gene_length_RPKM_controlled.txt', 'asrs_yes')

asrs_yes$gene_type = 'asRS_yes'
asrt_yes$gene_type = 'asRT_yes'

combine_dat_yes = rbind(asrs_yes, asrt_yes)


ggplot(head(asrs_yes %>% mutate(ocurrences_in_query = factor(ocurrences_in_query)), n=15), aes(x=name_1006, fill = ocurrences_in_query, y=(enrichment_p_values))) + geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip() + scale_fill_brewer(palette = 'Reds')



write.table(head(asrs_yes %>% mutate(ocurrences_in_query = factor(ocurrences_in_query)), n=15), file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig10a.txt', sep = '\t', row.names = F, quote = F)

#ggplot(head(asrs_yes,n = 10) , aes(x=name_1006, fill = gene_type, y=(enrichment_p_values))) + geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip()


#asrt_no = make_GO_plot('../GO_RESULTS/ALL_state3_bg_GENE_LENGTH_RPKM_CONTROLLED_THROW_OUT_GENES_noGO_enrichment_res_gene_length_RPKM_controlled.txt', 'asrt_no')
#asrs_no = make_GO_plot('../GO_RESULTS/ALL_state1_state2_bg_GENE_LENGTH_RPKM_CONTROLLED_THROW_OUT_GENES_noGO_enrichment_res_gene_length_RPKM_controlled.txt', 'asrs_no')

#asrs_no$gene_type = 'asRS_no'
#asrt_no$gene_type = 'asRT_no'

#combine_dat_no = rbind(asrs_no, asrt_no)

#ggplot(combine_dat_no , aes(x=name_1006, fill = gene_type, y=(enrichment_p_values))) + geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip()


#asrs_asrt_bkg = make_GO_plot('../GO_RESULTS/ALL_state1_state2_bg_state3_THROW_OUT_GENES_yes_GO_enrichment_res.txt', 'asrs_asrt_yes')
#ggplot(asrs_asrt_bkg , aes(x=name_1006, y=(enrichment_p_values))) + geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip() + theme_bw()

#
#load('../GO_RESULTS/ALL_state1_state2_bg_GENE_LENGTH_RPKM_CONTROLLED_THROW_OUT_GENES_yes_query_GO_table_gene_length_RPKM_controlled.RData')

#test = query_GO_table %>% filter(name_1006 == 'positive regulation of nuclear-transcribed mRNA catabolic process, deadenylation-dependent decay')
#cat('job completed\n')
#query_GO_table %>% filter(name_1006 %like% 'immu') %>% filter(hgnc_symbol %in% test$hgnc_symbol)

#pdf('../figures/GO_barplot_asrs_asrt_throw_out_unmatched_queries_RPKM_gene_length_matched_bkg.pdf', height = 12)
pdf('../figures/GO_barplot_asrs_asrt_add_intronic_throw_out_unmatched_queries_RPKM_gene_length_matched_bkg.pdf', height = 12)
ggplot(combine_dat_yes , aes(x=name_1006, fill = gene_type, y=(enrichment_p_values))) + geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip() + theme_bw()
dev.off()




## MAKE PRETTY FIGS NOW
plot_dat = head(asrs_yes %>% mutate(ocurrences_in_query = factor(ocurrences_in_query)),n=15)
#plot_dat$name_1006 = word_wrap(plot_dat$name_1006,15)

p = ggplot(plot_dat, aes(x=name_1006, fill = ocurrences_in_query, y=(enrichment_p_values))) + 
  geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip() + scale_fill_brewer(palette = 'Reds') + theme_classic(base_size = 18) + xlab('GO term') + ylab("-log10(enrichment p-value)")



write.table(plot_dat, file = '../../../SOURCE_DATA/MAIN_FIGS/fig5a.txt', sep = '\t', row.names = F, quote = F)

pdf('../figures/asrs_remove_unmatched_yes_RPKM_gene_length_bkg.pdf', width = 15)


pdf('../MANUSCRIPT_FIGURES/fig5/fig5a.pdf', width = 15)
p
dev.off()



plot_dat2 = head(asrt_yes %>% mutate(ocurrences_in_query = factor(ocurrences_in_query)),n=20)
#plot_dat$name_1006 = word_wrap(plot_dat$name_1006,15)

p = ggplot(plot_dat2, aes(x=name_1006, fill = ocurrences_in_query, y=(enrichment_p_values))) + 
  geom_bar(stat='identity', color = 'black', alpha = 0.67) + coord_flip() + scale_fill_brewer(palette = 'Blues') + theme_classic(base_size = 18) + xlab('GO term') + ylab("-log10(enrichment p-value)")


#pdf('../figures/asrt_remove_unmatched_yes_RPKM_gene_length_bkg.pdf', width = 15)
pdf('../figures/asrt_add_intronic_remove_unmatched_yes_RPKM_gene_length_bkg.pdf', width = 15)
p 
dev.off()

pdf('../MANUSCRIPT_FIGURES/extended_fig9/extended_fig9b.pdf', width = 15)
p
dev.off()

p 
dev.off()

