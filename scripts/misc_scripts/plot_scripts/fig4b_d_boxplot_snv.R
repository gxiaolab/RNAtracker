rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(ComplexUpset)
library(gtools)
library(RColorBrewer)
library(cowplot)



## PLOT SELECT EXAMPLE SNVs WITH CUSTOM AXIS RANGES



#setwd('~/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/deepseq/RNASEQ/scripts')
#args = commandArgs(trailingOnly=TRUE)
#CELL_LINE = args[1]
#WGS_LABEL = args[2]
#SNV_ID = ""
#gene_name = ""

#CELL_LINE = 'GM12878'
#WGS_LABEL = 'GM12878gDNA_S3_L006'
#SNV_ID = 'chr17:81672520:+'
#gene_name = 'CCDC137'

#CELL_LINE = 'HCT116'
#WGS_LABEL = 'HCT116gDNA_S2_L006'
#SNV_ID = 'chr21:41888378:-'
#gene_name = 'C2CD2'

CELL_LINE = 'MCF-7'
WGS_LABEL = 'MCF7gDNA_S1_L006'
SNV_ID = 'chr2:17781251:+'
gene_name = 'GEN1'





#read_dat = function(WGS_LABEL){
#  in_fn = paste0('../REAL_DATA_RECOVERED_SNVS_10_2/', WGS_LABEL, '.ALL_TIMEPOINTS.RECOVER_SNVs_INCLUDE_ZERO.RData')
 # load(in_fn)
 # return(all_timepoints_complete)
#}

read_dat = function(WGS_LABEL, TIMEPOINT){
  in_fn = paste0('../REAL_DATA_RESULTS/pickle_files_FIX_MAJOR_ALLELE_GATK_VAR_10_2/', WGS_LABEL, '/', TIMEPOINT, '/T10_M2_eiTrue_ARFalse_nsigma1.95996_p0.05_apfdr_bh_FIXED_MAJOR_ALLELE_ASE')
  gene_dat = fread(in_fn)
  gene_dat$TIMEPOINT = TIMEPOINT 
  return(gene_dat)
}



make_fig = function(CELL_LINE, WGS_LABEL, SNV_ID, gene_name){
  
#  ACTD_DAT = read_dat(WGS_LABEL)
  
  #out_data_fn =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV.RData')
  #out_data_fn =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV_ALL_STATE_GENES_PVAL_STRICT_TEST.RData')
  
  out_data_fn =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV_ALL_STATE_GENES.RData')
  load(out_data_fn)
  

  
  timepoints = c('0h', '2h', '8h', '24h')
  all_data = list()
  for (timepoint in timepoints){
    print(timepoint)
    dat = read_dat(WGS_LABEL, timepoint )
    all_data[[timepoint]] = dat
  }
  
  snv_level_dat = do.call("rbind", all_data)  %>% tidyr::separate_rows(SNVs, sep = ';') %>% select(Gene, ASEornonASE, SNVs, TIMEPOINT) %>% filter(SNVs != "")  %>% tidyr::separate(SNVs, into = c('coord', 'region', 'AR', 'ref_allelic_ratio', 'delta_AR', 'sigma', 'major', 'minor', 'p_orig', 'p_adj', 'nraw', 'nnorm'), sep = ',', convert = TRUE) %>% tidyr::separate(coord, into = c('chr', 'pos', 'strand'), sep = ':')
  snv_level_dat = mutate(snv_level_dat, id = paste0(chr, ':', pos, ':', strand))
  
  all_dat = all_results %>% filter(validation_status) %>% filter(gene_state %in% c('state1', 'state2'))
  clean_ids = gsub(':[A-Z]>[A-Z]', '', all_dat$id)
  nonASE_0h = snv_level_dat %>% filter(id %in% clean_ids) %>% filter(TIMEPOINT == '0h') %>% filter(p_adj > 0.05)
  ASE_later = snv_level_dat %>% filter(id %in% nonASE_0h$id) %>% filter(TIMEPOINT != '0h') %>% filter(p_adj < 0.05)
  top_hits = ASE_later %>% filter(p_adj < 0.05) %>% select(id, TIMEPOINT) %>% count(id) %>% arrange(desc(n))
  
  bru_results_dir = paste0('../../../HIER_MODEL_V2/REAL_DATA_RESULTS/pickle_files_FIX_MAJOR_ALLELE_GATK_VAR_merge_10_2/', CELL_LINE, '/')
  bru_timepoints = c('zero_hr', 'two_hr', 'six_hr')
  all_bru_list = list()
  for(timepoint in bru_timepoints){
    
    fn = paste0(bru_results_dir, timepoint, '/T10_M2_eiTrue_ARFalse_nsigma1.95996_p0.05_apfdr_bh_FIXED_MAJOR_ALLELE_ASE')
    in_dat = fread(fn)
    in_dat$TIMEPOINT = timepoint
    in_dat = in_dat %>% tidyr::separate_rows(SNVs, sep = ';') %>% select(Gene, ASEornonASE, SNVs, TIMEPOINT) %>% filter(SNVs != "")  %>% tidyr::separate(SNVs, into = c('coord', 'region', 'AR', 'ref_allelic_ratio', 'delta_AR', 'sigma', 'major', 'minor', 'p_orig', 'p_adj', 'nraw', 'nnorm'), sep = ',', convert = TRUE) %>% tidyr::separate(coord, into = c('chr', 'pos', 'strand'), sep = ':')
    in_dat = mutate(in_dat, id = paste0(chr, ':', pos, ':', strand))
    all_bru_list[[timepoint]] = in_dat
  }
  
  all_bru_df = do.call("rbind", all_bru_list)
  all_bru_df = all_bru_df %>% tidyr::separate(nnorm, into = c('ref1_norm', 'alt1_norm', 'ref2_norm', 'alt2_norm'), sep = '\\/')
  all_bru_df$ref1_norm = as.numeric(all_bru_df$ref1_norm)
  all_bru_df$alt1_norm = as.numeric(all_bru_df$alt1_norm)
  all_bru_df$ref2_norm = as.numeric(all_bru_df$ref2_norm)
  all_bru_df$alt2_norm = as.numeric(all_bru_df$alt2_norm)
  
  all_bru_df = all_bru_df %>% mutate(major_rep1 = pmax(ref1_norm, alt1_norm), major_rep2 = pmax(ref2_norm, alt2_norm), minor_rep1 = pmin(ref1_norm, alt1_norm), minor_rep2 = pmin(ref1_norm, alt1_norm))
  
  SPEC_SNV = SNV_ID
  
  gene_name = snv_level_dat %>% filter(id == SPEC_SNV) %>% select(Gene) %>% unique() %>% unlist(use.names = F)
  plot_dat = snv_level_dat %>% filter(id == SPEC_SNV) %>% as.data.frame() %>% tidyr::separate(nnorm, into = c('norm_ref1', 'norm_alt1', 'norm_ref2', 'norm_alt2', 'norm_ref3', 'norm_alt3'), sep = '\\/') %>% select(id, contains("norm"), TIMEPOINT) #%>% melt(id.vars = c('id', 'TIMEPOINT'))
  plot_dat$norm_ref1 = as.numeric(plot_dat$norm_ref1)
  plot_dat$norm_alt1 = as.numeric(plot_dat$norm_alt1)
  plot_dat$norm_ref2 = as.numeric(plot_dat$norm_ref2)
  plot_dat$norm_alt2 = as.numeric(plot_dat$norm_alt2)
  plot_dat$norm_ref3 = as.numeric(plot_dat$norm_ref3)
  plot_dat$norm_alt3 = as.numeric(plot_dat$norm_alt3)
  
  
  
  m_dat = melt(plot_dat, id.vars = c('id', 'TIMEPOINT'))
  m_dat$value = as.numeric(m_dat$value)
  m_dat$allele = gsub('[1-3]', '', m_dat$variable)
  m_dat$allele = gsub('norm_', '', m_dat$allele)
  m_dat$TIMEPOINT = factor(m_dat$TIMEPOINT, levels = c('0h', '2h', '8h', '24h'))
  
  
  stars_dat = snv_level_dat %>% filter(id == SPEC_SNV) %>% as.data.frame()  %>% select(TIMEPOINT, p_adj)
  stars_dat$sig = stars.pval(stars_dat$p_adj)
  
  max_pos = max(m_dat$value) -3
  #m_dat$allele = factor(m_dat$allele, levels = c('minor', 'major'))
  
 
  p = ggplot(m_dat, aes(x = TIMEPOINT,  y = value, fill = allele)) + geom_boxplot(aes(fill = allele))  + theme_classic() + ylab('Normalized counts') + xlab('Time post ActD treatment')  + geom_text(data = stars_dat, aes(x = TIMEPOINT, label = sig), y = max_pos, inherit.aes = FALSE, size = 7,  fontface = 'bold') + scale_fill_brewer(palette = 'BuPu') #+ ylim(c(0, max_pos))
  p1 = p + labs(title = paste0(SPEC_SNV, ' (', gene_name, ')')) + geom_point(position=position_dodge(width=0.75),aes(group=allele)) + geom_vline(xintercept = c(1.5,2.5,4.5), linetype = 'dashed', color = 'grey80')
  #print(final_p)
  
  # now get the bru dat
  spec_bru_dat = all_bru_df %>% filter(id == SPEC_SNV) %>% select(id, ref1_norm, alt1_norm, ref2_norm, alt2_norm, TIMEPOINT)
  bru_stars_dat = all_bru_df %>% filter(id == SPEC_SNV) %>% select(TIMEPOINT, p_adj)
  bru_stars_dat$sig = stars.pval(bru_stars_dat$p_adj)
  m_bru_dat = melt(spec_bru_dat, id.vars = c('id', 'TIMEPOINT'))
  m_bru_dat$allele = gsub("[1-2]_norm", "", m_bru_dat$variable)
  m_bru_dat$TIMEPOINT[m_bru_dat$TIMEPOINT == 'zero_hr'] = '0h'
  m_bru_dat$TIMEPOINT[m_bru_dat$TIMEPOINT == 'two_hr'] = '2h'
  m_bru_dat$TIMEPOINT[m_bru_dat$TIMEPOINT == 'six_hr'] = '6h'
  #m_bru_dat$allele = factor(m_bru_dat$allele, levels = c('minor', 'major'))
  bru_stars_dat$TIMEPOINT[bru_stars_dat$TIMEPOINT == 'zero_hr'] = '0h'
  bru_stars_dat$TIMEPOINT[bru_stars_dat$TIMEPOINT == 'two_hr'] = '2h'
  bru_stars_dat$TIMEPOINT[bru_stars_dat$TIMEPOINT == 'six_hr'] = '6h'
  p2 = ggplot(m_bru_dat, aes(x = TIMEPOINT, fill = allele, y = value)) + geom_boxplot() + ylab('Normalized counts') + xlab('Bru timepoint') + theme_classic() + scale_fill_brewer(palette = 'BuPu')
  bru_max_pos = max(m_bru_dat$value)
  p2 = p2 + geom_text(data = bru_stars_dat, aes(x = TIMEPOINT, label = sig), size = 7,  fontface = 'bold', inherit.aes = FALSE, y = bru_max_pos) + geom_point(position=position_dodge(width=0.75),aes(group=allele)) + geom_vline(xintercept = c(1.5,2.5,4.5), linetype = 'dashed', color = 'grey80')
  final_p = plot_grid(p1 + theme(legend.position = 'bottom'), p2 + theme(legend.position = 'bottom') + ggtitle(""), rel_widths = c(0.6, 0.4))
  print(final_p)
  #return(final_p)
  
  SNV_ID = gsub(':', '_', SNV_ID)
  
  
  out_fig_fn = paste0('../figures/', CELL_LINE, '_', gene_name, '_', SNV_ID, '_boxplot_validated_ASE_SNV_level_comparison.pdf')
  pdf(out_fig_fn, width = 8)
  final_p
  dev.off()
  m_dat = left_join(m_dat, stars_dat, by = 'TIMEPOINT')
  m_bru_dat = left_join(m_bru_dat, bru_stars_dat, by = 'TIMEPOINT') 
  out = list(m_dat, m_bru_dat, final_p) 
  return(out)
}



fig1 = make_fig('MCF-7', 'MCF7gDNA_S1_L006', 'chr2:17781251:+', 'GEN1')
fig2 = make_fig('GM12878', 'GM12878gDNA_S3_L006', 'chr17:81672520:+', 'CCDC137')
fig3 = make_fig('HCT116', 'HCT116gDNA_S2_L006', 'chr21:41888378:-', 'C2CD2')



final_fig = plot_grid(fig1[[3]], fig2[[3]], fig3[[3]], nrow = 1)


write.table(fig1[[1]], file = '../../../../SOURCE_DATA/MAIN_FIGS/fig4b_actD.txt', sep = '\t', quote = F, row.names = F)

write.table(fig1[[2]], file = '../../../../SOURCE_DATA/MAIN_FIGS/fig4b_Bru.txt', sep = '\t', quote = F, row.names = F)


write.table(fig2[[1]], file = '../../../../SOURCE_DATA/MAIN_FIGS/fig4c_actD.txt', sep = '\t', quote = F, row.names = F)

write.table(fig2[[2]], file = '../../../../SOURCE_DATA/MAIN_FIGS/fig4c_Bru.txt', sep = '\t', quote = F, row.names = F)



write.table(fig3[[1]], file = '../../../../SOURCE_DATA/MAIN_FIGS/fig4d_actD.txt', sep = '\t', quote = F, row.names = F)

write.table(fig3[[2]], file = '../../../../SOURCE_DATA/MAIN_FIGS/fig4d_Bru.txt', sep = '\t', quote = F, row.names = F)





pdf('../../../HIER_MODEL_V2/MANUSCRIPT_FIGURES/fig4/fig4b_to_d.pdf', width = 11, height = 4)
final_fig
dev.off()

pdf('../../../HIER_MODEL_V2/MANUSCRIPT_FIGURES/fig4/fig4b.pdf', width = 5, height = 2)
fig1[[3]] 
dev.off()

pdf('../../../HIER_MODEL_V2/MANUSCRIPT_FIGURES/fig4/fig4c.pdf', width = 11, height = 4)
fig2

dev.off()

pdf('../../../HIER_MODEL_V2/MANUSCRIPT_FIGURES/fig4/fig4d.pdf', width = 11, height = 4)
fig3

dev.off()



pdf('../../../HIER_MODEL_V2/MANUSCRIPT_FIGURES/fig4/fig4b_to_d_with_lines.pdf', width = 11, height = 4)
final_fig
dev.off()
