rm(list=ls())
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)

#setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts/')


ADD_INTRONIC_BG = TRUE
out_fn = paste0("../results/ADD_INTRONIC_BG_", ADD_INTRONIC_BG, "_gtex_eQTL_gene_overlap.RData")
load(out_fn)

match_gtex = fread('../match_gtex_cell_line')
match_gtex = match_gtex %>% filter(!cell_line %in% c('K562', 'Panc1', 'PC-3', 'PC-9', 'Caco-2'))


match_gtex$gtex[match_gtex$gtex == 'Breast_Mammary_Tissue'] = 'Breast tissue'
match_gtex$gtex[match_gtex$gtex == 'Cells_EBV-transformed_lymphocytes'] = 'EBV-transformed lymphocytes'
match_gtex$gtex[match_gtex$gtex == 'Colon_Sigmoid'] = 'Sigmoid colon'
match_gtex$gtex[match_gtex$gtex == 'Colon_Transverse'] = 'Transverse colon'
match_gtex$gtex[match_gtex$gtex == 'Muscle_Skeletal'] = 'Skeletal Muscle'
match_gtex$gtex[match_gtex$gtex == 'Whole_Blood'] = 'Whole blood'


## convert Inf to max value
max_value = max(setdiff(all_results$OR, 'Inf'))
all_results$OR[all_results$OR == 'Inf'] = max_value


# code according to gtools
p_val_func = function(j, i, x, y, w, h, fill) {
  if (is.na(p_mat[i,j])){
    grid.text("", x,y)
  }
  else if(p_mat[i, j] <= 0.001 & OR_mat[i,j] > 1) {
    grid.text("***", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  } else if(p_mat[i, j] <= 0.01 & OR_mat[i,j] > 1) {
    grid.text("**", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  } else if (p_mat[i,j] <= 0.05 & OR_mat[i,j] > 1){
    grid.text("*", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  } else if (p_mat[i,j] <= 0.1 & OR_mat[i,j] > 1){
    grid.text(".", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  }  else {
    grid.text("", x,y)
  }
  
  if(paste0(i,j) %in% index_pairs){
    grid.rect(x,y, width = 1/10, height = 1/16, gp = gpar(lwd = 1, fill = "transparent"))
  }
}


make_heatmap = function(all_results) {
  #p_val_func = plot_func
  all_results = all_results %>% tidyr::separate(comparison, into = c('cell line', 'GTEx tissue'), sep = '\\/')
  all_results = all_results %>% filter(!`cell line` %in% c('K562', 'Panc1'))
  all_results = all_results %>% filter(`cell line` != 'all')
  all_results$`GTEx tissue`[all_results$`GTEx tissue` == 'Breast_Mammary_Tissue'] = 'Breast tissue'
  all_results$`GTEx tissue`[all_results$`GTEx tissue` == 'Cells_EBV-transformed_lymphocytes'] = 'EBV-transformed lymphocytes'
  all_results$`GTEx tissue`[all_results$`GTEx tissue` == 'Colon_Sigmoid'] = 'Sigmoid colon'
  all_results$`GTEx tissue`[all_results$`GTEx tissue` == 'Colon_Transverse'] = 'Transverse colon'
  all_results$`GTEx tissue`[all_results$`GTEx tissue` == 'Muscle_Skeletal'] = 'Skeletal Muscle'
  all_results$`GTEx tissue`[all_results$`GTEx tissue` == 'Whole_Blood'] = 'Whole blood'
  all_results$`cell line`[all_results$`cell line` == 'mec'] <- 'HMEC'
  all_results$`cell line`[all_results$`cell line` == 'ecuv'] <- 'HUVEC'
  all_results$`cell line`[all_results$`cell line` == 'PC9'] <- 'PC-9'
  all_results$`cell line`[all_results$`cell line` == 'MCF_7'] <- 'MCF-7'
  all_results$`cell line`[all_results$`cell line` == 'MCF10A'] <- 'MCF 10A'
  all_results$`cell line`[all_results$`cell line` == 'OCILY7'] <- 'OCI-LY7'
  all_results$`cell line`[all_results$`cell line` == 'PC_3'] <- 'PC-3'
  all_results = all_results %>% arrange(`cell line`)
  OR_dat = all_results %>% select(`cell line`, `GTEx tissue`, OR)
  p_dat = all_results %>% select(`cell line`, `GTEx tissue`, p_value)
  
  OR_dat2 = dcast(data = OR_dat, formula = `cell line`~`GTEx tissue`, value.var = "OR")
  p_dat2 = dcast(data = p_dat, formula = `cell line` ~ `GTEx tissue`, value.var = 'p_value')
  
  OR_mat = as.matrix(OR_dat2[,-1])
  row.names(OR_mat) = OR_dat2$`cell line`
  p_mat = as.matrix(p_dat2[,-1])
  row.names(p_mat) = p_dat2$`cell line`
  #print(p_mat)
  #max_OR = round(max(OR_mat))
  #min_OR = round(min(OR_mat))
  #col_fun = colorRamp2(c(1,0.5,0), magma(3))
  #out = Heatmap(OR_mat, name = 'Enrichment', cell_fun = p_val_func, col = rev(plasma(6)))
  #return(out)
  
  ## get the indexes of matching cell line/tissue pairs
  index_list = rep("", nrow(match_gtex))
  for (i in 1:nrow(match_gtex)){
    cell_line = match_gtex$cell_line[i]
    tissue = match_gtex$gtex[i]
    x = which(colnames(OR_mat) == eval(quote(tissue)))
    y = which(row.names(OR_mat) == eval(quote(cell_line)))
    index_list[i] = paste0(y,x)
  }
  out = list(OR_mat, p_mat, index_list)
}


all_results = all_results %>% filter(!comparison %like% 'PC-9') %>% filter(!comparison %like% 'PC-3') %>% filter(!comparison %like% 'Caco-2')

plot_dat = make_heatmap(all_results)

col_fun = colorRamp2(c(0, 1.5, 3, 4.5, 6), brewer.pal(5, "YlOrRd"))
#col_fun = colorRamp2(c(0, 0.5, 1, 1.5, 2), rev(brewer.pal(5, "RdYlBu")))
#col_fun = colorRamp2(c(0,1,3), plasma(3))
#col_fun2 = colorRamp2(c(0,1.5,3), plasma(3))


## get the index pairs 
p_mat = plot_dat[[2]]
index_pairs = plot_dat[[3]]
col_fun = colorRamp2(c(0, 1, 2, 3, 4, 5,6), rev(brewer.pal(11, "RdYlBu")[2:8]))
OR_mat = plot_dat[[1]]
h1 = Heatmap(OR_mat, name = 'asRS\nEnrichment', cell_fun = p_val_func, col = col_fun, cluster_rows = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,  column_names_gp = grid::gpar(fontsize = 7),
  row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(labels_gp = gpar(fontsize = 8)))
h1



write.table(OR_mat, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig7b_OR.txt', sep = '\t', row.names = F, quote = F)

write.table(p_mat, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig7b_pval.txt', sep = '\t', row.names = F, quote = F)


pdf('../figures/asrs_asrt_gtex_eqtl_genes_overlap_heatmap.pdf')
h1
dev.off()


pdf('../MANUSCRIPT_FIGURES/extended_fig7/extended_fig7c.pdf', width = 3.5, height = 4.5)
h1
dev.off()

