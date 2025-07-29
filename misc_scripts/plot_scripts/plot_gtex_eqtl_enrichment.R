library(data.table)
library(dplyr)
library(gtools)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)

#setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/gene_strand_fix/scripts')
#setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts/')

match_gtex = fread('../match_gtex_cell_line')
match_gtex = filter(match_gtex, !cell_line %in% c('PC-3', 'PC-9', 'Caco-2', 'Panc1', 'K562'))

ADD_INTRONIC_BG = "TRUE"


match_gtex$gtex[match_gtex$gtex == 'Breast_Mammary_Tissue'] = 'Breast tissue'
match_gtex$gtex[match_gtex$gtex == 'Cells_EBV-transformed_lymphocytes'] = 'EBV-transformed lymphocytes'
match_gtex$gtex[match_gtex$gtex == 'Colon_Sigmoid'] = 'Sigmoid colon'
match_gtex$gtex[match_gtex$gtex == 'Colon_Transverse'] = 'Transverse colon'
match_gtex$gtex[match_gtex$gtex == 'Muscle_Skeletal'] = 'Skeletal Muscle'
match_gtex$gtex[match_gtex$gtex == 'Whole_Blood'] = 'Whole blood'

match_gtex$cell_line[match_gtex$cell_line == 'mec'] <- "HMEC"
match_gtex$cell_line[match_gtex$cell_line == 'ecuv'] <- "HUVEC"
match_gtex$cell_line[match_gtex$cell_line == 'PC9'] <- 'PC-9'
match_gtex$cell_line[match_gtex$cell_line == 'MCF_7'] <- 'MCF-7'
match_gtex$cell_line[match_gtex$cell_line == 'MCF10A'] <- 'MCF10A'
match_gtex$cell_line[match_gtex$cell_line == 'OCILY7'] <- 'OCI-LY7'
match_gtex$cell_line[match_gtex$cell_line == 'PC_3'] <- 'PC-3'

load(paste0("../results/MATCH_GENIC_REGION_REMOVE_5_CELL_LINES_ADD_INTRONIC_BG_", ADD_INTRONIC_BG, "_gtex_eQTL_asrs_enrichment.RData"))
load(paste0("../results/MATCH_GENIC_REGION_REMOVE_5_CELL_LINES_ADD_INTRONIC_BG_", ADD_INTRONIC_BG, "_gtex_eQTL_asrt_enrichment.RData"))

max_asrt_OR = max(setdiff(asrt_results$OR, 'Inf'))
max_asrs_OR = max(setdiff(asrs_results$OR, 'Inf'))

#?stars.pval
match_gtex = match_gtex %>% filter(!cell_line %in% c('K562', 'Panc1', 'PC-3', 'PC-9', 'Caco-2'))

make_heatmap = function(all_results) {
	#p_val_func = plot_func
	all_results = all_results %>% tidyr::separate(comparison, into = c('cell line', 'GTEx tissue'), sep = '\\/')
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


asrs_dat = make_heatmap(asrs_results)
asrt_dat = make_heatmap(asrt_results)

#col_fun = colorRamp2(c(0,3,6), plasma(3))
#col_fun2 = colorRamp2(c(0,1.5,3), plasma(3))


col_fun = colorRamp2(c(0, 1, 2, 3, 4, 5,6), rev(brewer.pal(11, "RdYlBu")[2:8]))
#col_fun = colorRamp2(c(0, 1.5, 3, 4.5, 6), brewer.pal(5, "YlOrRd"))
#col_fun2 = colorRamp2(c(0, 0.75, 1.5, 2.25, 3), rev(brewer.pal(5, "RdYlBu")))

## get the index pairs and p-value matrices
asrs_p_mat = asrs_dat[[2]]
asrs_index_pairs = asrs_dat[[3]]
asrt_p_mat = asrt_dat[[2]]
asrt_index_pairs = asrt_dat[[3]]


# define cell functions
# code according to gtools
asrs_p_val_func = function(j, i, x, y, w, h, fill) {
  if (is.na(asrs_p_mat[i,j])){
    grid.text("", x,y)
  }
  else if(asrs_p_mat[i, j] <= 0.001 & asrs_dat[[1]][i, j] > 1) {
    grid.text("***", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  } else if(asrs_p_mat[i, j] <= 0.01 & asrs_dat[[1]][i, j] > 1) {
    grid.text("**", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  } else if (asrs_p_mat[i,j] <= 0.05 & asrs_dat[[1]][i, j] > 1){
    grid.text("*", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  # remove symbols for p<=0.1
  ##} else if (asrs_p_mat[i,j] <= 0.1){
   # grid.text(".", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  }  else {
    grid.text("", x,y)
  }
  
  if(paste0(i,j) %in% asrs_index_pairs){
    grid.rect(x,y, width = 1/10, height = 1/16, gp = gpar(lwd = 0.5, fill = "transparent"))
  }
}


asrt_p_val_func = function(j, i, x, y, w, h, fill) {
  if (is.na(asrt_p_mat[i,j])){
    grid.text("", x,y)
  }
  else if(asrt_p_mat[i, j] <= 0.001 & asrt_dat[[1]][i, j] > 1 ) {
    grid.text("***", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  } else if(asrt_p_mat[i, j] <= 0.01 & asrt_dat[[1]][i, j] > 1 ) {
    grid.text("**", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  } else if (asrt_p_mat[i,j] <= 0.05 & asrt_dat[[1]][i, j] > 1 ){
    grid.text("*", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  #} else if (asrt_p_mat[i,j] <= 0.1){
  #  grid.text(".", x, y, gp = gpar(fontface = "bold", fontsize = 8))
  }  else {
    grid.text("", x,y)
  }
  
  if(paste0(i,j) %in% asrt_index_pairs){
    grid.rect(x,y, width = 1/10, height = 1/16, gp = gpar(lwd = 0.5, fill = "transparent"))
  }
}
asrs_dat[[1]][asrs_dat[[1]] == Inf] = max_asrs_OR
asrt_dat[[1]][asrt_dat[[1]] == Inf] = max_asrt_OR

h1 = Heatmap(asrs_dat[[1]], name = 'asRS Enrichment', cell_fun = asrs_p_val_func, col = col_fun, show_row_dend = FALSE, show_column_dend = FALSE)
h1
Heatmap(asrs_dat[[1]], name = 'asRS Enrichment', cell_fun = asrs_p_val_func, show_row_dend = FALSE, show_column_dend = FALSE)



h2 = Heatmap(asrt_dat[[1]], name = 'asRT Enrichment', cell_fun = asrt_p_val_func, col = col_fun, show_row_dend = FALSE, show_column_dend = FALSE )
h2

#Heatmap(asrt_dat[[1]], name = 'asRT Enrichment', cell_fun = asrt_p_val_func, show_row_dend = FALSE, show_column_dend = FALSE )



legend_par = list(labels_gp = gpar(fontsize = 6), title_gp = gpar(fontsize = 8))


h3 = Heatmap(asrs_dat[[1]], name = 'asRS Enrichment', cell_fun = asrs_p_val_func, col = col_fun,  show_row_dend = FALSE, show_column_dend = FALSE, column_names_gp = grid::gpar(fontsize = 7),
  row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = legend_par)
h3

h4 = Heatmap(asrt_dat[[1]], name = 'asRT Enrichment', cell_fun = asrt_p_val_func, col = col_fun, show_row_dend = FALSE, show_column_dend = FALSE, column_names_gp = grid::gpar(fontsize = 7),
  row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = legend_par)
#pdf('../figures/eQTL_asrs_asrt_enrichment.pdf', width = 11)


out_fn = paste0('../figures/match_by_genic_region_ADD_INTRONIC_BG_', ADD_INTRONIC_BG, '_eQTL_asrs_asrt_enrichment.pdf')



write.table(asrs_dat[[1]], file = '../../../SOURCE_DATA/MAIN_FIGS/fig2a_asrs_enrichment.txt', sep = '\t', row.names = F, quote = F)

write.table(asrt_dat[[1]], file = '../../../SOURCE_DATA/MAIN_FIGS/fig2a_asrt_enrichment.txt', sep = '\t', row.names = F, quote = F)


write.table(asrs_dat[[2]], file = '../../../SOURCE_DATA/MAIN_FIGS/fig2a_asrs_pval.txt', sep = '\t', row.names = F, quote = F)
write.table(asrt_dat[[2]], file = '../../../SOURCE_DATA/MAIN_FIGS/fig2a_asrt_pval.txt', sep = '\t', row.names = F, quote = F)



pdf(out_fn, width = 11)
h3 + h4
dev.off()


pdf('../MANUSCRIPT_FIGURES/fig2/fig2a.pdf', width = 5.5, height = 3)
h3 + h4
dev.off()


# remove the all cell line comparison
asrs_results = asrs_results %>% filter(!comparison %like% "all")
asrt_results = asrt_results %>% filter(!comparison %like% "all")

## now see prop of asrs that overlap eQTLs
asrs_results$total_asrs = asrs_results$number_hits + asrs_results$number_non_hits
asrs_results$asrs_prop = asrs_results$number_hits/asrs_results$total_asrs

## now see prop of asrt that overlap eQTLs
asrt_results$total_asrt = asrt_results$number_hits + asrt_results$number_non_hits
asrt_results$asrt_prop = asrt_results$number_hits/asrt_results$total_asrt

prop_dat = data.frame(asRS = asrs_results$asrs_prop, asRT = asrt_results$asrt_prop)
m_prop_dat = melt(prop_dat)
names(m_prop_dat) = c('SNV_type', 'prop')

library(ggpubr)
p0 = ggplot(m_prop_dat, aes(x = SNV_type, y = prop)) + 
  geom_violin(aes(fill = SNV_type), width = 0.7) + 
  geom_boxplot(width=0.2, color="black", alpha = 0.5, aes(fill = SNV_type)) + 
  scale_fill_manual(values = c("#875f9a","#48d1cc", "#a5a5a5")) + stat_compare_means(size = 4, hjust = -0.5) + theme_classic() + ylab('Proportion of variants overlapping eQTLs') + xlab('Variant type')


out_fn2 = paste0('../figures/match_by_genic_region_ADD_INTRONIC_BG_', ADD_INTRONIC_BG, '_prop.overlap.gtex.eQTL.pdf')


pdf(out_fn2)
p0
dev.off()

prop_dat = data.frame(comparison = asrs_results$comparison, asRS = asrs_results$asrs_prop, asRT = asrt_results$asrt_prop)
#m_prop_dat = melt(prop_dat)
#names(m_prop_dat) = c('SNV_type', 'prop')


write.table(prop_dat, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig7a.txt', sep = '\t', row.names = F, quote = F)




pdf('../MANUSCRIPT_FIGURES/extended_fig7/extended_fig7b.pdf', width = 4, height = 4)
p0 + theme(legend.position = 'none')
dev.off()

