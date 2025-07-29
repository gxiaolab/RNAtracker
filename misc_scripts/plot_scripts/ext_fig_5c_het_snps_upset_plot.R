library(data.table)
library(dplyr)
library(VGAM)
library(ComplexUpset)
library(ComplexHeatmap)
library(ggplot2)
source("~/project-gxxiao4/softwares/loadRData.R")



COVERAGE_CUTOFF = 10
MISMATCH_CUTOFF = 2
out_suffix = paste0('GATK_VAR_merge_', COVERAGE_CUTOFF, '_', MISMATCH_CUTOFF)
P_ADJ_CUTOFF = 0.05
SNP_COUNT_FILTER = 2

rdd_input_dir = '/u/home/h/huange7/project-gxxiao4/ASE_ASRS/REVISION1/WASP/GATK_VCF/RDD_format'

chromosomes = paste0('chr', 1:22)

read_rdd = function(CELL_LINE){
    ALL_SNVS = list()
    for (chrom in chromosomes){
        fn = paste0(rdd_input_dir, '/', CELL_LINE, '.MERGED.PASS_ONLY.het_snps.', chrom, '.rdd')
        chrom_dat = fread(fn)
        chrom_dat = mutate(chrom_dat, id = paste0(V1, ':', V2 + 1, ':', V3, ':', V4))
        ALL_SNVS[[chrom]] = chrom_dat
    }
    all_chrom_dat = do.call("rbind", ALL_SNVS)
    return(all_chrom_dat$id)
}

cell_lines = fread('../cell_lines_14_final', header = F)
cell_lines = cell_lines$V1
all_snv_list = list()

#for(cell_line in cell_lines){
#    all_snv_list[[cell_line]] = read_rdd(cell_line)
#}

#names(all_snv_list) = cell_lines

#snv_mat = list_to_matrix(all_snv_list)
#snv_tbl = as_tibble(snv_mat)

#save(snv_mat, snv_tbl, file = '../results/make_all_het_snps_upset.RData')

load('../results/make_all_het_snps_upset.RData')


number_upset_intersections = 25
PLOT_PREFIX_TITLE = 'all het'
PLOT_COLOR = 'darkblue'
snv_tbl = snv_tbl %>% select(-c('Caco-2', 'PC-3', 'PC-9'))
p2 = upset(snv_tbl, names(snv_tbl), n_intersections = number_upset_intersections, name = paste0(PLOT_PREFIX_TITLE, " SNVs"), 
             base_annotations = list('Intersection size' = intersection_size(counts = FALSE, mapping = aes(fill = 'bars_color'), width = 0.6) +
                                       scale_fill_manual(values = c('bars_color' = PLOT_COLOR), guide = 'none')),
             width_ratio = 0.2, height_ratio = 1, themes=upset_default_themes(text=element_text(size=7)), matrix=(intersection_matrix(geom=geom_point(size=1.2))),   set_sizes=(upset_set_size()
                                                                                                                                                                                + theme(axis.text.x=element_text(angle=90))))




#capture.output(all_snv_list, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig5c.txt')

pdf('../figures/upset_all_het_snps_test.pdf',width = 4.5, height = 3)
p2
dev.off()        

pdf('../MANUSCRIPT_FIGURES/extended_fig5/extended_fig5c.pdf', width = 3.5, height = 2.5)
p2
dev.off()

                                                                                                                                                                                
