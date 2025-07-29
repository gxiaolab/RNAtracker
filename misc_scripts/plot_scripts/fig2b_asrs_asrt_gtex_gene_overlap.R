rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(gtools)
library(sjmisc)


#setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts/')

ADD_INTRONIC_BG = "TRUE"
load(paste0("../results/ADD_INTRONIC_BG_", ADD_INTRONIC_BG, "_gtex_eQTL_gene_overlap.RData"))

# just plot the combined cell lines results
all_results = all_results %>% filter(comparison %like% "all")
all_results$comparison = gsub('all/', '', all_results$comparison)
all_results$comparison[all_results$comparison == 'Breast_Mammary_Tissue'] = 'Breast tissue'
all_results$comparison[all_results$comparison == 'Cells_EBV-transformed_lymphocytes'] = 'EBV-transformed lymphocytes'
all_results$comparison[all_results$comparison == 'Colon_Sigmoid'] = 'Sigmoid colon'
all_results$comparison[all_results$comparison == 'Colon_Transverse'] = 'Transverse colon'
all_results$comparison[all_results$comparison == 'Muscle_Skeletal'] = 'Skeletal Muscle'
all_results$comparison[all_results$comparison == 'Whole_Blood'] = 'Whole blood'


all_results$comparison = word_wrap(all_results$comparison, 12)


# % eQTL that overlaps asRS or asRT
plot_dat1 = all_results %>% dplyr::select(comparison, asrs_eqtl_prop, asrt_eqtl_prop) %>% melt(id.vars = 'comparison')
plot_dat1$variable = as.character(plot_dat1$variable)
plot_dat1$variable[plot_dat1$variable == 'asrs_eqtl_prop'] = 'asRS'
plot_dat1$variable[plot_dat1$variable == 'asrt_eqtl_prop'] = 'asRT'
names(plot_dat1) = c('tissue', 'gene type', 'prop')

# % asRS or asRT that overlaps eQTL
plot_dat2 = all_results %>% select(comparison, asrs_prop, asrt_prop) %>% melt(id.vars = 'comparison')
plot_dat2$variable = as.character(plot_dat2$variable)
plot_dat2$variable[plot_dat2$variable == 'asrs_prop'] = 'asRS'
plot_dat2$variable[plot_dat2$variable == 'asrt_prop'] = 'asRT'
names(plot_dat2) = c('tissue', 'gene type', 'prop')

pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc')
plot_dat1$perc = plot_dat1$prop * 100
p1 = ggplot(plot_dat1, aes(fill = `gene type`, y = perc, x = tissue)) + 
  geom_bar( stat="identity", alpha = 0.67, color = 'black', width = 0.5, position = position_dodge(width = 0.6)) + scale_fill_manual(values = pal) + ylab("% eQTL genes")


plot_dat2$perc = plot_dat2$prop * 100 
plotting_order = plot_dat2 %>% filter(`gene type` == 'asRS') %>% arrange(prop)
plot_dat2$tissue = factor(plot_dat2$tissue, levels = plotting_order$tissue)
p2_tmp = ggplot(plot_dat2, aes(fill = `gene type`, y = perc, x = tissue)) + 
  geom_bar(width = 0.5, position = position_dodge(width = 0.6), stat="identity", color = 'black', size = 0.25) + scale_fill_manual(values = pal) + ylab("% of genes overlapping eGenes")

write.table(plot_dat2 %>% mutate(tissue = gsub("\n", "", plot_dat2$tissue)), file = '../../../SOURCE_DATA/MAIN_FIGS/fig2b.txt', sep = '\t', row.names = F, quote = F)

stars_dat_tmp = plot_dat2 %>% group_by(tissue) %>% slice_max(perc) # this is just for getting the position to plot stars
stars_dat = inner_join(stars_dat_tmp, all_results %>% select(comparison, p_value) %>% rename(tissue = comparison), by = 'tissue')
stars_dat$stars = stars.pval(stars_dat$p_value)

p2 = p2_tmp + geom_text(data = stars_dat,aes(label = stars), size = 5, nudge_y = 5, angle = 90) + xlab('Tissue') + coord_flip() #+ theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

out_fn = paste0("../figures/ADD_INTRONIC_BG_", ADD_INTRONIC_BG, "_asrs_asrt_gtex_eqtl_genes_overlap.pdf")


#write.table(stars_dat %>% as.data.frame(), file = '../../../SOURCE_DATA/MAIN_FIGS/fig2b.txt', sep = '\t', row.names = F, quote = F)

out = left_join(plot_dat2, stars_dat %>% as.data.frame() %>% select(-c(prop, perc)), by = 'tissue') %>% select(-perc)  %>% mutate(tissue = gsub("\n", "", tissue))
write.table(out, file = '../../../SOURCE_DATA/MAIN_FIGS/fig2b.txt', sep = '\t', row.names = F, quote = F)


pdf(out_fn, width = 10, height = 6.5)
p1
p2 + theme_bw()
dev.off()



pdf('../MANUSCRIPT_FIGURES/fig2/fig2b.pdf', width = 3, height = 4)
p2 + theme_classic() + theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8)) + guides(fill=guide_legend(title="Gene Type"))  +
  theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6), legend.position = 'bottom') #change legend text font size
dev.off()

