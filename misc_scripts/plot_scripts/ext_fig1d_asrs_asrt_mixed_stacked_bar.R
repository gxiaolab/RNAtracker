rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexUpset)
library(ComplexHeatmap)

#setwd('~/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')

all_snvs_fn = '../bed_files/all_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed'

all_snvs_dat = fread(all_snvs_fn)

names(all_snvs_dat) = c('chrom', 'pos0', 'pos1', 'gene_name', 'gene_state', 'gene_strand', 'type_snp', 'cell_line', 'genomic_region')


## get mean prev for mixed genes
combine_dat %>% mutate(prev = state_genes/total_genes) %>% filter(general_label == 'mixed') %>% select(prev) %>% ungroup() %>% summarize(mean = mean(prev))


## add mixed intronic asRT
all_intronic_res = fread('../bed_files/all_snvs_intronic_0h_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed')
names(all_intronic_res) = c('chrom', 'pos0', 'pos1', 'gene_name', 'gene_state', 'gene_strand', 'type_snp', 'cell_line', 'genomic_region')
all_intronic_res$gene_state[all_intronic_res$gene_state == 1] = 'intronic asRT 0h'
all_snvs_dat = rbind(all_snvs_dat, all_intronic_res)


# color pal
pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'intron-based asRT' = '#48d1cc40', 'mixed' ='#cbc0d3')


all_gene_label_dat = all_snvs_dat %>% select(gene_name, gene_state, cell_line) %>% unique()

## remove polyploid cell lines
all_gene_label_dat = all_gene_label_dat %>% filter(!cell_line %in% c('PC-9', 'PC-3', 'Caco-2'))

#pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'Mixed' = 'grey96') #'#8DA0CB')

all_gene_label_dat$general_label[all_gene_label_dat$gene_state %in% c('state1', 'state2')] = 'asRS'
all_gene_label_dat$general_label[all_gene_label_dat$gene_state == 'state3'] = 'asRT'
all_gene_label_dat$general_label[all_gene_label_dat$gene_state %in% c('state4', 'state5', 'state6')] = 'mixed'
all_gene_label_dat$general_label[all_gene_label_dat$gene_state %in% c('intronic asRT 0h')] = 'intron-based asRT'

# prev info
#total_counts = all_gene_label_dat %>% group_by(cell_line) %>% count() %>% rename(total_genes = n)
#state_counts = all_gene_label_dat %>% group_by(general_label, cell_line) %>% count() %>% rename(state_genes = n)
#combine_dat = inner_join(total_counts, state_counts, by = 'cell_line')


#combine_dat %>% filter(general_label %in% c('asRS', 'mixed')) %>% group_by(cell_line) %>% summarize(total_state_gnes = sum(state_genes), total_genes = total_genes) %>% mutate(prev = total_state_gnes/total_genes) %>% select(prev) %>% unique() %>% arrange(desc(prev))

#combine_dat %>% filter(general_label %in% c('asRS', 'mixed')) %>% group_by(cell_line) %>% summarize(total_state_gnes = sum(state_genes), total_genes = total_genes) %>% mutate(prev = total_state_gnes/total_genes) %>% select(prev) %>% unique() %>% arrange(desc(prev))


#combine_dat %>% filter(general_label %in% c('asRS', 'mixed')) %>% group_by(cell_line) %>% summarize(total_state_gnes = sum(state_genes), total_genes = total_genes) %>% mutate(prev = total_state_gnes/total_genes) %>% select(prev) %>% unique() %>% arrange(desc(prev))



plot_dat = filter(all_gene_label_dat, general_label %in% c('asRS', 'asRT', 'mixed', 'intron-based asRT'))

total_gene_count = plot_dat %>% group_by(cell_line) %>% count(cell_line)  %>% arrange(desc(n))
plot_dat$cell_line = factor(plot_dat$cell_line, levels = total_gene_count$cell_line)
p = ggplot(plot_dat, aes(fill = general_label, x = cell_line)) + 
  geom_bar(position="stack", stat="count", color = 'black', linewidth = 0.25) + scale_fill_manual(values = pal, name = 'Gene type', breaks = c('asRS', 'asRT', 'mixed', 'intron-based asRT')) + theme_bw() 


pdf("../figures/asrs_asrt_ambig_gene_counts.pdf", width = 10, height = 6.5)
p + ylab('Number of Genes') + xlab('Cell line') + theme_classic() +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dev.off()

p = ggplot(plot_dat %>% filter(general_label != 'intron-based asRT'), aes(fill = general_label, x = cell_line)) +
  geom_bar(position="stack", stat="count", color = 'black', linewidth = 0.25, width = 0.5) + scale_fill_manual(values = pal, name = 'Gene type', breaks = c('asRS', 'asRT', 'mixed')) + theme_classic()


write.table(plot_dat %>% filter(general_label != 'intron-based asRT'), file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig1b.txt', sep = '\t', row.names = F, quote = F)


pdf('../MANUSCRIPT_FIGURES/extended_fig1/extended_fig1d.pdf', width = 7, height = 2.5)
p + ylab('Number of Genes') + xlab('') + theme_classic() + theme(legend.position = 'top') #+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
#p + ylab('Number of Genes') + xlab('Cell line') + theme(
#  legend.box.spacing = unit(1, "pt")) + theme(axis.text=element_text(size=7),
#                                              axis.title=element_text(size=8)) +
#  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
#        legend.key.height = unit(0.5, 'cm'), #change legend key height
#        legend.key.width = unit(0.5, 'cm'), #change legend key width
#        legend.title = element_text(size=8), #change legend title font size
#        legend.text = element_text(size=7)) + #change legend text font size
#  theme(axis.text.x = element_text(angle = 30))
#dev.off()



pdf('../figures/asrs_asrt_ambig_gene_counts_no_intron_based.pdf', width = 7, height = 2.5)
p + ylab('Number of Genes') + xlab('Cell line') + theme(
  legend.box.spacing = unit(1, "pt")) + theme(axis.text=element_text(size=7),
                                              axis.title=element_text(size=8)) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=7)) + #change legend text font size
  theme(axis.text.x = element_text(angle = 30)  )
dev.off()



