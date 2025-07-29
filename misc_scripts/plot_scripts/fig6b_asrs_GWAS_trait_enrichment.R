rm(list = ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(R.utils)
library(unikn)
library(ggrepel)
library(gtools)
library(cowplot)
library(gwasrapidd)


#setwd('~/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/GWAS_CATALOG/scripts')
#load('asrs_snps_gwas_pval_dat.RData')

ADJ_P = 0.05
GENE_CUTOFF = 3
#SNV_CUTOFF = 5

#load('../results/asrs_snps_gwas_add_direct_overlaps_adj_p_0.1_pval_dat.RData')
#load('../results/asrt_snps_gwas_add_direct_overlaps_adj_p_0.1_pval_dat.RData')

# get list of immune system traits
EFO_id = 'EFO_0000540' # immune system disease EFO ID

#get_immune_traits = get_child_efo(EFO_id)
#save(get_immune_traits, file = '../IMD_ref_files/immune_EFO.RData')

load('../../../../gene_strand_fix/GWAS_CATALOG/gwas_catalog_IMD_fix/IMD_ref_files/immune_EFO.RData')

#enrichment_results_fn = "../results/empirical.test.enrichment.results10000perms.add.direct.and.total.counts.ANNOTATED"
enrichment_results_fn = '../results/empirical.test.enrichment.results10000perms.add.direct.and.total.counts'
enrichment_results = fread(enrichment_results_fn)
enrichment_results = enrichment_results %>% filter(trait != 'Acute myeloid leukemia')
enrichment_results$is_immune = enrichment_results$efo_id %in% get_immune_traits[[1]]
write.table(enrichment_results, file = "../results/empirical.test.enrichment.results10000perms.add.direct.and.total.counts.ANNOTATED", quote = F, sep = '\t', row.names = F)

#enrichment_results$trait = capitalize(enrichment_results$trait)
#immune_traits = enrichment_results %>% filter(is_immune)

## need to take out acute myeloid leukemia because it came from that sketchy study 


#top_results = traits_of_interest %>% rename(OR = asrs_OR_vec)
#top_results$sig = stars.pval(top_results$asrs_adj.p)
#top_results = mutate(top_results, Trait = paste0(Trait, '\n', asrs_asrs_count))
#top_results =  top_results %>% mutate(Trait = fct_reorder(Trait, desc(OR)))
#top_results$asrt_compare = top_results$delta_RR > 0
#title1 = paste0('asRS results, adj.p < ', adj_p, ' gene count >= ', gene_count)

load('../results/asrs_snps_gwas_add_direct_overlaps_ALL_DISEASE.RData')

asrs_snps_gwas_pval_dat$trait = capitalize(asrs_snps_gwas_pval_dat$trait)
asrs_trait_counts = asrs_snps_gwas_pval_dat %>% select(trait, gene_name) %>% unique() %>% count(trait) %>% rename(asrs_gene_count = n)

#load('../IMD_ref_files/immune_EFO.RData')
load('../../../../gene_strand_fix/GWAS_CATALOG/gwas_catalog_IMD_fix/IMD_ref_files/immune_EFO.RData')
enrichment_results_fn = "../results/empirical.test.enrichment.results10000perms.add.direct.and.total.counts.ANNOTATED"
enrichment_results = fread(enrichment_results_fn)
enrichment_results = enrichment_results %>% filter(trait != 'acute myeloid leukemia')
enrichment_results$trait = capitalize(enrichment_results$trait)
enrichment_results = inner_join(enrichment_results, asrs_trait_counts, by = 'trait')

immune_traits = enrichment_results %>% filter(is_immune)

traits_of_interest = enrichment_results %>% filter(asrs_FDR < ADJ_P) %>% filter(is_disease) %>% filter(asrs_gene_count >= GENE_CUTOFF)# %>% filter(is_immune)
traits_of_interest = traits_of_interest %>% filter(trait != 'acute myeloid leukemia') 
sig_traits = traits_of_interest$trait
traits_of_interest$log10p = -log10(traits_of_interest$asrs_FDR)


boxplot_dat = asrs_snps_gwas_pval_dat %>% filter(trait %in% traits_of_interest$trait)
boxplot_dat$is_immune = boxplot_dat$trait %in% immune_traits$trait  
boxplot_dat = boxplot_dat %>% select(gwas_SNP, neg_log10p, trait) %>% unique()


## get trait order based on mean GWAS p-value


median_disease_logp = asrs_snps_gwas_pval_dat %>% select(gwas_SNP, neg_log10p, trait) %>% unique() %>% group_by(trait) %>% summarize(median_log = median(neg_log10p)) %>% arrange(desc(median_log)) 
#disease_trait_mean_logp =  mean_disease_logp %>% summarize(mean_mean_log = mean(mean_log))
#sig_traits_mean_logp = filter(mean_disease_logp, trait %in% trait_list$trait)


trait_order = median_disease_logp %>%  filter(trait %in% traits_of_interest$trait) %>% arrange(median_log)
boxplot_dat$trait = factor(boxplot_dat$trait, levels = trait_order$trait)
traits_of_interest$trait = factor(traits_of_interest$trait, levels = trait_order$trait)

# get colors in the right order
n <- length(sig_traits[sig_traits %in% immune_traits$trait])
col_vector = brewer.pal(n, name = 'Set3')
color_lookup = data.frame(trait = c(sig_traits[sig_traits %in% immune_traits$trait], sig_traits[!sig_traits %in% immune_traits$trait]), color = c(col_vector, rep('grey90', length(sig_traits[!sig_traits %in% immune_traits$trait]))))

## replace autoimmune thyroid gray color
color_lookup$color[color_lookup$trait == 'Autoimmune thyroid disease'] = brewer.pal(11, name = 'Set3')[11]

color_lookup$trait = factor(color_lookup$trait, levels = trait_order$trait)
color_lookup = color_lookup %>% arrange(trait)

save(color_lookup, file = 'trait.color.pal.RData')


p1 = ggplot(boxplot_dat, aes(x = trait, y = neg_log10p)) + 
  geom_boxplot(aes(fill = trait)) + coord_flip() + scale_fill_manual(values = color_lookup$color) + theme_classic() + ylab(expression("-log"["10"]~"(GWAS"["p"]~")")) + xlab('Trait') + 
  theme(legend.position = 'none')  


p2 = ggplot(data = traits_of_interest %>% rename(`# of asRS genes` = asrs_gene_count), aes(x = trait, y = log10p, group=1)) +
  geom_line()+
  geom_point(aes(size = `# of asRS genes`, fill = trait), color = 'black', pch = 21) + theme_classic() + coord_flip() + ylab(expression("-log"["10"]~"(enrich"["p"]~")")) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + scale_fill_manual(values = color_lookup$color) + 
  #scale_size_continuous(range = c(3,10))+ 
  guides(fill = "none") 

# need to break apart p1
#p_tmp1 = p1 + scale_y_continuous(limits = c(0, 130)) 
#p_tmp2 = p1 + scale_y_continuous(limits = c(290, 310), breaks = c(290, 300, 310)) + xlab("") +  theme(axis.ticks = element_blank(), axis.text.y = element_blank()) 

#p1_complete = plot_grid(p_tmp1, p_tmp2, ncol = 2, rel_widths = c(1,20/130))

#final_out = plot_grid(p1_complete, p2 , rel_widths = c(2,0.8))
final_out = plot_grid(p1, p2 , rel_widths = c(2,0.8))


## maybe later will change color scale so each immune disease has its own color

pdf('../figures/asrs_empirical_test_enrichment_boxplot_and_lineplot_combined.pdf', height = 11, width = 13)
final_out
dev.off()
          
pdf('../../MANUSCRIPT_FIGURES/fig6/fig6b.pdf', width = 8.5, height = 5)
final_out  + theme(
            legend.box.spacing = unit(1, "pt")) + theme(axis.text=element_text(size=6),
           axis.title=element_text(size=8)) + 
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) #change legend text font size
dev.off()

