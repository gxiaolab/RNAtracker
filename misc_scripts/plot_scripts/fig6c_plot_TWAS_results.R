library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')

TWAS_res = fread('../results/post_reiview_TWAS_all.tsv')
asrs_list = fread('../snv_types/ASRS_snvs_0.95_1_0.333333333333333_0.95_1_0.333333333333333_PASS_PROP_THRESHOLD_0.5.bed')

# before filtering out the diseases i'm not interested in anymore
TWAS_res = TWAS_res %>% mutate(label = paste0(gene_name, '_', disease, '_', tissue))

genic_res = filter(TWAS_res, type == 'genic') #%>% select(TWAS.P, label ) %>% filter(!is.na(TWAS.P)) %>% rename(genic_TWAS.P = TWAS.P) 
cis_res = filter(TWAS_res, type == 'cis') #%>% select(TWAS.P, label )  %>% filter(!is.na(TWAS.P)) %>% rename(cis_TWAS.P = TWAS.P)


sig_genic_res = filter(TWAS_res, type == 'genic') %>% filter(p.adj < 0.05)
sig_cis_res = filter(TWAS_res, type == 'cis') %>% filter(p.adj < 0.05) # #%>% select(TWAS.P, label )  %>% filter(!is.na(TWAS.P)) %>% rename(cis_TWAS.P = TWAS.P)


intersect(sig_genic_res$label, sig_cis_res$label) %>% length() # 23, which matches ryo's venn diagram. so he is plotting gene/disease/tissue pairs
setdiff(sig_genic_res$description, sig_cis_res$description) %>% length() # 19
setdiff(sig_cis_res$description, sig_genic_res$description) %>% length() # 9

## after re-running psoriasis with skin, the numbers above are 21, 17, 10 




# get rid of psoriasis and T1D, also ankylosing spondylitis is no longer a trait of interest
TWAS_res = TWAS_res %>% filter(!disease %in% c('psoriasis', 'T1D', 'ankylosing_spondylitis'))
genic_res = filter(TWAS_res, type == 'genic') #%>% select(TWAS.P, label ) %>% filter(!is.na(TWAS.P)) %>% rename(genic_TWAS.P = TWAS.P) 
cis_res = filter(TWAS_res, type == 'cis') #%>% select(TWAS.P, label )  %>% filter(!is.na(TWAS.P)) %>% rename(cis_TWAS.P = TWAS.P)


sig_genic_res = filter(TWAS_res, type == 'genic') %>% filter(p.adj < 0.05)
sig_cis_res = filter(TWAS_res, type == 'cis') %>% filter(p.adj < 0.05) # #%>% select(TWAS.P, label )  %>% filter(!is.na(TWAS.P)) %>% rename(cis_TWAS.P = TWAS.P)


intersect(sig_genic_res$label, sig_cis_res$label) %>% length() # 23
setdiff(sig_genic_res$description, sig_cis_res$description) %>% length() # 15
setdiff(sig_cis_res$description, sig_genic_res$description) %>% length() # 9

## after re-running psoriasis with skin, the numbers above are 21, 13, 10


setdiff(sig_genic_res$label, sig_cis_res$label)


genic_res = filter(TWAS_res, type == 'genic') %>% select(p.adj, label ) %>% filter(!is.na(p.adj)) %>% rename(genic_TWAS.P = p.adj) 
cis_res = filter(TWAS_res, type == 'cis') %>% select(p.adj, label )  %>% filter(!is.na(p.adj)) %>% rename(cis_TWAS.P = p.adj)


inner_join(genic_res, cis_res, by = 'label') %>% mutate(genic_more_sig = genic_TWAS.P < cis_TWAS.P) %>% filter(genic_TWAS.P < 0.05) %>% count(genic_more_sig)
inner_join(genic_res, cis_res, by = 'label') %>% mutate(genic_more_sig = genic_TWAS.P < cis_TWAS.P) %>% filter(genic_TWAS.P < 0.05) %>% filter(cis_TWAS.P > 0.05)

## which of these genes had not been identified from the GWAS analysis?

load('../GWAS_CATALOG/results/asrs_snps_gwas_add_direct_overlaps_ALL_DISEASE.RData')



asrs_trait_counts = asrs_snps_gwas_pval_dat %>% select(trait, gene_name) %>% unique() %>% count(trait) %>% rename(asrs_gene_count = n)

#load('../IMD_ref_files/immune_EFO.RData')
load('../../../gene_strand_fix/GWAS_CATALOG/gwas_catalog_IMD_fix/IMD_ref_files/immune_EFO.RData')
enrichment_results_fn = "../GWAS_CATALOG/results/empirical.test.enrichment.results10000perms.add.direct.and.total.counts.ANNOTATED"
enrichment_results = fread(enrichment_results_fn)

enrichment_results = inner_join(enrichment_results, asrs_trait_counts, by = 'trait')

immune_traits = enrichment_results %>% filter(is_immune)

ADJ_P = 0.05
GENE_CUTOFF = 3
traits_of_interest = enrichment_results %>% filter(asrs_FDR < ADJ_P) %>% filter(is_disease) %>% filter(asrs_gene_count >= GENE_CUTOFF) %>% filter(is_immune)
traits_of_interest = traits_of_interest %>% filter(trait != 'acute myeloid leukemia')
sig_traits = traits_of_interest$trait
traits_of_interest$log10p = -log10(traits_of_interest$asrs_FDR)
GWAS_genes = filter(asrs_snps_gwas_pval_dat, trait %in% traits_of_interest$trait)

unique_TWAS_genes = setdiff(sig_genic_res$gene_name, GWAS_genes$gene_name)

GWAS_TWAS_overlap = intersect(sig_genic_res$gene_name, GWAS_genes$gene_name)

filter(sig_genic_res, gene_name %in% unique_TWAS_genes) %>% arrange(p.adj)
filter(sig_genic_res, gene_name %in% GWAS_TWAS_overlap) %>% arrange(p.adj)

#### check if any of these results had not previously been reported in GWAS analysis
gwas_EFO_sig = fread('~/project-gxxiao4/ASE_ASRS/gene_strand_fix/GWAS_CATALOG/gwas_catalog_IMD_fix/IMD_ref_files/gwas.catalog.sig.EFO.MAPPED.URI.FIVE.IMD.CLEANED.GWAS.CATALOG.INFO')
traits = c('inflammatory bowel disease', "Crohn's disease", "systemic lupus erythematosus", "multiple sclerosis",  "rheumatoid arthritis")




twas_sum_clean = TWAS_res
p_cutoff_frame = twas_sum_clean %>%
  filter(p.adj > 0.05) %>%
  summarise(p_cutoff = min(TWAS.P))

p_cutoff = p_cutoff_frame$p_cutoff

twas_sum_clean$bp_cum = as.numeric(twas_sum_clean$bp_cum)

# first, among the genic results only, how many significant gene-disease associations are there? 
twas_sum_clean %>% filter(type == 'genic') %>% filter(p.adj < 0.05) %>% select(gene_name, disease) %>% unique() # 29 gene-disease associations, 23 unique genes
twas_sum_clean %>% filter(type == 'genic') %>% select(gene_name) %>% unique()

# now compaare gene/cis models

plot_dat = twas_sum_clean %>% filter(type == 'genic') %>% 
  #filter(!tissue == 'Colon_Sigmoid') %>%
  #filter(gene_type == 'asrs') %>%
  mutate(sig = TWAS.P < p_cutoff) %>%
  mutate(disease = ifelse(sig, disease, 'non-sig')) %>%
  filter(!is.na(disease))


## load color palette 
load('../GWAS_CATALOG/scripts/trait.color.pal.RData')


plot_dat$disease[plot_dat$disease == 'Crohns_Disease'] = "Crohn's disease"
plot_dat$disease[plot_dat$disease == 'inflammatory_bowl_disease'] = "Inflammatory bowel disease"
plot_dat$disease[plot_dat$disease == 'multiple_sclerosis'] = "Multiple sclerosis"
plot_dat$disease[plot_dat$disease == 'rheumatoid_arthritis'] = "Rheumatoid arthritis"
plot_dat$disease[plot_dat$disease == 'systemic_lupus_erythematosus'] = "Systemic lupus erythematosus"

color_lookup = filter(color_lookup, color != 'grey90') %>% filter(trait %in% plot_dat$disease)

# add shape per tissue
## add shape for each tissue
shape_lookup = twas_sum_clean %>% dplyr::select(tissue, disease) %>% unique() %>% arrange(disease)
shape_lookup$shapes = c(21, 22, 21, 22, 21, 22, 21, 23, 22, 23, 23, 21, 22 )
shape_lookup = shape_lookup %>% select(tissue, shapes) %>% unique()
plot_dat$tissue = factor(plot_dat$tissue, levels = shape_lookup$tissue)


color_lookup$trait = as.character(color_lookup$trait)
color_lookup = rbind(color_lookup, c('non-sig', 'grey92'))



plot_dat$disease = factor(plot_dat$disease, levels = color_lookup$trait)


axis_set <- plot_dat %>%
  #filter(!is.na(TWAS.P)) %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

g1 = ggplot(plot_dat,
            aes(x = bp_cum, y = -log10(as.numeric(TWAS.P) ),
                fill = disease,shape = tissue
            )) +
  geom_hline(yintercept = -log10(p_cutoff), color = "grey40", linetype = "dashed") +
  geom_point(aes(fill = disease, size = disease)) +
  geom_text_repel(aes(label = gene_name), size = 4, min.segment.length = 0,
                  data = plot_dat %>% filter(sig))  +
  #filter(TWAS.P < p_cutoff)) + 
  #scale_size_continuous(range = c(0.5,3))  + 
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  labs(x = NULL,
       y = "-log10(p)"
  ) +
  # guides(size = "none") +
  theme_classic() +
  #  theme(axis.text.x=element_text(size=6))  
  scale_shape_manual(values = shape_lookup$shapes) + scale_fill_manual(values = color_lookup$color) + scale_color_manual(values = color_lookup$color) + scale_size_manual(values = c(rep(3.5,5),2))

g2 = g1 +
  guides(fill=guide_legend(override.aes=list(shape=21, size = 3)))








pdf('../MANUSCRIPT_FIGURES/fig6/fig6d.pdf', width = 9.5, height = 4.5)
g2 + theme(legend.position = 'none') + xlab('Chromosome') # + theme(axis.text=element_text(size=7),
                                             #                                             axis.title=element_text(size=8))
                                             #pdf('../../../MANUSCRIPT_FIGURES/fig5/fig5d.pdf', width = 12, height = 5)
                                             #g2 + theme(axis.text=element_text(size=6),
                                            #axis.title=element_text(size=8)) + xlab('Chromosome')

dev.off()


