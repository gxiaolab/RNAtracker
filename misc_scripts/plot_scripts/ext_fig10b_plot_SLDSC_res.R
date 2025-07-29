library(data.table)
library(dplyr)
library(ggplot2)


setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')

res_frame = fread('../results/SLDSC_result.tsv')

table(res_frame$disease)
plot_dat = filter(res_frame, !disease %in% c('ankylosing_spondylitis', 'T1D')) %>% filter(type != 'all')

plot_dat$disease[plot_dat$disease == 'Crohns_Disease'] = "Crohn's disease"
plot_dat$disease[plot_dat$disease == 'inflammatory_bowl_disease'] = "Inflammatory bowel disease"
plot_dat$disease[plot_dat$disease == 'multiple_sclerosis'] = "Multiple sclerosis"
plot_dat$disease[plot_dat$disease == 'rheumatoid_arthritis'] = "Rheumatoid arthritis"
plot_dat$disease[plot_dat$disease == 'psoriasis'] = "Psoriasis"
plot_dat$disease[plot_dat$disease == 'asystemic_lupus_erythematosus'] = "Systemic lupus erythematosus"

plot_dat$type[plot_dat$type == 'ASRS'] = 'asRS'
plot_dat$type[plot_dat$type == 'ASRT'] = 'asRT'
plot_dat$type[plot_dat$type == 'non_ASE'] = 'non-ASE'



pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'non-ASE' = 'grey')


p = ggplot(plot_dat, aes(color = type)  )+
  geom_point(aes(y  = type, x = Enrichment)) +
  geom_errorbar(aes(xmin=Enrichment - Enrichment_std_error,
                      xmax=Enrichment + Enrichment_std_error,
                      y = type), width = 0.5) +
  facet_wrap(~disease, scales = 'free') +
  coord_flip() + 
  theme_classic() + ylab('Variant set') + scale_color_manual(values = pal)


pdf('../MANUSCRIPT_FIGURES/extended_fig10/extended_fig10b.pdf', height = 4, width = 7)
p + theme(legend.position = 'none') + xlim(c(-2,8))
dev.off()
