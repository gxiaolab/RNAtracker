library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(sjmisc)


#setwd('/Users/elainehuang/Desktop/project-gxxiao4.2/ASE_ASRS/REVISION2/scripts')

ldsc_files = list.files('../SUMSTATS/GWAS_SUMSTATS/ldsc', pattern = ".results", full.names = TRUE)

all_ldsc = list()
for (i in 1:length(ldsc_files)){
	dat = fread(ldsc_files[[i]])
	dat$label = basename(ldsc_files[[i]])
	all_ldsc[[i]] = dat
}


all_ldsc_df = do.call("rbind", all_ldsc)
all_ldsc_df$id = gsub(".results", "", all_ldsc_df$label)

## asrs significant
sig_asrs = filter(all_ldsc_df, Category == 'L2_0') %>% filter(Enrichment_p < 0.05) 

interesting_results = filter(all_ldsc_df, label %in% sig_asrs$label) 
# weird cases
# 29769526-GCST006938-EFO_0009585.results, 34737426-GCST90041799-EFO_1001214: no genome-wide significant snps
# GCST90041799: "A generalized linear mixed model association tool for biobank-scale data."



# match ldsc results with study details
dat = fread('../EFO_0000540_studies_export.tsv')
harmonised_links = fread('../harmonised_list.txt', header = FALSE)
dat$study_folder = gsub('http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/', '', dat$summaryStatistics)
harmonised_links$study_folder = gsub("/harmonised/.*", "", gsub("\\./", "", harmonised_links$V1))
immune_harmonised_links = left_join(dat, harmonised_links, by = 'study_folder')
immune_harmonised_links$id = gsub(".h.tsv.gz", "", basename(immune_harmonised_links$V1))
interesting_results = left_join(all_ldsc_df, immune_harmonised_links, by = 'id')

## check sample size
interesting_results$N = unlist(lapply(stringr::str_extract_all(interesting_results$discoverySampleAncestry, "\\d+"),  function(x) sum(as.numeric(x))))

#interesting_results %>% filter(N>=5000) %>% filter(genotypingTechnologies == 'Genome-wide genotyping array') %>% filter(Category == 'L2_0')  %>% filter(Enrichment >= Enrichment_std_error) %>% filter(Enrichment >= 1) %>% arrange(desc(Enrichment)) %>% filter(title!= "A generalized linear mixed model association tool for biobank-scale data. ") %>% select(Enrichment, Enrichment_p, reportedTrait, efoTraits) %>% filter(!efoTraits %like% ",") %>% count(efoTraits) %>% arrange(desc(n))
interesting_results$type = ""
interesting_results$type[interesting_results$Category == "L2_0"] = "asRS"
interesting_results$type[interesting_results$Category == "L2_1"] = "asRT"
interesting_results$type[interesting_results$Category == "L2_2"] = "non-ASE"
interesting_results$type[interesting_results$Category == "L2_3"] = "all"

pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'non-ASE' = 'grey')


interesting_study_info = interesting_results %>% select(title, journal, reportedTrait, pubmedId,id) %>% unique()
write.table(interesting_study_info, file = '../sldsc_study_info.txt', sep = '\t', row.names = F, quote = F)

interesting_results$reportedTrait = paste0(interesting_results$reportedTrait, '\n', interesting_results$id)
p1 = ggplot(interesting_results, aes(color = type)  )+
  geom_point(aes(y  = type, x = Enrichment)) +
  geom_errorbar(aes(xmin=Enrichment - Enrichment_std_error,
                      xmax=Enrichment + Enrichment_std_error,
                      y = type), width = 0.5) +
  facet_wrap(~reportedTrait) +
  coord_flip() +
  theme_classic() + ylab('Variant set') + scale_color_manual(values = pal)

  p2 = ggplot(interesting_results %>% filter(type!='all'), aes(color = type)  )+
  geom_point(aes(y  = type, x = Enrichment)) +
  geom_errorbar(aes(xmin=Enrichment - Enrichment_std_error,
                      xmax=Enrichment + Enrichment_std_error,
                      y = type), width = 0.5) +
  facet_wrap(~reportedTrait) +
  coord_flip() +
  theme_classic() + ylab('Variant set') + scale_color_manual(values = pal)

pdf('../figures/ldsc_enrichment.pdf', height = 11, width = 15)
p1 + theme(legend.position = 'none') #+ xlim(c(-2,8))
p2 + theme(legend.position = 'none') #+ xlim(c(-2,8))
dev.off()



## plot like a heatmap

# apply filters
SAMPLE_SIZE = 5000



asrs_res = interesting_results %>% filter(type=='asRS') %>% filter(N>=SAMPLE_SIZE) %>% filter(Enrichment >= Enrichment_std_error) %>% filter(Enrichment_p <= 0.1) %>% filter(!efoTraits%like%",")
asrs_res$Category = 'asRS'
write.table(asrs_res, file = '../SUPP_TABLE_SX2_ASRS_SLDSC.txt', sep = '\t', row.names = F, quote = F)

write.table(asrs_res, file = '../../SOURCE_DATA/MAIN_FIGS/fig6b.txt', sep = '\t', row.names = F, quote = F)


#asrs_res$Enrichment_p = gtools::stars.pval(asrs_res$Enrichment_p)
enrich_dat = asrs_res %>% arrange(desc(Enrichment)) %>%  select(efoTraits, id, type, Enrichment) 
p_dat = asrs_res %>% arrange(desc(Enrichment)) %>% select(efoTraits, id, type, Enrichment_p)

enrich_dat2 = dcast(data = enrich_dat, formula = id ~ efoTraits, value.var = "Enrichment")
p_dat2 = dcast(data = p_dat, formula = id ~ efoTraits, value.var = "Enrichment_p") 


enrich_mat = as.matrix(enrich_dat2[,-1])
row.names(enrich_mat) = gsub("_", " ", enrich_dat2$id)
row.names(enrich_mat) =  word_wrap(row.names(enrich_mat), wrap = 15)
p_mat = as.matrix(p_dat2[,-1])
row.names(p_mat) = p_dat2$id


# code according to gtools
p_val_func = function(j, i, x, y, w, h, fill) {
  if (is.na(p_mat[i,j])){
    grid.text("", x,y)
  }
  else if(p_mat[i, j] <= 0.001 & enrich_mat[i,j] > 1) {
    grid.text("***", x, y, gp = gpar(fontface = "bold", fontsize = 11))
  } else if(p_mat[i, j] <= 0.01 & enrich_mat[i,j] > 1) {
    grid.text("**", x, y, gp = gpar(fontface = "bold", fontsize = 11))
  } else if (p_mat[i,j] <= 0.05 & enrich_mat[i,j] > 1){
    grid.text("*", x, y, gp = gpar(fontface = "bold", fontsize = 11))
  #} else if (p_mat[i,j] <= 0.1 & enrich_mat[i,j] > 1){
  #  grid.text(".", x, y, gp = gpar(fontface = "bold", fontsize = 11))
  }  else {
    grid.text("", x,y)
  }
  
}


p_mat = t(p_mat)
enrich_mat = t(enrich_mat)
col_fun = colorRamp2(range(enrich_mat, na.rm = TRUE), hcl_palette = "Sunset", reverse = TRUE)

rank_diseases = sort(rowMeans(enrich_mat, na.rm = TRUE), decreasing = TRUE)
enrich_mat = enrich_mat[names(rank_diseases),]
p_mat = p_mat[names(rank_diseases),]



pdf('../figures/ldsc_asrs_enrichments_for_immune_hits_gwas.pdf', width = 9, height = 3)
Heatmap(enrich_mat, name = 'Heritability Enrichment', cell_fun = p_val_func, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE)
dev.off()











