library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
counts = fread('../results/perbase_splice_spec_counts.txt')

counts$timepoint = factor(counts$timepoint, levels = c('0h', '2h', '8h', '24h'))

# we only want to look at spliced reads
counts = filter(counts, read_type == 'spliced')
SIG_POS = c(173364016, 834745 )
count_dat = filter(counts, POS %in% SIG_POS)

#sig_vars = c("CDCA7:chr2:173364016:+:T.to.C", "CD151:chr11:834745:+:G.to.T")

count_dat$variant_id = ""
count_dat$variant_id[count_dat$POS == 173364016] = "CDCA7:chr2:173364016:T>C"
count_dat$variant_id[count_dat$POS == 834745] = "CD151:chr11:834745:G>T"
count_dat$variant_id = factor(count_dat$variant_id, levels = c('CDCA7:chr2:173364016:T>C', 'CD151:chr11:834745:G>T'))
count_dat$alt_allelic_ratio = 1-count_dat$ref_allelic_ratio

baseline_dat = count_dat %>% filter(timepoint == '0h') %>% group_by(variant_id) %>% summarize(baseline = mean(alt_allelic_ratio))
count_dat = left_join(count_dat, baseline_dat, by = 'variant_id')
my_comparisons = list(c('0h', '2h'), c('0h', '8h'), c('0h', '24h'))

count_dat = count_dat %>% select(variant_id, timepoint, alt_allelic_ratio,baseline)

write.table(count_dat, file = "../../../SOURCE_DATA/MAIN_FIGS/fig4h.txt", sep = '\t', row.names = F, quote = F)



p1 = ggplot(count_dat, aes(x = timepoint, y = alt_allelic_ratio, fill = timepoint)) +
  geom_boxplot(width = 0.5) + facet_grid(.~variant_id) + stat_compare_means(comparisons = my_comparisons, method = 't.test', label = "p.signif", hide.ns = TRUE, vjust = 0.3, tip.length = 0.015) + 
  theme_classic() + ylab('Alternative Allelic Ratio') + xlab('Timepoint') + scale_fill_manual(values = c('grey', brewer.pal(3, 'Purples'))) + geom_hline(aes(yintercept = baseline), linetype = 'dashed', color = 'grey') #+ geom_point(position=position_dodge(width=0.75),aes(group=allele))

pdf('../figures/validated_2_variants_allelic_ratio_perbase_spliced_counts.pdf', width = 5, height = 3)
p1 + theme(legend.position = 'none')
dev.off()
