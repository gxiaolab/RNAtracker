library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)


CELL_LINE = 'K562'
all_data = fread('../Bru-vs-NoBru_from_Mats.txt')


m_dat = all_data %>% select(ensembl, name, rpkm1, rpkm2) %>% unique() %>% melt(id.vars = c('ensembl', 'name'))
m_dat$variable = as.character(m_dat$variable)
m_dat$variable[m_dat$variable == 'rpkm1'] = 'Bru'
m_dat$variable[m_dat$variable == 'rpkm2'] = 'unlabeled'
names(m_dat) = c('ensembl_id', 'gene_name', 'condition', 'RPKM')

all_data = all_data %>% filter(log(rpkm1)>0) %>% filter(log(rpkm2)>0)



test =cor.test(all_data$rpkm1, all_data$rpkm2)
print(test)

label = paste0("\nPearson's Correlation = ", round(test$estimate, 3), "\nP <  2.2e-16")
linear_model = lm(log(rpkm1) ~ log(rpkm2), all_data)
all_data$residuals = as.data.frame(linear_model$residuals)
all_data$outlier = abs(all_data$residuals) >= 1
cor.test(all_data$rpkm1, all_data$rpkm2)

out = all_data %>% select(name, ensembl, rpkm1, rpkm2, outlier) %>% rename(bru_RPKM = rpkm1, unlabeled_RPKM = rpkm2)
write.table(out, file = '../../../SOURCE_DATA/EXTENDED_DATA_FIGS/extended_fig1a.txt', sep = '\t', row.names = F, quote = F)



p = ggplot(all_data, aes(x = log(rpkm1), y = log(rpkm2))) + geom_point(aes(color = outlier)) + theme_classic() + 
  xlab('Bru-labeled RPKM') + ylab('unlabeled RPKM') + geom_smooth() + scale_color_manual(values = c('black', 'red')) + ggtitle(label)

  #geom_label_repel(data = filter(all_data, outlier), aes(label = name)) + ggtitle(label)

pdf('../figures/ljungman_lab_K562_bru_vs_nobru_scatter.pdf', width = 9)
p + theme(legend.position = "none")
dev.off()

pdf('../MANUSCRIPT_FIGURES/extended_fig1/extended_fig1a.pdf', width = 4, height = 3.5)
p + theme(legend.position = "none")
dev.off()

