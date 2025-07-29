library(data.table)
library(dplyr)
library(ggplot2)

setwd('~/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')

load('../results/ASRS_complete_functional_barplot.RData')

p1 = ggplot(all_val_summary, aes(x = general_label)) + geom_bar(aes(fill = label), color = 'black') + theme_bw() + scale_fill_brewer(palette = 'Dark2') + ylab('Number of SNVs') + xlab('Functional validation method')

p2 = ggplot(all_val_summary %>% select(-id) %>% unique(), aes(x = general_label)) + geom_bar(aes(fill = label), color = 'black') + theme_bw() + scale_fill_brewer(palette = 'Dark2') + ylab('Number of genes') + xlab('Functional validation method')


library(cowplot)


pdf('../figures/complete_functiona_validation_methods_barplot.pdf', width = 10, height = 6)
plot_grid(p1, p2)
dev.off()

# make donut chart

count_dat = donut_plot_dat_max %>% group_by(label) %>% count()
count_dat = rename(count_dat, `Validation Method` = label)
count_dat$fraction = count_dat$n/sum(count_dat$n)
# Compute the cumulative percentages (top of each rectangle)
count_dat$ymax = cumsum(count_dat$fraction)

# Compute the bottom of each rectangle
count_dat$ymin = c(0, head(count_dat$ymax, n=-1))

# Compute label position and make label
count_dat$labelPosition <- (count_dat$ymax + count_dat$ymin) / 2
count_dat$label = paste0(count_dat$`Validation Method`, ": ", count_dat$n)
#count_dat$`Validation Method` = factor(count_dat$`Validation Method`, levels = c('ASB', 'MPRAu', 'MapUTR'))



p3 = ggplot(count_dat, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill = `Validation Method`)) +
  geom_rect(color = 'black') +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette="Set2") +
  coord_polar(theta="y", start = -90) +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

pdf('../MANUSCRIPT_FIGURES/fig4/fig4f.pdf')
p3
dev.off()



