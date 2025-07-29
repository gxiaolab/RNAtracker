library(data.table)
library(dplyr)
library(ggplot2)


dat = fread('../results/asrs_functional_variants_simplified_info.strand.fixed.txt')
dat$variant_type = 'asRS'

dat$val_method[dat$val_method == 'MPRAu_HEK293FT'] <- 'MPRAu'
dat$val_method[dat$val_method == 'MPRAu_HEPG2'] <- 'MPRAu'
dat$val_method[dat$val_method == 'MPRAu_HMEC'] <- 'MPRAu'
dat$val_method[dat$val_method == 'MPRAu_K562'] <- 'MPRAu'
dat$val_method[dat$val_method == 'MPRAu_SKNSH'] <- 'MPRAu'
dat$val_method[dat$val_method == 'MPRAu_GM12878'] <- 'MPRAu'
dat$val_method[dat$val_method == 'old_MAPUTR_HEK293'] <- 'old_MAPUTR'
dat$val_method[dat$val_method == 'old_MAPUTR_HELA'] <- 'old_MAPUTR'


dat$val_method[dat$val_method == 'new_MAPUTR'] <- 'MapUTR'
dat$val_method[dat$val_method == 'old_MAPUTR'] <- 'MapUTR'

dat = rename(dat, `Validation Method` = val_method)


p = ggplot(data = dat %>% filter(variant_type == 'asRS') %>% unique(), aes(x = variant_type, fill = `Validation Method`)) + geom_bar(position = "stack")  + labs(title = 'Functionally validated as-RS SNVs')  + scale_fill_brewer(palette = "Accent") + xlab("validation method") + theme_bw()


# make donut chart
dat = dat %>% select(-gene_name) %>% unique()
count_dat = dat %>% group_by(`Validation Method`) %>% count()
count_dat$fraction = count_dat$n/sum(count_dat$n)
# Compute the cumulative percentages (top of each rectangle)
count_dat$ymax = cumsum(count_dat$fraction)

# Compute the bottom of each rectangle
count_dat$ymin = c(0, head(count_dat$ymax, n=-1))

# Compute label position and make label
count_dat$labelPosition <- (count_dat$ymax + count_dat$ymin) / 2
count_dat$label = paste0(count_dat$`Validation Method`, ": ", count_dat$n)
count_dat$`Validation Method` = factor(count_dat$`Validation Method`, levels = c('ASB', 'MPRAu', 'MapUTR'))

p2 = ggplot(count_dat, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill = `Validation Method`)) +
  geom_rect() +
  geom_label( x=3.5, fill = 'white', aes(y=labelPosition, label=label), size=10) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

p3 = ggplot(count_dat, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill = `Validation Method`)) +
  geom_rect(color = 'black') +
  geom_label( x=3.5, aes(y=labelPosition, label=n, fill = `Validation Method`), size=6) +
  scale_fill_brewer(palette="Dark2") +
  coord_polar(theta="y", start = -90) +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

count_dat = count_dat %>% 
  arrange(desc(`Validation Method`)) %>%
  mutate(prop = n / sum(count_dat$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

count_dat$ypos[3] = 0

library(PNWColors)
pal=pnw_palette("Starfish",3, type = "discrete")
p4 = ggplot(count_dat, aes(x = "", y = prop, fill = `Validation Method`)) + geom_bar(alpha = 0.67, stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_label(aes(label = label), position = position_stack(vjust=0.5), fill = "white", size=6) +
  scale_fill_manual(values = pal)


count_dat = count_dat %>% mutate(pos = mean(ymax, ymin))

ggplot(count_dat, aes(x = "", y = prop, fill = `Validation Method`)) + geom_bar(alpha = 0.67, stat="identity", width=1, color="black") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 8) + coord_polar("y", start=0) +
  theme_void()   +
  scale_fill_manual(values = pal)



pdf('../figures/func_validated_asrs.pdf')
p
p2
p3
dev.off()


