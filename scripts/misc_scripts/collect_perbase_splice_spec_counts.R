library(data.table)
library(dplyr)


perbase_files = list.files('../perbase_splice_spec/', full.names = TRUE)

counts = list()

for (i in 1:length(perbase_files)){
    filename = perbase_files[i]
    print(filename)
    dat = fread(filename, sep = '\t')
    dat$label = basename(filename)
    counts[[i]] = dat
}

full_counts_df = do.call("rbind", counts)


final_out = list()
# get allelic ratio 
REF_ALLELE = 'T'
ALT_ALLELE = 'C'
counts_df = full_counts_df %>% filter(POS == 173364016)

counts_df = counts_df %>% tidyr::separate(label, into = c('label1', 'label2', 'timepoint', 'replicate_number'), sep = '-') %>% tidyr::separate(replicate_number, into = c('replicate_number', 'sample_number'), sep = '_') %>% tidyr::separate(sample_number, into = c('sample_number', 'read_type', 'extra'), sep = '\\.')

counts_df = counts_df %>% mutate(ref_allelic_ratio = T/(T+C))
final_out[[1]] = counts_df

# maybe try to do it for the other variants too

## EPHA3
REF_ALLELE = 'C'
ALT_ALLELE = 'A'
counts_df = full_counts_df %>% filter(POS == 89480907)

counts_df = counts_df %>% tidyr::separate(label, into = c('label1', 'label2', 'timepoint', 'replicate_number'), sep = '-') %>% tidyr::separate(replicate_number, into = c('replicate_number', 'sample_number'), sep = '_') %>% tidyr::separate(sample_number, into = c('sample_number', 'read_type', 'extra'), sep = '\\.')

counts_df = counts_df %>% mutate(ref_allelic_ratio = C/(C+A))
final_out[[2]] = counts_df

## RNF114
REF_ALLELE = 'G'
ALT_ALLELE = 'T'
counts_df = full_counts_df %>% filter(POS == 49952392)

counts_df = counts_df %>% tidyr::separate(label, into = c('label1', 'label2', 'timepoint', 'replicate_number'), sep = '-') %>% tidyr::separate(replicate_number, into = c('replicate_number', 'sample_number'), sep = '_') %>% tidyr::separate(sample_number, into = c('sample_number', 'read_type', 'extra'), sep = '\\.')

counts_df = counts_df %>% mutate(ref_allelic_ratio = G/(G+T))
final_out[[3]] = counts_df

# CD151
REF_ALLELE = 'G'
ALT_ALLELE = 'T'
counts_df = full_counts_df %>% filter(POS == 834745)

counts_df = counts_df %>% tidyr::separate(label, into = c('label1', 'label2', 'timepoint', 'replicate_number'), sep = '-') %>% tidyr::separate(replicate_number, into = c('replicate_number', 'sample_number'), sep = '_') %>% tidyr::separate(sample_number, into = c('sample_number', 'read_type', 'extra'), sep = '\\.')

counts_df = counts_df %>% mutate(ref_allelic_ratio = G/(G+T))
final_out[[4]] = counts_df

REF_ALLELE = 'C'
ALT_ALLELE = 'T'
counts_df = full_counts_df %>% filter(POS == 838672)

counts_df = counts_df %>% tidyr::separate(label, into = c('label1', 'label2', 'timepoint', 'replicate_number'), sep = '-') %>% tidyr::separate(replicate_number, into = c('replicate_number', 'sample_number'), sep = '_') %>% tidyr::separate(sample_number, into = c('sample_number', 'read_type', 'extra'), sep = '\\.')

counts_df = counts_df %>% mutate(ref_allelic_ratio = C/(C+T))
final_out[[5]] = counts_df



all_counts_df = do.call("rbind",final_out)
write.table(all_counts_df, file = '../results/perbase_splice_spec_counts.txt', sep = '\t', row.names = FALSE, quote = F)
