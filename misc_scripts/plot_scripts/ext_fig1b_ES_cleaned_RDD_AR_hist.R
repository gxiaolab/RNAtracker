library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
#COVERAGE_CUTOFF = args[1]
#MISMATCH_CUTOFF = args[2]



COVERAGE_CUTOFF = 10
MISMATCH_CUTOFF = 2
out_suffix = paste0('GATK_VAR_merge_', COVERAGE_CUTOFF, '_', MISMATCH_CUTOFF)
out_dir = paste0("../REAL_DATA_RECOVERED_SNVS_", COVERAGE_CUTOFF, "_", MISMATCH_CUTOFF, "/")



get_dat_from_RDD = function(CELL_LINE){
	final_out = list()
	timepoints = c('zero_hr', 'two_hr', 'six_hr')

	for(timepoint in timepoints) {
		if(timepoint == 'zero_hr'){
		    timepoint_numeric = '0h'
		}

		if(timepoint == 'two_hr'){
		    timepoint_numeric = '2h'
		}

		if(timepoint == 'six_hr'){
		    timepoint_numeric = '6h'
		}


		out_fn1 = paste0('../ES_CN_CLEANED_RDD/', CELL_LINE, '/', timepoint, '/', CELL_LINE, '.MERGED.PASS_ONLY.het_snps.final.all_snvs.',CELL_LINE, '_', timepoint_numeric,'_rep1.rdd.annot')
		out_fn2 = paste0('../ES_CN_CLEANED_RDD/', CELL_LINE, '/', timepoint, '/', CELL_LINE, '.MERGED.PASS_ONLY.het_snps.final.all_snvs.',CELL_LINE, '_', timepoint_numeric,'_rep2.rdd.annot')


		dat1 = fread(out_fn1)
		dat2 = fread(out_fn2)
		dat1 = dat1 %>% tidyr::separate(V5, into = c('ref', 'alt', 'other'), sep = ':', convert = TRUE) %>% mutate(total = ref+alt) %>% mutate(AR = ref/total)

		dat2 = dat2 %>% tidyr::separate(V5, into = c('ref', 'alt', 'other'), sep = ':', convert = TRUE) %>% mutate(total = ref+alt) %>% mutate(AR = ref/total)


		dat1 = dat1 %>% mutate(id = paste0(V1, ':', V2, ':', V3, ':', V6))
		dat2 = dat2 %>% mutate(id = paste0(V1, ':', V2, ':', V3, ':', V6))
		out_dat = rbind(dat1 %>% select(id, ref, alt, total, AR) %>% unique() %>% mutate(rep = 'rep1'), dat2 %>% select(id, ref, alt, total, AR) %>% unique() %>% mutate(rep = 'rep2'))
		out_dat$timepoint = timepoint
		out_dat$timepoint_numeric = timepoint_numeric
		final_out[[timepoint]] = out_dat
	}
	final_df = do.call("rbind", final_out)
	final_df$cell_line = CELL_LINE
	p = ggplot(final_df, aes(x = AR, fill = timepoint_numeric)) + geom_histogram(color = 'black', alpha = 0.67) + theme_bw() + ggtitle(CELL_LINE)

	#number_genes = final_df %>% filter(ref >=2 & alt >=2 & total>=5) %>% select(V7) %>% unique() %>% length()
	p2 = ggplot(final_df %>% filter(ref >=2 & alt >=2 & total>=5), aes(x = AR, fill = timepoint_numeric)) + geom_histogram(color = 'black', alpha = 0.67) + theme_bw() + ggtitle(paste0(CELL_LINE, ' C>=5, M>=2'))	
	print(p)
	print(p2)
	return(final_df)
}

get_dat = function(CELL_LINE){
	out_fn = paste0(out_dir, CELL_LINE, '.ALL_TIMEPOINTS.RECOVER_SNVs_INCLUDE_ZERO.RData')
	load(out_fn) # all_timepoints_complete
	zero_hr = all_timepoints_complete %>% select(id, contains('0h')) 
	two_hr = all_timepoints_complete %>% select(id, contains('2h')) 
	six_hr = all_timepoints_complete %>% select(id, contains('6h')) 
	zero_hr$timepoint_numeric = '0h'
	two_hr$timepoint_numeric = '2h'
	six_hr$timepoint_numeric = '6h'
	names(zero_hr) = c('id', 'p_orig', 'p_adj', 'MODEL_COUNT', 'total_counts', 'timepoint', 'timepoint_numeric')
	names(two_hr) = c('id', 'p_orig', 'p_adj', 'MODEL_COUNT', 'total_counts', 'timepoint', 'timepoint_numeric')
	names(six_hr) = c('id', 'p_orig', 'p_adj', 'MODEL_COUNT', 'total_counts', 'timepoint', 'timepoint_numeric')
	all_dat = rbind(zero_hr, two_hr, six_hr) %>% unique()
	all_dat = mutate(all_dat, AR = MODEL_COUNT/total_counts)
	p = ggplot(all_dat, aes(x = AR, fill = timepoint_numeric)) + geom_histogram(color = 'black', alpha = 0.67) + theme_classic() + ggtitle(CELL_LINE) + theme(legend.position = 'none')
	#print(p)
	return(p + geom_vline(xintercept = 0.5, linetype = 'dashed', color = 'red'))
}

get_dat_legend = function(CELL_LINE){
  out_fn = paste0(out_dir, CELL_LINE, '.ALL_TIMEPOINTS.RECOVER_SNVs_INCLUDE_ZERO.RData')
  load(out_fn) # all_timepoints_complete
  zero_hr = all_timepoints_complete %>% select(id, contains('0h')) 
  two_hr = all_timepoints_complete %>% select(id, contains('2h')) 
  six_hr = all_timepoints_complete %>% select(id, contains('6h')) 
  zero_hr$timepoint_numeric = '0h'
  two_hr$timepoint_numeric = '2h'
  six_hr$timepoint_numeric = '6h'
  names(zero_hr) = c('id', 'p_orig', 'p_adj', 'MODEL_COUNT', 'total_counts', 'timepoint', 'timepoint_numeric')
  names(two_hr) = c('id', 'p_orig', 'p_adj', 'MODEL_COUNT', 'total_counts', 'timepoint', 'timepoint_numeric')
  names(six_hr) = c('id', 'p_orig', 'p_adj', 'MODEL_COUNT', 'total_counts', 'timepoint', 'timepoint_numeric')
  all_dat = rbind(zero_hr, two_hr, six_hr) %>% unique()
  all_dat = mutate(all_dat, AR = MODEL_COUNT/total_counts)
  p = ggplot(all_dat, aes(x = AR, fill = timepoint_numeric)) + geom_histogram(color = 'black', alpha = 0.67) + theme_bw() + ggtitle(CELL_LINE) 
  #print(p)
  return(p)
}




#CELL_LINES = fread('../cell_lines_noPanc1', header = F)
CELL_LINES = fread('../cell_lines_11_final', header = F)
all_dat = lapply(CELL_LINES$V1, get_dat)

pdf('../figures/AR_histo_CN_ES_filtered.pdf', width = 18, height = 11)

plot_grid(plotlist = all_dat)

dev.off()


pdf('../MANUSCRIPT_FIGURES/extended_fig1/extended_fig1c.pdf', width = 16, height = 6.5)
plot_grid(plotlist = all_dat, nrow = 2)
dev.off()

#pdf('../figures/AR_histo_CN_ES_CACO2_BERG.pdf')
#get_dat('CACO2_BERG')
#dev.off()














