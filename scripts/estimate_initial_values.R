library(data.table)
library(dplyr)
library(VGAM)
source("~/project-gxxiao4/softwares/bb_mle.R") # this script can be obtained from: https://rdrr.io/github/andreaskapou/scMET/src/R/bb_mle.R
#' This function is from
#' [StackOverflow post 5577221](https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file)
loadRData <- function(file_name){
  load(file_name)
  get(ls()[ls() != "file_name"])
}



args = commandArgs(trailingOnly=TRUE)
COVERAGE_CUTOFF = args[1]
MISMATCH_CUTOFF = args[2]
out_suffix = paste0('GATK_VAR_merge_', COVERAGE_CUTOFF, '_', MISMATCH_CUTOFF)
P_ADJ_CUTOFF = 0.05
SNP_COUNT_FILTER = 2

output_dir = paste0('../REAL_DATA_RESULTS/pickle_files_FIX_MAJOR_ALLELE_', out_suffix, '/')


# load 0h BEAPR results


## function begin
read_data = function(CELL_LINE){
    timepoints = c('zero_hr', 'two_hr', 'six_hr')
    dat_list = list()
    for(timepoint in timepoints){
        out_fn = paste0(output_dir,CELL_LINE, '/', timepoint, '.filtered_testable_snvs_with_padj_info.txt')
        dat = fread(out_fn)
        dat$timepoint = timepoint 
        dat_list[[timepoint]] = dat
    }
    all_dat = do.call("rbind", dat_list)
    return(all_dat)
}

read_BEAPR_res = function(CELL_LINE){
    timepoints = c('zero_hr', 'two_hr', 'six_hr')
    dat_list = list()
    for(timepoint in timepoints){
        out_fn = paste0(output_dir, CELL_LINE,"/", timepoint, "/T", COVERAGE_CUTOFF, "_M", MISMATCH_CUTOFF, "_eiTrue_ARFalse_nsigma1.95996_p0.05_apfdr_bh_FIXED_MAJOR_ALLELE_ASE")
        dat = fread(out_fn)
        dat$timepoint = timepoint 
        dat_list[[timepoint]] = dat
    }
    all_dat = do.call("rbind", dat_list)
    return(all_dat)

}



make_estimates = function(CELL_LINE, include_timepoints, FILTER_BY_0h_TESTABLE = TRUE){
    print(CELL_LINE)
    all_data = read_data(CELL_LINE)
    all_data = filter(all_data, timepoint %in% include_timepoints)
    all_data = all_data %>% rename(id = V1)

    BEAPR_res = read_BEAPR_res(CELL_LINE)
    zero_hr_BEAPR_res = filter(BEAPR_res, timepoint == 'zero_hr')
    BEAPR_res = filter(BEAPR_res, timepoint %in% include_timepoints)


    ASE_SNV_level_res = BEAPR_res %>% tidyr::separate_rows(SNVs, sep = ';') %>% select(Gene, ASEornonASE, SNVs, timepoint) %>% filter(SNVs != "")  %>% tidyr::separate(SNVs, into = c('coord', 'region', 'AR', 'ref_allelic_ratio', 'delta_AR', 'sigma', 'major', 'minor', 'p_orig', 'p_adj', 'nraw', 'nnorm'), sep = ',', convert = TRUE) %>% tidyr::separate(coord, into = c('chr', 'pos', 'strand'), sep = ':')
    ASE_SNV_level_res = mutate(ASE_SNV_level_res, id = paste0(chr, ':', pos, ':', strand))
    sig_ASE = filter(ASE_SNV_level_res, p_adj < 0.05)

    # some snvs in all_data are not in the gene-level BEAPR results; these were removed later due to being intronic
    combine_dat = left_join(ASE_SNV_level_res %>% select(timepoint, id, Gene) %>% unique(), all_data, by = c('timepoint', 'id'))
    all_data = combine_dat %>% rename(gene_name = Gene)
    # require snvs to be at all timepoints
    #NUMBER_OF_TIMEPOINTS = length(unique(all_data$timepoint))
    #original_snv_count = length(unique(all_data$id))
    #pass_snvs = all_data %>% select(id, timepoint) %>% unique() %>% count(id) %>% filter(n==NUMBER_OF_TIMEPOINTS)
    #print(paste('originally we have', original_snv_count, 'snvs'))
    #all_data = filter(all_data, id %in% pass_snvs$id)
    #print(paste('after requiring snv to be at all include_timepoints, we have', nrow(pass_snvs), 'snvs'))

    # implement SNP count filter
    NUMBER_OF_TIMEPOINTS = length(unique(all_data$timepoint))
    print(NUMBER_OF_TIMEPOINTS)
    gene_snp_counts =  all_data %>% select(timepoint, gene_name, id) %>% unique() %>% group_by(timepoint, gene_name) %>% count(gene_name) %>% filter(n>=SNP_COUNT_FILTER)
    # require the gene to pass SNP_COUNT_FILTER for all timepoints 
    universal_genes = gene_snp_counts %>% group_by(gene_name) %>% count(gene_name) %>% filter(n==NUMBER_OF_TIMEPOINTS)
    all_data = filter(all_data, gene_name %in% universal_genes$gene_name)
    print(paste("after SNP_COUNT_FILTER, we have", length(unique(all_data$id)), 'snvs'))

    # gene must be testable at 0h
    if(FILTER_BY_0h_TESTABLE){
         all_data = all_data %>% filter(gene_name %in% zero_hr_BEAPR_res$Gene)
        print(paste("after requiring gene to have testable SNVs at 0h, we have", length(unique(all_data$id)), 'snvs'))
    }
   

   #all_data$rep1_major_allele_counts = pmax(all_data$rep1_ref_raw_counts, all_data$rep1_alt_raw_counts)
   #all_data$rep2_major_allele_counts = pmax(all_data$rep2_ref_raw_counts, all_data$rep2_alt_raw_counts)

    ## now define ASE and nonASE snvs based on the LOESS p-value
    all_data = all_data %>% select(-gene_name) %>% unique()
    ASE_SNVs = all_data %>% filter(adjusted_P < P_ADJ_CUTOFF)
    nonASE_SNVs = all_data %>% filter(adjusted_P >= P_ADJ_CUTOFF)

    #ASE_MAJOR_COUNT = c(ASE_SNVs$rep1_major_allele_counts,ASE_SNVs$rep2_major_allele_counts) # this still fails; probably because the range of allelic ratios is still between 0.5 and 0.9
    ASE_REF_COUNT = c(ASE_SNVs$rep1_ref_raw_counts,ASE_SNVs$rep2_ref_raw_counts)
    ASE_TOTAL_COUNT = c(ASE_SNVs$rep1_ref_raw_counts + ASE_SNVs$rep1_alt_raw_counts, ASE_SNVs$rep2_ref_raw_counts + ASE_SNVs$rep2_alt_raw_counts)
    NONASE_REF_COUNT = c(nonASE_SNVs$rep1_ref_raw_counts, nonASE_SNVs$rep2_ref_raw_counts)
    NONASE_TOTAL_COUNT = c(nonASE_SNVs$rep1_ref_raw_counts + nonASE_SNVs$rep1_alt_raw_counts, nonASE_SNVs$rep2_ref_raw_counts + nonASE_SNVs$rep2_alt_raw_counts)


    #fit_ASE_MAJOR = bb_mle(cbind(ASE_MAJOR_COUNT, ASE_TOTAL_COUNT),  n_starts = 10)
    fit_nonASE_REF = bb_mle(cbind(NONASE_TOTAL_COUNT, NONASE_REF_COUNT), w = c(1.1,1.1), n_starts = 10)
    print('done estimating nonASE hyperparameters')
    #print(warnings())
    fit_ASE_REF = bb_mle(cbind(ASE_TOTAL_COUNT, ASE_REF_COUNT), w = c(1.1,1.1), n_starts = 10)
     print('done estimating ASE hyperparameters')
    #print(warnings())
    out = data.frame(cell_line = CELL_LINE, a0 = fit_nonASE_REF$alpha, b0 = fit_nonASE_REF$beta, a1 = fit_ASE_REF$alpha, b1 = fit_ASE_REF$beta, non_ASE_conv_status = fit_nonASE_REF$is_conv, ASE_conv_status = fit_ASE_REF$is_conv)   
    return(out)
}

#cell_lines = c("A673", "Caco-2", "Calu3", "GM12878", "HCT116", "HepG2", "HMEC", "HUVEC", "IMR-90", "K562", "MCF10A", "MCF-7", "OCI-LY7", "PC-3", "PC-9") # no Panc1 because no passing genes

cell_lines = c("A673", "Caco-2", "Calu3", "GM12878", "HCT116", "HepG2", "HMEC", "HUVEC", "IMR-90", "MCF10A", "MCF-7", "OCI-LY7", "PC-3", "PC-9")


all_estimates_only_0h_no_0h_filter = do.call("rbind", lapply(cell_lines, make_estimates, c('zero_hr'), FALSE))
estimates_fn = paste0('../estimated_initial_values_RECOVERED_SNVs_INCLUDE_ZERO_', out_suffix, '_ONLY_0h_FIXED_ALLELE_NO_0h_TESTABLE_REQUIREMENT')
write.table(all_estimates_only_0h_no_0h_filter, file = estimates_fn, sep = '\t', row.names = F, quote = F)



all_estimates_df = do.call("rbind", lapply(cell_lines, make_estimates, c('two_hr', 'six_hr')))


estimates_fn = paste0('../estimated_initial_values_RECOVERED_SNVs_INCLUDE_ZERO_', out_suffix, '_REMOVE_0h_FIXED_ALLELE')
#estimates_fn = paste0('../estimated_initial_values_RECOVERED_SNVs_INCLUDE_ZERO_', out_suffix)
write.table(all_estimates_df, file = estimates_fn, sep = '\t', row.names = F, quote = F)


all_estimates_include_0h_df = do.call("rbind", lapply(cell_lines, make_estimates, c('zero_hr', 'two_hr', 'six_hr')))

estimates_fn = paste0('../estimated_initial_values_RECOVERED_SNVs_INCLUDE_ZERO_', out_suffix, '_INCLUDE_0h_FIXED_ALLELE')
#estimates_fn = paste0('../estimated_initial_values_RECOVERED_SNVs_INCLUDE_ZERO_', out_suffix)
write.table(all_estimates_include_0h_df , file = estimates_fn, sep = '\t', row.names = F, quote = F)


## now do the estimates without requiring the snv to be testable at 0h
all_estimates_no_0h_filter = do.call("rbind", lapply(cell_lines, make_estimates, c('two_hr', 'six_hr'), FALSE))
estimates_fn = paste0('../estimated_initial_values_RECOVERED_SNVs_INCLUDE_ZERO_', out_suffix, '_REMOVE_0h_FIXED_ALLELE_NO_0h_TESTABLE_REQUIREMENT')
write.table(all_estimates_no_0h_filter , file = estimates_fn, sep = '\t', row.names = F, quote = F)

all_estimates_inclde_0h_no_0h_filter = do.call("rbind", lapply(cell_lines, make_estimates, c('zero_hr', 'two_hr', 'six_hr'), FALSE))
estimates_fn = paste0('../estimated_initial_values_RECOVERED_SNVs_INCLUDE_ZERO_', out_suffix, '_INCLUDE_0h_FIXED_ALLELE_NO_0h_TESTABLE_REQUIREMENT')
write.table(all_estimates_inclde_0h_no_0h_filter, file = estimates_fn, sep = '\t', row.names = F, quote = F)

## now make estimates using 0h timepoint only
# TRUE OR FALSE for FILTER_BY_0H_TESTABLE parameter should yield same result
test1 = make_estimates('GM12878', 'zero_hr', TRUE)
test2 = make_estimates('GM12878', 'zero_hr', FALSE)
# confirmed identical

print(test1)
print(test2)

