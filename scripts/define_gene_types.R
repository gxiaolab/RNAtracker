library(data.table)
library(dplyr)
library(tidyr)

ASRS_conditions = c('0-0-1', '0-1-1', '0-NA-1', '0-1-NA')
ASRT_conditions = c('1-1-1', '1-NA-1', '1-1-NA', '1-NA-NA')
mixed_conditions = c('1-0-0', '1-0-1', '1-0-NA', '1-NA-0', '1-1-0')

celltypes=c('GM12878', 'K562', 'HepG2','PC_3','Panc1','HCT116', 'MCF_7', 'IMR90', 'A673', 'Caco2', 'Calu3', 'MCF10A', 'OCILY7', 'PC9', 'HMEC', 'HUVEC')

load('../input_files/GENCODE_CNV_filtered.RData') # file with predicted cell line regions per cell line 

# change to wherever the output files from  RNAtracker_ase.py are saved
ase_dir = '../input_files/'

make_gene_label_table = function(cell_line){
	zero_hr_fn = paste0(ase_dir, cell_line, "/zero_hr/T10_M2_eiTrue_ARFalse_nsigma1.95996_p0.05_apfdr_bh_ASE")
    two_hr_fn = paste0(ase_dir, cell_line, "/two_hr/T10_M2_eiTrue_ARFalse_nsigma1.95996_p0.05_apfdr_bh_ASE")
    six_hr_fn = paste0(ase_dir, cell_line, "/six_hr/T10_M2_eiTrue_ARFalse_nsigma1.95996_p0.05_apfdr_bh_ASE")
	remove_genes = all_gencode_CNV_filtered_info_dat  %>% filter(dataset == cell_line) %>% filter(CNV_filtered == 'TRUE')
	zero_hr = fread(zero_hr_fn)
	two_hr = fread(two_hr_fn)
	six_hr = fread(six_hr_fn)
	zero_hr$time = '0h'
	two_hr$time = '2h'
	six_hr$time = '6h'
	all_ase = rbind(zero_hr, two_hr, six_hr)
	all_ase = filter(all_ase, !Gene %in% remove_genes$Gene)
	all_ase_full = all_ase
	all_ase_full$dataset = cell_line
	all_ase$ASEornonASE[all_ase$ASEornonASE == 'ASE'] <- 1
	all_ase$ASEornonASE[all_ase$ASEornonASE == 'nonASE'] <- 0
	ase_transformed = dcast(data = all_ase, formula = Gene ~ time, value.var = "ASEornonASE")
	ase_transformed$label = paste(ase_transformed$`0h`, ase_transformed$`2h`, ase_transformed$`6h`, sep = '-')
	ASRS_events = ase_transformed %>% filter(label %in% ASRS_conditions)
	ASRT_events = ase_transformed %>% filter(label %in% ASRT_conditions)
	mixed_events = ase_transformed %>% filter(label %in% mixed_conditions)
	# label genes by gene type
	all_ase_full$gene_label = 'not categorized'
	all_ase_full$gene_label[all_ase_full$Gene %in% ASRS_events$Gene] = 'asRS'
	all_ase_full$gene_label[all_ase_full$Gene %in% ASRT_events$Gene] = 'asRT'
	all_ase_full$gene_label[all_ase_full$Gene %in% mixed_events$Gene] = 'mixed'
	# if you want SNV level info
	#all_ase_full = all_ase_full %>% separate_rows(SNVs, sep = ';') %>% filter(SNVs != "")  %>% separate(SNVs, into = c('coord', 'region', 'AR', 'ref_allelic_ratio', 'snv_AR_delta', 'snv_AR_sigma', 'major', 'minor', 'p_orig', 'p_adj', 'nraw', 'nnorm'), sep = ',', convert = TRUE) %>% separate(coord, into = c('chr', 'pos', 'strand'), sep = ':')
	return(all_ase_full)
}

for (i in 1:length(celltypes)){
  cell_line = celltypes[i]
  print(cell_line)
  gene_label_dat = make_gene_label_table(cell_line)  
  assign(paste0(cell_line, '_gene_info'), gene_label_dat)
}

all_gene_info_list = lapply(ls(pattern = "*_gene_info"), get)
all_gene_info_dat = do.call("rbind", all_gene_info_list)

save(all_gene_info_dat, file = '../all_RNAtracker_results.RData')














