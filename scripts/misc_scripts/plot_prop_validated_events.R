library(data.table)
library(dplyr)
library(ggplot2)


get_complete_results = function(CELL_LINE){
    in_data =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV_ALL_STATE_GENES.RData')
    load(in_data)
    all_results$cell_line = CELL_LINE 
    return(all_results)
}

get_complete_results_old = function(CELL_LINE){
    in_data =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV.RData')
    load(in_data)
    all_dat$cell_line = CELL_LINE 
    return(all_dat)
}

HCT116 = get_complete_results('HCT116')
GM12878 = get_complete_results('GM12878')
MCF7 = get_complete_results('MCF-7')

all_res = rbind(HCT116, GM12878, MCF7)




all_res %>% filter(gene_state %in% c('state1', 'state2')) %>% select(gene_name) %>% unique() #214
all_res %>% filter(gene_state %in% c('state1', 'state2')) %>% filter(validation_status) %>% select(gene_name) %>% unique() #170, after effect size cutoff: 159


all_res %>% filter(gene_state %in% c('state1', 'state2')) %>% select(id) %>% unique() #1195
all_res %>% filter(gene_state %in% c('state1', 'state2')) %>% filter(validation_status) %>% select(id) %>% unique() #591, after effect size cutoff: 525


## 
counts_per_status = all_res %>% select(id, cell_line, gene_name, gene_state, validation_status) %>% unique() %>% group_by(gene_name, gene_state, validation_status, cell_line) %>% count() %>% rename(snvs_per_status = n)
total_counts = all_res %>% select(id, cell_line, gene_name, gene_state, validation_status) %>% unique() %>% group_by(gene_name, gene_state, cell_line) %>% count() %>% rename(total_snvs = n)


prop_validated = left_join(counts_per_status, total_counts, by = c('gene_name', 'gene_state', 'cell_line')) %>% filter(validation_status == TRUE) %>% mutate(prop = snvs_per_status/total_snvs) %>% filter(gene_state %in% c('state1', 'state2')) 
prop_validated %>% filter(prop >= 0.5)


tmp_dat = left_join(counts_per_status, total_counts, by = c('gene_name', 'gene_state', 'cell_line')) %>% mutate(prop = snvs_per_status/total_snvs) %>% filter(gene_state %in% c('state1', 'state2')) 

make_snv_summary_dat = function(CELL_LINE, keep_conditions){
  in_data =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV_ALL_STATE_GENES.RData')
  load(in_data)
  total_snvs = all_results %>% filter(general_label %in% keep_conditions) %>% select(id, gene_state, validation_status) %>% unique() %>% group_by(gene_state) %>% count() %>% rename(total_snvs_per_state = n) # total snvs testable in both ActD and bru data
  overall_validation_counts = all_results %>% filter(general_label %in% keep_conditions) %>% select(id, gene_state, validation_status) %>% unique() %>% group_by(gene_state, validation_status) %>% count() %>% rename(total_snvs_validated = n)
  overall_validation_counts = left_join(total_snvs, overall_validation_counts, by = 'gene_state') %>% mutate(prop_validated = total_snvs_validated/total_snvs_per_state) %>% filter(validation_status == TRUE) %>% as.data.frame()
  overall_validation_counts$general_gene_state = ""
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state0')] = 'nonASE'
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state1', 'state2')] = 'asRS'
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state3')] = 'asRT'
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state4', 'state5', 'state6')] = 'mixed'
  summarized_overall_validation_counts = overall_validation_counts %>% group_by(general_gene_state) %>% summarize(total_snvs_per_state = sum(total_snvs_per_state), total_snvs_validated = sum(total_snvs_validated)) %>% mutate(prop_validated = total_snvs_validated/total_snvs_per_state)
  summarized_overall_validation_counts$cell_line = CELL_LINE
  overall_validation_counts$cell_line = CELL_LINE
  out = list(overall_validation_counts, summarized_overall_validation_counts)
  return(out)
}


make_gene_summary_dat = function(CELL_LINE, keep_conditions){
  in_data =  paste0('../results/', CELL_LINE, '_validate_by_BEAPR_ASE_SNV_ALL_STATE_GENES.RData')
  load(in_data)
  total_genes = all_results %>% filter(general_label %in% keep_conditions) %>% select(gene_name, gene_state) %>% unique() %>% group_by(gene_state) %>% count() %>% rename(total_genes_per_state = n) # total genes testable in both ActD and bru data
  
  
  overall_validation_counts = all_results %>% filter(general_label %in% keep_conditions) %>% filter(validation_status) %>% select(gene_name, gene_state) %>% unique() %>% group_by(gene_state) %>% count() %>% rename(total_genes_validated = n)
  overall_validation_counts = left_join(total_genes, overall_validation_counts, by = 'gene_state') %>% mutate(prop_validated = total_genes_validated/total_genes_per_state) %>% as.data.frame()
  overall_validation_counts$general_gene_state = ""
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state0')] = 'nonASE'
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state1', 'state2')] = 'asRS'
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state3')] = 'asRT'
  overall_validation_counts$general_gene_state[overall_validation_counts$gene_state %in% c('state4', 'state5', 'state6')] = 'mixed'
  summarized_overall_validation_counts = overall_validation_counts %>% group_by(general_gene_state) %>% summarize(total_genes_per_state = sum(total_genes_per_state), total_genes_validated = sum(total_genes_validated)) %>% mutate(prop_validated = total_genes_validated/total_genes_per_state)
  summarized_overall_validation_counts$cell_line = CELL_LINE
  overall_validation_counts$cell_line = CELL_LINE
  out = list(overall_validation_counts, summarized_overall_validation_counts)
  return(out)
}

both_val_methods_snv = rbind(make_snv_summary_dat('GM12878', c('nonASE_0h', 'ASE_0h', 'not validated'))[[2]], 
                             make_snv_summary_dat('HCT116', c('nonASE_0h', 'ASE_0h', 'not validated'))[[2]],
                             make_snv_summary_dat('MCF-7', c('nonASE_0h', 'ASE_0h', 'not validated'))[[2]])

both_val_methods_gene = rbind(make_gene_summary_dat('GM12878', c('nonASE_0h', 'ASE_0h', 'not validated'))[[2]], 
                             make_gene_summary_dat('HCT116', c('nonASE_0h', 'ASE_0h', 'not validated'))[[2]],
                             make_gene_summary_dat('MCF-7', c('nonASE_0h', 'ASE_0h', 'not validated'))[[2]])



nonASE_0h_val_methods_snv = rbind(make_snv_summary_dat('GM12878', c('nonASE_0h', 'not validated'))[[2]], 
                             make_snv_summary_dat('HCT116', c('nonASE_0h', 'not validated'))[[2]],
                             make_snv_summary_dat('MCF-7', c('nonASE_0h',  'not validated'))[[2]])

nonASE_0h_val_methods_gene = rbind(make_gene_summary_dat('GM12878', c('nonASE_0h',  'not validated'))[[2]], 
                              make_gene_summary_dat('HCT116', c('nonASE_0h',  'not validated'))[[2]],
                              make_gene_summary_dat('MCF-7', c('nonASE_0h',  'not validated'))[[2]])



ASE_0h_val_methods_snv = rbind(make_snv_summary_dat('GM12878', c('ASE_0h',  'not validated'))[[2]], 
                                make_snv_summary_dat('HCT116', c('ASE_0h',  'not validated'))[[2]],
                                make_snv_summary_dat('MCF-7', c('ASE_0h',  'not validated'))[[2]])

ASE_0h_val_methods_gene = rbind(make_gene_summary_dat('GM12878', c('ASE_0h',  'not validated'))[[2]], 
                                   make_gene_summary_dat('HCT116', c('ASE_0h',  'not validated'))[[2]],
                                   make_gene_summary_dat('MCF-7', c('ASE_0h',  'not validated'))[[2]])

pal = c('asRS' = '#875f9a' , 'asRT' = '#48d1cc', 'mixed' = '#8DA0CB', 'nonASE' = 'grey')



pdf('../figures/barplot_prop_validated_genes_all_states.pdf', width = 9, height = 5)
ggplot(both_val_methods_snv, aes(x = cell_line, y = prop_validated, label = total_snvs_validated, fill = general_gene_state)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('Both validation methods (SNV level results)')
ggplot(both_val_methods_gene, aes(x = cell_line, y = prop_validated, label = total_genes_validated, fill = general_gene_state)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('Both validation methods (gene level results)')


ggplot(both_val_methods_snv %>% filter(general_gene_state %in% c('nonASE', 'asRS')), aes(x = cell_line, y = prop_validated, label = total_snvs_validated, fill = general_gene_state)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('Both validation methods (SNV level results)')
ggplot(both_val_methods_gene %>% filter(general_gene_state %in% c('nonASE', 'asRS')), aes(x = cell_line, y = prop_validated, label = total_genes_validated, fill = general_gene_state)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('Both validation methods (gene level results)')


ggplot(nonASE_0h_val_methods_snv, aes(x = cell_line, y = prop_validated, fill = general_gene_state, label = total_snvs_validated)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('nonASE 0h validation (SNV level results)')
ggplot(nonASE_0h_val_methods_gene, aes(x = cell_line, y = prop_validated, fill = general_gene_state, label = total_genes_validated)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('nonASE 0h validated (gene level results')

ggplot(ASE_0h_val_methods_snv, aes(x = cell_line, y = prop_validated, fill = general_gene_state, label = total_snvs_validated)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('ASE 0h validation (SNV level results)')
ggplot(ASE_0h_val_methods_gene, aes(x = cell_line, y = prop_validated, fill = general_gene_state, label = total_genes_validated)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_bw() + scale_fill_manual(values = pal) + geom_label(position = position_dodge(width=1)) + ggtitle('ASE 0h validated (gene level results)')

dev.off()

