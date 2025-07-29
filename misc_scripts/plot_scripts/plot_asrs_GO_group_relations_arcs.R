rm(list = ls())
library(tidyverse)
library(viridis)
library(patchwork)
#library(hrbrthemes)
library(igraph)
library(ggraph)
#library(colormap)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ComplexUpset)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(R.utils)

## code for plotting adapted from https://www.data-to-viz.com/story/AdjacencyMatrix.html

#setwd('/Users/elainehuang/Desktop/project-gxxiao4/ASE_ASRS/REVISION1/HIER_MODEL_V2/scripts')
#load('../GO_results/rrvgo_results_min_occur10.Rdata')


#load('../results/all_asrs_ensembl_genes.txt.GO_table.RData')

MIN_OCCUR = 5
SIM_THRESHOLD = 0.7
out_fn = paste0('../GO_RESULTS/REMOVE_5_CELL_LINES_rrvgo_results_min_occur', MIN_OCCUR, '_SIM_THRESHOLD_', SIM_THRESHOLD, '_GENE_LENGTH_RPKM_MATCHED_BKG.Rdata')

asrs_GO_table = '../GO_RESULTS/REMOVE_5_CELL_LINES_ALL_state1_state2_bg_GENE_LENGTH_RPKM_CONTROLLED_THROW_OUT_GENES_yes_query_GO_table_gene_length_RPKM_controlled.RData'

load(out_fn)
load(asrs_GO_table)

asrs_reducedTerms = left_join(asrs_reducedTerms_size %>% dplyr::select(-score), asrs_reducedTerms %>% dplyr::select(term, score), by = 'term')
INTEREST_GO_TERMS =  asrs_reducedTerms %>% dplyr::filter(parentTerm %in% c('innate immune response', 'response to Gram-positive bacterium', 'positive regulation of autophagy', 'proteolysis')) %>% dplyr::select(term) %>% unlist(use.names = FALSE)



MIN_TERMS_PER_CLUSTER = 2
bigClusters = table(asrs_reducedTerms$parent) %>% as.data.frame() %>% filter(Freq >= MIN_TERMS_PER_CLUSTER)
combine_dat = left_join(asrs_reducedTerms %>% filter(parent %in% bigClusters$Var1), query_GO_table %>% dplyr::rename(go = go_id), by = 'go', multiple = 'all')
#inner_join(combine_dat, traits_and_genes %>% dplyr::rename(hgnc_symbol = gene_name), by = 'hgnc_symbol', relationship = 'many-to-many') %>% dplyr::select(trait, hgnc_symbol, parentTerm) %>% unique() %>% group_by(trait, hgnc_symbol) %>%  summarise(parentTerm_collection = paste(parentTerm, collapse=', ')) %>% print(n = 58)

gene_list = list()
parents = unique(combine_dat$parentTerm)
for (i in 1:length(parents)){
  asrs_genes = combine_dat %>% filter(parentTerm == parents[i])
  asrs_genes = unique(asrs_genes$hgnc_symbol)
  asrs_genes = asrs_genes[asrs_genes != ""] # need to remove blank or the list_to_matrix function will not work
  gene_list[[i]] = asrs_genes
  names(gene_list)[i] = parents[i]
  #names(gene_list)[i] = paste0('group', i)
}




#names(gene_list) = paste0('group', c(1:length(gene_list)))
#ASRS_gene_mat = list_to_matrix(gene_list, universal_set = as.character(unlist(gene_list)))
ASRS_gene_mat = list_to_matrix(gene_list)
ASRS_tbl = as_tibble(ASRS_gene_mat)
ASRS_df = data.frame(ASRS_gene_mat)

interest_clusters = c('innate.immune.response', 'proteolysis', 'defense.response.to.Gram.positive.bacterium', 'positive.regulation.of.autophagy')
ASRS_df = ASRS_df %>% dplyr::select(interest_clusters) 


## need to make gene-to-gene object
from_vec = list()
for (i in 1:nrow(ASRS_df)){
  hit_rows = ASRS_df[i,]
  hit_terms = names(ASRS_df[i,])[ASRS_df[i,] == 1]
  if (length(hit_terms) <= 1){
    print('no intersections')
  }
  else {
    print(hit_terms)
    gene_name = row.names(hit_rows)
    add_rows = t(combn(paste0(gene_name, '.', hit_terms), 2))
    from_vec[[i]] = add_rows
  }
}

all_rows = do.call("rbind", from_vec)
all_rows = data.frame(all_rows)

all_rows = all_rows %>% tidyr::separate(X1, into = c('gene_name', 'group'), remove = FALSE, sep = '\\.', extra = 'merge') %>% arrange(group) %>% dplyr::select(-c(gene_name, group))


all_rows = dplyr::rename(all_rows, from = X1, to = X2)
all_rows$value = 1
connect = all_rows
# Number of connection per person
c( as.character(connect$from), as.character(connect$to)) %>%
  as.tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> coauth

colnames(coauth) <- c("name", "n")
#dim(coauth)
# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
# Find community
com <- walktrap.community(mygraph)
#max(com$membership)
#Reorder dataset and make the graph
coauth <- coauth %>% 
  mutate( grp = com$membership) %>%
  arrange(grp) %>%
  mutate(name=factor(name, name))

coauth = coauth %>% separate(name, into = c('gene_name', 'func_group'), sep = '\\.', remove = FALSE, extra = 'merge') %>% arrange(func_group, desc(n))
# keep only 10 first communities
#coauth <- coauth %>% 
#  filter(grp<16)
# keep only this people in edges
connect <- connect %>%
  filter(from %in% coauth$name) %>%
  filter(to %in% coauth$name)


# Add label angle
number_of_bar=nrow(coauth)
coauth$id = seq(1, nrow(coauth))
angle= 360 * (coauth$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
coauth$hjust <- ifelse(angle > 90 & angle<270, 1, 0)
coauth$angle <- ifelse(angle > 90 & angle<270, angle+180, angle)


# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
# prepare a vector of n color in the viridis scale
#mycolor <- colormap(colormap=colormaps$viridis, nshades=max(coauth$grp))
#mycolor <- sample(mycolor, length(mycolor))
pal = brewer.pal(8, "Dark2")
mycolor = c(pal[1], pal[3], pal[5], pal[6])



# Make the graph
p1 = ggraph(mygraph, layout="linear") + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.5, fold=TRUE) +
  geom_node_point(aes(size=n, color=as.factor(func_group), fill=func_group), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=gene_name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  theme_void() +
  theme(plot.margin=unit(c(0,0,0.4,0), "null"),
    panel.spacing=unit(c(0,0,3.4,0), "null")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 


p2 = ggraph(mygraph, layout="linear") + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.5, fold=TRUE) +
  geom_node_point(aes(size=n, color=as.factor(func_group), fill=func_group), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=gene_name), hjust=1, nudge_y = -1.1, size=2.3) +
  theme_void()  + coord_flip()  



p3 = ggraph(mygraph, layout="circle") + 
  geom_edge_link(edge_colour="black", edge_alpha=0.2, edge_width=0.5) +
  geom_node_point(aes(size=n, color=as.factor(func_group), fill=func_group), alpha=0.9) +
  scale_size_continuous(range=c(3,8)) +
  scale_color_manual(name = 'GO group', values=mycolor, labels = c('Defense response to Gram-positve bacterium', 'Innate immune response', 'Positive regulation of autophagy', 'Proteolysis')) +
  geom_node_text(aes(label=paste("    ",gene_name,"    "), angle=angle, hjust=hjust), size=2.3, color="black") +
  theme_void() +
  theme(
   # legend.position="none",
    plot.margin=unit(c(0,0,0,0), "null"),
    panel.spacing=unit(c(0,0,0,0), "null")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +  guides(size=FALSE, fill = FALSE, color = guide_legend(override.aes = list(size = 5))) 

#pdf('../figures/asrs_GO_result_arcdiagram_vertical.pdf', height = 13, width = 4.5)
#p2 + theme(legend.position = 'none')
#dev.off()

write.table(all_rows, file = '../../../SOURCE_DATA/MAIN_FIGS/fig5c.txt', sep = '\t', row.names = F, quote = F)

pdf('../figures/asrs_GO_result_arcdiagram_horizontal.pdf', width = 13, height = 13)
p3 + theme(legend.position="bottom",
           legend.box.spacing = unit(1, "pt"))
dev.off()


pdf('../MANUSCRIPT_FIGURES/fig5/fig5c.pdf', width = 4.5, height = 4.5)
p3 + theme(legend.position="none",
           legend.box.spacing = unit(1, "pt"))
dev.off()


