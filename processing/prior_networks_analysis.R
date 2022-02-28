rm(list=ls())

setwd("./secondtest/TCGA-BRCA/")

source("./BRCA_data_paired_wm/master_script_drugability_map.R")


load("workspace_drugmap_TCGA-BRCA.RData")

drug_set <- as.data.frame(t(combn(colnames(mat_dist_mat), 2)))
self <- data.frame(unique(drug_set[,1]), unique(drug_set[,1]))
colnames(self) <- colnames(drug_set)
drug_set <- rbind(drug_set, self)


cov_mat <- matrix(NA, nrow = length(colnames(mat_dist_mat)), ncol = length(colnames(mat_dist_mat)), dimnames = list(colnames(mat_dist_mat), colnames(mat_dist_mat)))
for(i in 1:nrow(drug_set)){
  cov_mat[drug_set[i,1], drug_set[i,2]] <- get_drug_input_coverage(Graph = conGraph, drugs_input = drug_set[i,], drug_target_df = filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% tolower(colnames(SMILES_dist2))),], drugs_col = 6, targets_col = 1)
}

which(cov_mat==min(cov_mat, na.rm = TRUE), arr.ind = TRUE)

cov_all <- get_drug_input_coverage(Graph = conGraph, drugs_input = colnames(mat_dist_mat), drug_target_df = filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% tolower(colnames(SMILES_dist2))),], drugs_col = 6, targets_col = 1)
tmp <- reshape::melt(cov_mat)

saveRDS(cov_mat, "coverage_mat_test.rds")

#acamprosate <- filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name=="acamprosate"),]

######

library(GA)
library(igraph)
library(tidyverse)

vec <- c()
df <- data.frame()

for(i in 1:length(colnames(mat_dist_mat))){
  prova <- get_drug_input_coverage(Graph = conGraph, drugs_input = colnames(mat_dist_mat)[i], drug_target_df = filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% tolower(colnames(SMILES_dist2))),], drugs_col = 6, targets_col = 1)
  vec <- c(prova[[1]], prova[[2]])
  df <- rbind(df, vec)
}

colnames(df) <- c("zscore", "pval")
rownames(df) <- colnames(mat_dist_mat)


# mim <- build.mim(t(inputinform$inform_mat), estimator = "pearson")
# adj_mat <- minet::aracne(mim)
# conGraph_aracne <- get_iGraph(adj_mat)


net_PPI <- read.csv("./drugmap/LUSC/functional_interaction_regulation_network_binary_only3.csv", header = TRUE, sep = "\t", quote = "")
rownames(net_PPI) <- net_PPI$X
net_PPI$X <- NULL
rowSums(net_PPI)

# conGraph_edges <- strsplit(as_ids(E(conGraph)), "|", fixed=TRUE)
# conGraph_edges <- do.call("rbind", conGraph_edges)
# 
# conGraph_edges <- data.frame(conGraph_edges)

net_PPI <- igraph::graph_from_adjacency_matrix(as.matrix(net_PPI), weighted = TRUE, mode = "undirected")

######
# net_PPI <- igraph::read_graph("./Parsimonious Composite Network (PCNet).cx.graphml", format = "graphml")
# net_PPI <- as.undirected(net_PPI, mode = "collapse")
###

#net_PPI <- read.table("./Parsimonious Composite Network (PCNet).cx", header = TRUE, sep = "\t", quote = "")
length(E(igraph::intersection(net_PPI, conGraph)))
length(intersect(as_ids(V(conGraph)), as_ids(V(net_PPI))))

prior_net <- igraph::intersection(conGraph, net_PPI)
non_prior_net <- igraph::difference(conGraph, net_PPI)
length(which(igraph::degree(prior_net)==0))
prior_net <-  delete.vertices(prior_net, which(igraph::degree(prior_net)==0))
non_prior_net <-  delete.vertices(non_prior_net, which(igraph::degree(non_prior_net)==0))

length(unique(filter_by_L1000_drugs$dat.drug.molecule_name))
length(unique(filter_by_L1000_drugs$dat.target.gene_info.symbol))


filter_drugs_on_prior <- filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.target.gene_info.symbol %in% as_ids(V(prior_net))),]
filter_drugs_on_non_prior <- filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.target.gene_info.symbol %in% as_ids(V(non_prior_net))),]

length(unique(filter_drugs_on_prior$dat.drug.molecule_name))
write.table(filter_drugs_on_prior, file = "filter_by_L1000_drugs_on_prior.txt", row.names=FALSE, sep = "\t", quote = FALSE)
length(unique(filter_drugs_on_non_prior$dat.drug.molecule_name))
length(unique(filter_drugs_on_prior$dat.target.gene_info.symbol))
length(unique(filter_drugs_on_non_prior$dat.target.gene_info.symbol))
intersect(filter_drugs_on_prior$dat.drug.molecule_name, colnames(mat_dist_mat))

summary(graph_rank[names(graph_rank) %in% as_ids(V(prior_net))]) ### median rank on all the nodes
summary(graph_rank[names(graph_rank) %in% as_ids(V(non_prior_net))])
summary(graph_rank[names(graph_rank) %in% filter_drugs_on_prior$dat.target.gene_info.symbol]) ### median rank on drug targets
summary(graph_rank[names(graph_rank) %in% filter_drugs_on_non_prior$dat.target.gene_info.symbol])

write.graph(prior_net, file = "prior_network_all_3_layers_alisa_BRCA.gml", format = "gml")


