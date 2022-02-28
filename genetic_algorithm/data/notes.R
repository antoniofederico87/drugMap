# graph_rank = readRDS("data/TCGA-BRCA/graph_rank_BRCA.rds")

# order: Integer giving the order of the neighborhood. = 1
compute_cov_score <- function(targets, Graph, order, graph_rank) {

  target_neighbour_list <- list()
  
  for(target in 1:length(targets)){
    # per ogni target, trova tutti i vicini a distanza 1 nel prior graph
    target_neighbour_list[[target]] <- igraph::ego(Graph, order = order, nodes = targets[target])
  }
  # estrai i nomi dei nodi di ogni vicinato
  target_neighbour_list <- lapply(target_neighbour_list, function(x) igraph::as_ids(x[[1]]))
  # ogni vicino di ogni target viene considerato una sola volta
  covered_nodes_overall <- unique(unlist(target_neighbour_list))
  # calcolo il rank mediano ti tutti i vicini considerati according to il ranking nel file graph_rank
  inform_rank <- median(which(names(graph_rank) %in% covered_nodes_overall))
  # calcola la proporzione di coverage della network di prior knowledge
  coverage_overall <- length(covered_nodes_overall)/length(V(Graph))
  # drug score???
  toninos_drug_score <- coverage_overall/inform_rank #chiedi a Giovanni se va bene, altrimenti usiamo solo il coverage
  return(toninos_drug_score)

}


get_drug_input_coverage <- function(Graph, drugs_input, drug_target_df, drugs_col, targets_col, graph_rank, deg_dist, order = 1) {
  toninos_drug_score <- c()

  # estrai i target delle drug selezionate dall'individuo corrente del GA
  targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% drugs_input)]))
  
  # cerca solo i target che sono annotati nel grafo di prior knowledge
  targets <- targets[targets %in% as_ids(V(Graph))]
  
  #!!!
  observed <- compute_cov_score(targets, Graph, order, graph_rank)

  # bootstrap per calcolare il pvalue
  simulated <- rep(0,100)
  for(i in 1:100){
    random_targ <- random_targets(targets, Graph, deg_dist)
    simulated[i] <- compute_cov_score(random_targ, Graph, order, graph_rank)
  }
  z_score = (observed - mean(simulated))/sd(simulated)
  list(z_score = z_score, pvalue = pnorm(-abs(z_score)))
}



coverage_sum <- function(filter_by_L1000_drugs,conGraph,SMILES_dist2, graph_rank, drugs_input) {
	# prendi da L1000 i target del set di 316 drug che stiamo ottimizzando
	drug_target_df = filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% tolower(colnames(SMILES_dist2))),]
	
	# costruisco un DF con tutti i nodi e i loro degree in ordine crescente
	deg_dist <- tibble(vertex=as_ids(V(conGraph)), dist = degree(conGraph)) %>% arrange(dist)
	# quantizzo la distribuzione dei digree in bins di almeno 200 elementi a bin
	binned <- binr::bins(deg_dist$dist, target.bins = 20,minpts = 100)
	
	# aggiunge una colonna al DF dove dice il degree a che bin appartiene e ordina per i bin in ordine crescente
	deg_dist <- deg_dist %>% dplyr::mutate(BIN = findInterval(dist,binr::bins.getvals(binned))) %>% arrange(BIN)
	# raggruppi tutti i nodi che appartengono allo stesso bin e ti registri delle statistiche (numero, minimo, massimo) in un nuovo DF
	deg_dist %>% group_by(BIN) %>% summarise(n=n(),min=min(dist), max=max(dist)) %>% as.data.frame()
	
	
	get_drug_input_coverage(Graph = conGraph, drugs_input = drugs_input,
	drug_target_df = drug_target_df, drugs_col = 6,
	targets_col = 1, graph_rank=graph_rank, deg_dist = deg_dist)

}

