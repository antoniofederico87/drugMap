### This is the starting point of the package ###
#group= as.factor(c(rep("T", length(tumor)), rep("N", length(normal))))

#' get_informative_genes is a function that allows to identify the informative genes (e.g. differentially expressed genes) on which a co-expressed network will be built.
#' 
#' @param expr_mat An expression matrix where rows are genes and columns are samples.
#' @param group A factor indicating the  biological conditions of the samples (e.g. group=as.factor(c(rep("T", length(tumor)), rep("N", length(normal))))).
#' @param test The statistical test to be performed for the identification of the informative genes. If test="wilcoxon" a Wilcoxon test is performed. If test="var.test" an F-test on the variance is performed. The default is test="wilcoxon + var.test" and both Wilcoxon and var.test are performed.
#' @param percentile An integer indicating the percentile of the ranked distribution of informative genes to be returned. The default is percentile=10.
#' @return An expression matrix for the informative genes to be used as input of the get_coexpression_network function.
#' @examples
#' get_informative_genes(expr_mat = expr_matrix, group = as.factor(c(rep("T", length(tumor)), rep("N", length(normal)))), test = "wilcoxon + var.test", percentile = 10).


get_informative_genes <- function(expr_mat, group, test="wilcoxon + var.test", percentile=10) {
  
  if(is.null(expr_mat)){
    
    stop("Error: please provide an expression matrix!")
    
  }
  
  if(is.null(group)){
    stop("Error: please provide groups!")
    
  }
  
  if(length(group)!=ncol(expr_mat)){
    stop("Error: the group length is different from the sample size!")
  }
  
  if(!isClass(percentile, Class = "numeric")){
    stop("Error:please provide an integer")
  }
  
  if(percentile < 1  | percentile >100) {
    
    stop("Error: please provide a value between 1 and 100!")
    
  }
    
  
  if(!test %in% c("wilcoxon","var.test","wilcoxon + var.test")){
    stop("Error: please provide one of the following tests: wilcoxon, var.test, wilcoxon + var.test")
  }
  
  
  if (test == "wilcoxon") {
    
    wilres <- apply(expr_mat, 1, function(x){wilcox.test(x ~ group, data = expr_mat)})
    pvalwil <- sapply(wilres, function(x){x$p.value})
    pvalwil[pvalwil==0] <- 0.000000000000001
    
    logFC <- log(rowMeans(expr_mat[, which(group==names(table(group)[2]))])/rowMeans(expr_mat[, which(group==names(table(group)[1]))])) ### Change this according to factors levels order
    
    wilcox_score <- abs(logFC*(-log(pvalwil)))
    wilcox_score <- stats::setNames(object = wilcox_score, nm = rownames(expr_mat))
    
    uqua <- length(expr_mat[,1])*percentile/100
    
    inputinform <- list(expr_mat[names(sort(wilcox_score, decreasing=TRUE))[1:uqua],], wilcox_score)
    names(inputinform) <- c("inform_mat", "wilcox_score")
    
  }else if(test == "var.test") {
    
    varres <- apply(expr_mat, 1, function(x){var.test(x ~ group, data = expr_mat)})
    pvalvar <- sapply(varres, function(x){x$p.value})
    pvalvar[pvalvar == 0] <- 0.000000000000001
    
    logFC <- log(rowMeans(expr_mat[, which(group==names(table(group)[2]))])/rowMeans(expr_mat[, which(group==names(table(group)[1]))])) ### Change this according to factors levels order
    
    var_score <- abs(logFC*(-log(pvalvar)))
    
    var_score <- stats::setNames(object = var_score, nm = rownames(expr_mat))
    
    uqua <- length(expr_mat[,1])*percentile/100
    
    
    inputinform <- list(expr_mat[names(sort(var_score, decreasing=TRUE))[1:uqua],], var_score)
    names(inputinform) <- c("inform_mat", "var_score")
    
    
  }else if(test == "wilcoxon + var.test"){
    
    wilres <- apply(expr_mat, 1, function(x){wilcox.test(x ~ group, data = expr_mat)})
    pvalwil <- sapply(wilres, function(x){x$p.value})
    pvalwil[pvalwil==0] <- 0.000000000000001
    
    
    varres <- apply(expr_mat, 1, function(x){var.test(x ~ group, data = expr_mat)})
    pvalvar <- sapply(varres, function(x){x$p.value})
    pvalvar[pvalvar == 0] <- 0.000000000000001
    
    
    logFC <- log(rowMeans(expr_mat[, which(group==names(table(group)[2]))])/rowMeans(expr_mat[, which(group==names(table(group)[1]))])) ### Change this according to factors levels order
    
    wilcox_score <- abs(logFC*(-log(pvalwil)))
    
    wilcox_score <- stats::setNames(object = wilcox_score, nm = rownames(expr_mat))
    
    var_score <- abs(logFC*(-log(pvalvar)))
    var_score <- stats::setNames(object = var_score, nm = rownames(expr_mat))
    
    uqua <- length(expr_mat[,1])*percentile/100
    
    #all_samples_filtered_wilcoxon <- expr_mat[names(sort(wilcox_score, decreasing=TRUE))[1:uqua],]
    #all_samples_filtered_var <- expr_mat[names(sort(var_score, decreasing=TRUE))[1:uqua],]
    
    inputinform <- list(expr_mat[union(names(sort(wilcox_score, decreasing=TRUE))[1:uqua], names(sort(var_score, decreasing=TRUE))[1:uqua]),], var_score, wilcox_score)
    names(inputinform) <- c("inform_mat", "wilcox_score", "var_score")
    
  }
  return(inputinform)
}



#' get_informative_genes is a function that allows to identify the informative genes (e.g. differentially expressed genes) on which a co-expressed network will be built.
#' 
#' @param Graph An igraph object.
#' @param ranked_vectors One or more ranked vectors to be considered in the ranking of the network (e.g. scores of biological significance), together with the centrality scores.
#' @param test The statistical test to be performed for the identification of the informative genes. If test="wilcoxon" a Wilcoxon test is performed. If test="var.test" an F-test on the variance is performed. The default is test="wilcoxon + var.test" and both Wilcoxon and var.test are performed.
#' @return A named vector representing the rank of the graph's nodes.
#' @examples
#' graph_rank <- get_graph_named_ranklist(Graph = conGraph, ranked_vectors = list(var_score, wilcox_score))


get_graph_named_ranklist <- function(Graph, ranked_vectors){
  
  if(!class(Graph)=="igraph") {
    stop("Error: please provide an object of class igraph!")
  }
  
  if(is.null(Graph)) {
    stop("Error: please provide an igraph object!")
  }
  
  if(is.null(ranked_vectors)) {
    stop("Error: please provide a list of one or more ranked vectors!")
  }
  
  if(!class(ranked_vectors)=="list") {
    stop("Error: please provide a list of ranked vectors!")
  }
  
  
  if(length(ranked_vectors)==1) {
    Graph=set.vertex.attribute(Graph, "ranked_vectors_1", value = ranked_vectors[[1]][V(Graph)$name])
    
  }else if(length(ranked_vectors)>1) {
    for (rank in 1:length(ranked_vectors)) {
      Graph=set.vertex.attribute(Graph, paste0("ranked_vectors_", rank), value = ranked_vectors[[rank]][V(Graph)$name])
      
    }
    
  }
  
  net_attr <- c("betweenness", "cc", "degree", "closeness", "eigenvector", paste0("ranked_vectors_", 1:length(ranked_vectors)))
  
  ranked.graph=get_ranked_gene_list(iGraph = Graph, rank_list_attr = net_attr, debug_output = TRUE)
  named_ranklist=setNames(c(1:length(ranked.graph)), ranked.graph)
  
  return(named_ranklist)
  
}

### In the following steps, we retrieve the OpenTargets data about drugs, targets and diseases, and we filter them for the drugs available in L1000
### and then we map such drugs and diseases on the coexpression network, on the significant pvalues obtained by permutation test performed on the jaccard indices
### on the sets of I and II degree of the drug trget nodes


#' Convenience function to parse OpenTargets JSON files
#' 
#' @param json_file The OpenTargets JSON file.
#' @return A dataframe representing the OpenTargets database.
#' @examples
#' opentargets_parsed <- parse_opentargets_json(json_file≈"")


parse_opentargets_json <- function(json_file) {
  
  if(is.null(json_file)) {
    stop("Error: please type in the OpenTargets JSON file to be processed!")
  }
 
  #json_file <- "19.02_evidence_data.json"
  out <- lapply(readLines(json_file), fromJSON) # load opentargets json file
  dat <- rbindlist(lapply(out, function(x) {as.list(unlist(x))}), fill=TRUE)
  final_opentargets_parsed <- data.frame(dat$target.gene_info.symbol, dat$target.gene_info.name, dat$target.gene_info.geneid, dat$disease.efo_info.label, dat$disease.name, dat$drug.molecule_name, dat$drug.id, dat$drug.molecule_type, dat$drug.max_phase_for_all_diseases.label, dat$unique_association_fields.pathway_id)
  
  return(final_opentargets_parsed)
}


#final_opentargets_parsed <- readRDS("/home/antonio/final_opentargets_parsed.rds")

#' Once obtained a gene network, the function maps the OpenTargets drugs on the target nodes.
#' 
#' @param Graph An igraph object.
#' @param parsed_opentargets A dataframe representing the OpenTargets database (e.g. obtained from the parse_opentargets_json).
#' @return A dataframe reporting the drugs and diseases mapped on the provided graph.
#' @examples
#' map_opentargets <- map_opentargets_on_graph(Graph = conGraph, parsed_opentargets = opentargets_parsed)
 

map_opentargets_on_graph <- function(Graph, parsed_opentargets) {
  
  if(is.null(Graph)) {
    stop("Error: please provide a valid igraph object!")
  }
  
  if(!class(Graph)=="igraph") {
    stop("Error: the Graph object is not of class igraph!")
  }
  
  if(is.null(parsed_opentargets)) {
    stop("Error: please provide the parsed OpenTargets database (e.g. obtained from the parse_opentargets_json function)")
  }
  
  
  if(!class(parsed_opentargets)=="data.frame") {
    stop("Error: please provide an object of class data.frame!")
  }
  
  map_opentargets <- parsed_opentargets[which(parsed_opentargets$dat.target.gene_info.symbol %in% names(V(Graph))),]
  map_opentargets$dat.drug.molecule_name <- tolower(map_opentargets$dat.drug.molecule_name)
  map_opentargets$dat.drug.molecule_name <- gsub(map_opentargets$dat.drug.molecule_name, pattern = " ", replacement = ".")
  map_opentargets$dat.disease.efo_info.label <- gsub(map_opentargets$dat.disease.efo_info.label, pattern = " ", replacement = ".")
  
  return(map_opentargets)
  
}


#' Compute a squared matrix reporting the Jaccard similarity indices among couples of drug targets neighborhoods.
#' 
#' @param Graph An igraph object.
#' @param target_nodes A vector of gene symbols likely to be drug targets.
#' @param neighborhood The order of neighborhood of the provided target nodes.
#' @return A squared matrix reporting the Jaccard similarity indices among couples of drug targets neighborhoods.
#' @examples
#' jacc_matrix <- get_jacc_drug_targets_neighbors(Graph = conGraph, target_nodes = c("GAPDH", "PPARG", "ELF1", "FOX1"), neighborhood=2)

get_jacc_drug_targets_neighbors <- function(Graph, target_nodes, neighborhood) {
  
  drug_target_mat_L1000 <- matrix("NA", ncol = length(target_nodes), nrow = length(target_nodes), dimnames = list(target_nodes, target_nodes))
  jacc_index <- c()
  
  
  for (i in 1:length(as.character(unique(target_nodes)))) {
    degree_i <- igraph::ego(conGraph, order = neighborhood, nodes = as.character(unique(target_nodes)[i]))
    for (k in 1:length(as.character(unique(target_nodes)))) {
      degree_k <- igraph::ego(conGraph, order = neighborhood, nodes = as.character(unique(target_nodes)[k]))
      jacc_index <- length(intersect(as.character(names(degree_i[[1]])), as.character(names(degree_k[[1]]))))/length(union(as.character(names(degree_i[[1]])), as.character(names(degree_k[[1]]))))
      drug_target_mat_L1000[i, k] <- jacc_index
    }
    
  }
  
  return(drug_target_mat_L1000)
  
}


#' Performs a permutation test to assess the significance of the jaccard scores among couples of drug targets (e.g. obtained from the get_jacc_drug_targets_neighbors function).
#' 
#' @param drug_target_mat A squared matrix containing similarity scores between drug targets pairs (e.g. obtained from the get_jacc_drug_targets_neighbors function).
#' @param Graph An igraph object.
#' @param nIter Number of permutations to be performed.
#' @return A squared matrix reporting the statistical significance of the similarity measures (e.g. Jaccard indices) among couples of drug targets neighborhoods.
#' @examples
#' pval_1000 <- test_similarity_significance(drug_target_mat = jacc_matrix, Graph = conGraph, nIter = 1000)

test_similarity_significance <- function(drug_target_mat, Graph, nIter=1000) {
  
  if(is.null(drug_target_mat)){
    stop("Error: please provide a squared drug target matrix!")
  }
  
  if(!class(drug_target_mat)=="matrix") {
    stop("Error: please provide an object of class matrix!")
  }
  
  if(is.null(Graph) | !class(Graph)=="igraph"){
    stop("Error: please provide an igraph object!")
  }
  
  if(!class(nIter)=="numeric"){
    stop("Error: please provide the number of iterations!")
  }
  
  dist <- c()
  jacc_original <- c()
  nIter = nIter
  PvalueMatL1000 = matrix(NA, nrow = length(rownames(drug_target_mat)), ncol = length(rownames(drug_target_mat)) )
  
  
  pb = txtProgressBar(min = 0, max = length(rownames(drug_target_mat)), style = 3)
  
  for (i in 1:length(rownames(drug_target_mat))) {
    for (k in 1:length(colnames(drug_target_mat))) {
      print(paste(i, k, sep=","))
      jacc_permutated <-  c()
      jacc_original <- as.numeric(drug_target_mat[i,k])
      
      # start_time <- Sys.time()
      
      pbin = txtProgressBar(min = 1, max=nIter, style=2)
      for(nperm in 1:nIter){
        dist <- igraph::distances(graph = Graph, v = colnames(drug_target_mat)[i], to = colnames(drug_target_mat)[k], mode = "all")
        
        pick_first_node <- base::sample(names(V(Graph))[-c(which(names(V(Graph))==colnames(drug_target_mat)[i]), which(names(V(Graph))==colnames(drug_target_mat)[k]))], size = 1)
        first_node <- igraph::ego(Graph, order = dist, nodes = pick_first_node)
        pick_second_node <- base::sample(as.character(names(first_node[[1]])), size = 1)
        second_node <- igraph::ego(Graph, order = dist, nodes = pick_second_node)
        
        new_jacc <- length(intersect(as.character(names(first_node[[1]])), as.character(names(second_node[[1]]))))/length(union(as.character(names(first_node[[1]])), as.character(names(second_node[[1]]))))
        
        jacc_permutated <-  c(jacc_permutated, new_jacc)
        setTxtProgressBar(pbin,nperm)
      }
      close(pbin)
      # end_time = Sys.time()
      # end_time - start_time
      
      pvalue <- 1-sum(jacc_original > jacc_permutated)/nIter
      PvalueMatL1000[i,k] = pvalue
    }
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(PvalueMatL1000)
  #saveRDS(PvalueMatL1000, file = "PvalueMatL1000.rds")
  
}


#' Produces expression matrices from L1000 data in tab delimited format.
#' 
#' @param L1000_data An object of class gct or gctx from the LINCS1000 study (e.g. parsed by the parse.gct function of the cmapR package).
#' @param L1000_metadata Object containing the metadata relative to the L1000 dataset (e.g file GSE70138_Broad_LINCS_sig_info_2017-03-06.txt).
#' @param include_controls A logical value indicating whether to include the plate-specific control samples in the expression matrix.
#' @param L1000_data_CTRL An object of class gct or gctx from the LINCS1000 study for control samples (e.g. DMSO).
#' @param L1000_metadata_CTRL Object containing the metadata relative to the L1000 dataset for control samples.
#' @param dest_path Folder to be created as a destination folder for the expression matrices. Default to current directory.
#' @param inputgenes A character vector of gene symbols to be included in the expression matrix (e.g. L1000 landmark genes). Default all the genes will be considered.
#' @return A squared matrix reporting the statistical significance of the similarity measures (e.g. Jaccard indices) among couples of drug targets neighborhoods.
#' @examples
#' pval_1000 <- test_similarity_significance(drug_target_mat = jacc_matrix, Graph = conGraph, nIter = 1000)

#landmark_genes <- read.delim("GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", header=TRUE, quote="")
#expr_matrices <- "/home/MOKA/antonio/BRCA_data_paired_wm/L1000_expression_matrices"
#L1000_data = entire_gctx, L1000_metadata = col_meta_filtered_all


build_drug_expression_matrices <- function(L1000_data, L1000_metadata, include_controls = FALSE, L1000_data_CTRL=NULL, L1000_metadata_CTRL=NULL, dest_path = ".", inputgenes) {
  
  if(is.null(L1000_data)) {
    stop("Error: please provide a L1000 expression matrix!")
  }
  
  if(is.null(L1000_metadata)) {
    stop("Error: please provide a L1000 metadata file!")
  }
  
  if(!class(include_controls)=="logical") {
    stop("Error: please provide a logical value to the include_controls argument!")
  }
  
  if(include_controls==TRUE & is.null(L1000_data_CTRL)) {
    stop("Error: please provide a L1000 expression matrix for control samples!")
  }
  
  if(include_controls==TRUE & is.null(L1000_metadata_CTRL)) {
    stop("Error: please provide a L1000 metadata file for control samples!")
  }
  
  if (!class(dest_path)=="character") {
    stop("Error: please provide a valid path!")
  }
  
  if (is.null(inputgenes)) {
    inputgenes <- rownames(L1000_data@mat)
  }
  
  # if (!class(inputgenes)=="character") {
  #   stop("Error: please provide a vector of genes symbols!")
  # }
  
  if (!file.exists(dest_path)){
    dir.create(file.path(dest_path))
  } 
  
  
  distil_ids <- c()
  distil_ids_CTRL <- c()
  plate <- c()
  CTRL_expr <- data.frame()
  finalmat <- data.frame()
  
  for (drug in 1:length(unique(L1000_metadata$pert_iname))) {
    distil_ids <- L1000_metadata[which(L1000_metadata$pert_iname==L1000_metadata$pert_iname[drug]), "distil_id"]
    if(length(distil_ids)>1){
    drugs_expr <- L1000_data@mat[which(rownames(L1000_data@mat) %in% inputgenes), distil_ids] # filtra landmark genes
    rownames(drugs_expr) <- inputgenes[which(inputgenes %in% rownames(L1000_data@mat))]
    }
    
      if(include_controls == TRUE) {
        plate <- sapply(strsplit(distil_ids, ":", fixed=TRUE), "[", 1)
        CTRL_expr_all_plates <- data.frame(row.names = 1:length(rownames(drugs_expr)))
        for (distil_id in 1:length(plate)) {
          distil_ids_CTRL <- L1000_metadata_CTRL[grep(plate[distil_id], L1000_metadata_CTRL$distil_id, fixed = TRUE), "distil_id"]
          CTRL_expr <- L1000_data_CTRL@mat[which(rownames(L1000_data_CTRL@mat) %in% inputgenes), distil_ids_CTRL] # filtra landmark genes
        
          CTRL_expr_all_plates <- cbind(CTRL_expr_all_plates, CTRL_expr)
      }
      
        rownames(CTRL_expr_all_plates) <- inputgenes[which(inputgenes %in% rownames(L1000_data_CTRL))]
        colnames(CTRL_expr_all_plates) <- sapply(colnames(CTRL_expr_all_plates), function(x) paste0("CTRL_", x))
      
        finalmat <- cbind(drugs_expr, CTRL_expr_all_plates)
        write.table(finalmat, file = paste0(dest_path, "/", "Expression_matrix_", unique(L1000_metadata$pert_iname)[drug], ".txt"), quote = FALSE, sep = "\t", row.names=TRUE, col=NA)
      
    }else{
      write.table(drugs_expr, file = paste0(dest_path, "/", "Expression_matrix_", unique(L1000_metadata$pert_iname)[drug], ".txt"), quote = FALSE, sep = "\t", row.names=TRUE, col=NA)
    }
  
}
  
  }




### Infer coexpression networks only with treatments samples, excluding the drugs for which
### less than 5 samples are available
### To be ran as batch job ###

#' Infer coexpression networks excluding the drugs for which less than minSamples samples are available
#' 
#' @param path The path where the expression matrices are stored.
#' @param pattern An optional regular expression. Only file names which match the regular expression will be returned.
#' @param minSamples Only the matrices with minSamples number of samples (columns) will be considered.
#' @param ncores Number of cores to infer the coexpression networks.
#' @return A list of adjacency matrices.
#' @examples
#' netlist <- get_multiple_coexpression_networks(path = "/home/user/expression_mat", pattern = "expression_mat*", minSamples = 5, ncores = 12)

get_multiple_coexpression_networks <- function(path, pattern=NULL, minSamples=5, ncores = 12) {
  
  if (is.null(path)) {
    path <- "."
  }
  
  myfiles <- list.files(path = path, pattern = pattern)
  netlist <- rep(list(NULL), length(myfiles))
  inputinform <- data.frame()
  
  for (exprmatrix in 1:length(myfiles)) {
    
  #for (exprmatrix in 1:5) {  
    print(paste("Processing: ", myfiles[exprmatrix]))
    print(paste("Drug", exprmatrix, "of", length(myfiles)))
    
    inputinform  <- read.delim(paste0(path, myfiles[exprmatrix]), header = TRUE, quote = "")
    
    if (is.null(minSamples)) {
      minSamples <- length(colnames(inputinform))
    }
    
    if (length(colnames(inputinform))>=minSamples) {    # Filter out the matrices having less than minSamples samples
      
      #controls <- inputinform[, grep(pattern = "CTRL", x = colnames(inputinform))]
      
      generatematrices=get_ranked_consensus_matrix(gx_table = inputinform, iMethods = c("clr","aracne","mrnet"),
                                                   iEst = c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
                                                   iDisc=c("none","equalfreq","equalwidth","globalequalwidth"), ncores = ncores, debug_output = TRUE, updateProgress = TRUE)
      
      
      rankMat.parsed=parse_edge_rank_matrix(edge_rank_matrix = generatematrices, edge_selection_strategy = "default",
                                            mat_weights = "rank", topN = 10, debug_output = TRUE, updateProgress = TRUE)
      
      #conGraph <- get_iGraph(rankMat.parsed$bin_mat)
      
      
      netlist[[exprmatrix]] <- rankMat.parsed$bin_mat
      
      
    }else{print(paste("Number of samples is too low for ", myfiles[exprmatrix], ": skipping!"))}
    
    
  }
  
  return(netlist)
}


#names(netlist) <- sapply(myfiles, function(x) strsplit(strsplit(x, "_")[[1]][2], ".", fixed=TRUE)[[1]][1])


#' Creates conversion tables between compounds names, CIDs, CASs and SMILEs
#' 
#' @param compounds_name A character vector of compound names.
#' @param XML_dest path where to store the XML files.
#' @param SDF_dest path where to store the SDF files.
#' @return A dataframe with colnames CMAP, CID, CAS, SMILE.
#' @examples
#' drugs_conv <- download_SDF_and_cid_from_compounds_names(compounds_name = names(L1000_networks), XML_dest = "/home/User/XML_folder", SDF_dest = "/home/User/SDF_folder")


download_SDF_and_cid_from_compounds_names = function(compounds_name, XML_dest = ".", SDF_dest = "."){
  
  if(is.null(compounds_name)) {
    stop("Error: please indicate a character vector of compound names!")
  }
  
  if(!class(compounds_name)=="character") {
    stop("Error: please provide an object of class character!")
  }
  
  
  info_mat = c()
  
  for(i in 1:length(compounds_name)){
    #for(i in 1:10){
    
    cid = NA
    CAS = NA
    result = tryCatch({
      
      # #This automatically download SDF files from compound names
      download.file(url = paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", compounds_name[i], "/SDF/?record_type=3d", sep=""),
                    destfile = paste(SDF_dest, "/", compounds_name[i],".sdf", sep=""))
      #
      #This automatically download XML files from compound names; I can extract from this the CID code
      
      download.file(url = paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", compounds_name[i], "/XML", sep=""),
                    destfile = paste(XML_dest,"/", compounds_name[i], ".xml", sep=""))
      data = xmlParse(paste(XML_dest, "/", compounds_name[i], ".xml", sep=""))
      xml_data <- xmlToList(data)
      
      cid = xml_data$`PC-Compound`$`PC-Compound_id`$`PC-CompoundType`$`PC-CompoundType_id`$`PC-CompoundType_id_cid`
      
      #This download a bigger XML files containing also the CAS number, I’m still not able to convert it but, I’ve found a website where you
      #can do a bulk conversion between CID and CAS: http://cts.fiehnlab.ucdavis.edu/conversion/batchConvert
      
      #download.file(url=paste(“https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/“,cid,“/XML/?response_type=display”,sep=“”)
      #              ,destfile = paste(“TGGATE_XML/“,compounds_name[i],“_2.xml”,sep=“”))
      
      library(RCurl)
      Myurl = getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/", cid, "/XML/?response_type=display", sep=""))
      html_data = htmlParse(Myurl)
      #xpat command that search for all the node namend “name” that contains the text “CAS”. For each of these nodes we take the following sibling named “stringvalue” that contains the CAS number
      nset = getNodeSet(html_data,"//name[text()='CAS']/following-sibling::*[1]")
      
      CAS = unique(unlist(lapply(nset,FUN = function(i){xmlValue(i)})))[1]
      
      nset = getNodeSet(html_data,"//tocheading[text()='Canonical SMILES']/following-sibling::*[3]")
      SMILE = unique(unlist(lapply(nset,FUN = function(i){xmlValue(i)})))[1]
      
      info_mat = rbind(info_mat, c(compounds_name[i],cid,CAS,SMILE))
      cat(i,"\n")
    }, warning = function(w) {
      #warning-handler-code
    }, error = function(e) {
      cat("Problems with ", compounds_name[i])
      
      #info_mat = rbind(info_mat, c(compounds_name[i],NA))
      
      #error-handler-code
    }, finally = {
      #cleanup-code
    })
    
  }
  
  info_mat = unique(info_mat)
  colnames(info_mat) = c("CMAP", "CID", "CAS", "SMILE")
  return(info_mat)
}


# drug_target_df = filter_by_L1000_drugs
# drugs_col = dat.drug.molecule_name
# targets_col = dat.target.gene_info.symbol


#' Creates a squared matrix of average shortes paths between couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param Graph An object of class "igraph".
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @return A squared matrix of average shortest paths between couples of drug targets.
#' @examples
#' 


get_avg_shortest_path <- function(Graph, drug_target_df, drugs_col, targets_col) {
  
  if(is.null(Graph) | !class(Graph)=="igraph"){
    stop("Error: please provide an igraph object!")
  }
  
  if(is.null(drug_target_df)) {
    stop("Error: please provide a dataframe with at least two columns: one for the drugs and another for the corresponding targets!")
  }
  
  if(is.null(drugs_col)) {
    stop("Error: please provide the column containing the drugs!")
  }
  
  if(is.null(targets_col)) {
    stop("Error: please provide the column containing the targets!")
  }
  
  if(!class(drugs_col)=="numeric") {
    stop("Error: please provide the index of the drug column!")
  }
  
  if(!class(targets_col)=="numeric") {
    stop("Error: please provide the index of the targets column!")
  }
  
  shortpath <- matrix(NA, nrow = length(unique(tolower(drug_target_df[, drugs_col]))), ncol = length(unique(tolower(drug_target_df[, drugs_col]))), dimnames = list(unique(tolower(drug_target_df[, drugs_col])), unique(tolower(drug_target_df[, drugs_col]))))
  
  
  for (drug1 in 1:length(unique(drug_target_df[, drugs_col]))) {
    print(paste("Computing avg shortest path for drug ", drug1 , ": ", unique(drug_target_df[, drugs_col])[drug1]))
    gene1 <- unique(as.vector(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% unique(tolower(drug_target_df[, drugs_col]))[drug1])]))
    
    for (drug2 in 1:length(unique(drug_target_df[, drugs_col]))){
      spathtotal <- c()
      gene2 <- unique(as.vector(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% unique(tolower(drug_target_df[, drugs_col]))[drug2])]))
      
      for(genes in 1:length(gene1)){
        spathgene <- igraph::get.shortest.paths(Graph, from = as.character(gene1)[genes], to = as.character(gene2))
        spathtotal <- c(spathtotal, mean(sapply(spathgene$vpath, function(x) length(x)-1)))
      }
      spathmean <- mean(spathtotal)
      shortpath[drug1,drug2] <- spathmean
    }
    
  }
  
  return(shortpath)
}


#' Creates a squared matrix of average jaccard indices between subsets of neighbors of a certain order of couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param Graph An object of class "igraph".
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param order Order of the neighborhood to be considered.
#' @return A squared matrix of average shortest paths between couples of drug targets.
#' @examples
#' 


get_avg_jacc_indices <- function(Graph, drug_target_df, drugs_col, targets_col, order = 1) {
  
  if(is.null(Graph) | !class(Graph)=="igraph") {
    stop("Error: please provide an igraph object!")
  }
  
  if(is.null(drug_target_df)) {
    stop("Error: please provide a dataframe with at least two columns: one for the drugs and another for the corresponding targets!")
  }
  
  if(is.null(drugs_col)) {
    stop("Error: please provide the column containing the drugs!")
  }
  
  if(is.null(targets_col)) {
    stop("Error: please provide the column containing the targets!")
  }
  
  if(!class(drugs_col)=="numeric") {
    stop("Error: please provide the index of the drug column!")
  }
  
  if(!class(targets_col)=="numeric") {
    stop("Error: please provide the index of the targets column!")
  }
  
  if(!class(order)=="numeric") {
    stop("Error: please provide the neighborhood order as a numeric value!")
  }
  
  jacc_drugs <- matrix(NA, nrow = length(unique(tolower(drug_target_df[, drugs_col]))), ncol = length(unique(tolower(drug_target_df[, drugs_col]))), dimnames = list(unique(tolower(drug_target_df[, drugs_col])), unique(tolower(drug_target_df[, drugs_col]))))
  
  for (drug1 in 1:length(unique(drug_target_df[, drugs_col]))) {
    
    #target drug1
    print(paste("Computing avg Jaccard index for drug ", drug1 , ": ", unique(drug_target_df[, drugs_col])[drug1]))
    gene1 <- unique(as.vector(drug_target_df[, targets_col][which(tolower(drug_target_df[, drugs_col]) %in% unique(tolower(drug_target_df[, drugs_col]))[drug1])]))
    
    for (drug2 in 1:length(unique(drug_target_df[, drugs_col]))){
      
      #targets drug2
      gene2 <- unique(as.vector(drug_target_df[, targets_col][which(tolower(drug_target_df[, drugs_col]) %in% unique(tolower(drug_target_df[, drugs_col]))[drug2])]))
      
      jacc_vec <- c()
      #for each target of drug1 
      for(genes1 in 1:length(gene1)){
        # find close genes
        gene1_first_degree <- igraph::ego(Graph, order = order, nodes = gene1[genes1])
        #for each target of drug2
        for (genes2 in 1:length(gene2)) {
          #find close genes
          gene2_first_degree <- igraph::ego(Graph, order = order, nodes = gene2[genes2])
          
          # jaccard index between the genes
          jacc_index <- length(intersect(as.character(names(gene1_first_degree[[1]])), as.character(names(gene2_first_degree[[1]]))))/length(union(as.character(names(gene1_first_degree[[1]])), as.character(names(gene2_first_degree[[1]]))))
          
          jacc_vec <- c(jacc_vec, jacc_index)
          
        }
      }
      jaccmean <- mean(jacc_vec)
      jacc_drugs[drug1, drug2] <- jaccmean
      
    }
    
  }
  return(jacc_drugs)
}


#' Creates a squared matrix of average jaccard indices between subsets of neighbors of a certain order of couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param graph_rank An ordered named vector of a graph's nodes. Typically obtained from the get_graph_named_ranklist function.
#' @return A rank of drugs based on the average ranks of their targets.
#' @examples
#' 

get_drug_avg_rank_from_graph <-  function(drug_target_df, drugs_col, targets_col, graph_rank) {
  
  if(is.null(drug_target_df)) {
    stop("Error: please provide a dataframe with at least two columns: one for the drugs and another for the corresponding targets!")
  }
  
  if(is.null(drugs_col)) {
    stop("Error: please provide the column containing the drugs!")
  }
  
  if(is.null(targets_col)) {
    stop("Error: please provide the column containing the targets!")
  }
  
  if(!class(drugs_col)=="numeric") {
    stop("Error: please provide the index of the drug column!")
  }
  
  if(!class(targets_col)=="numeric") {
    stop("Error: please provide the index of the targets column!")
  }
  
  if(is.null(graph_rank)) {
    stop("Error: please provide a named ranked vector of nodes!")
  }
    
    
   drug_inform_vec <- c()
  
  for(drugs in 1:length(unique(drug_target_df[, drugs_col]))) {
    targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col]==unique(drug_target_df[, drugs_col])[drugs])]))
    inform_rank <- mean(which(names(graph_rank) %in% targets))
    drug_inform_vec <- c(drug_inform_vec, inform_rank)
    
  }
  
  names(drug_inform_vec) <- unique(drug_target_df[, drugs_col])
  
  return(drug_inform_vec)
  
}


#' Creates a squared matrix of average jaccard indices between subsets of neighbors of a certain order of couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param Graph An object of class igraph.
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param graph_rank An ordered named vector of a graph's nodes. Typically obtained from the get_graph_named_ranklist function.
#' @param order Order of the neighborhood to be considered.
#' @return A rank of drugs based on the average ranks of their targets.
#' @examples
#' 

get_drug_coverage_score <- function(Graph, drug_target_df, drugs_col, targets_col, graph_rank, order = 1) {
  darios_drug_score <- c()
  target_jacc_list <- list()
  ji <- c()
  
  for (drugs in 1:length(unique(drug_target_df[, drugs_col]))) {
    target_jacc_list <- list()
    targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col]==unique(drug_target_df[, drugs_col])[drugs])]))
    print(paste("Computing graph coverage for", unique(drug_target_df[, drugs_col])[drugs], "based on", length(targets), "drug targets", sep = " "))
    inform_rank <- median(which(names(graph_rank) %in% targets))
    if(length(targets)==1) {
      ji <- 1
    }else{
      for (target in 1:length(targets)) {
        
        target_jacc_list[[target]] <- igraph::ego(Graph, order = order, nodes = targets[target])
        
      }
      
      for(jacc1 in 1:length(target_jacc_list)) {
        for(jacc2 in 1:length(target_jacc_list)) {
          ji <- c(ji, length(intersect(target_jacc_list[[jacc1]], target_jacc_list[[jacc2]]))/length(union(target_jacc_list[[jacc1]], target_jacc_list[[jacc2]])))
          
        }
        
      }
      
      
    }
    
    darios_drug_score <- c(darios_drug_score, length(targets)/inform_rank*median(ji))
    
  }
  
  names(darios_drug_score) <- unique(drug_target_df[, drugs_col])
  return(darios_drug_score)
  
}


#' Computes the tonino's coverage score based on the network coverage of input drugs and the network rank of target genes
#' 
#' @param Graph An object of class igraph.
#' @param drugs_input A vector of drugs' names.
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param graph_rank An ordered named vector of a graph's nodes. Typically obtained from the get_graph_named_ranklist function.
#' @param order Order of the neighborhood to be considered.
#' @return A rank of drugs based on the average ranks of their targets.
#' @examples
#' 

# get_drug_input_coverage <- function(Graph, drugs_input, drug_target_df, drugs_col, targets_col, graph_rank, order = 1) {
#   toninos_drug_score <- c()
#   target_neighbour_list <- list()
#   
#   targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% drugs_input)]))
#   for(target in 1:length(targets)){
#     target_neighbour_list[[target]] <- igraph::ego(Graph, order = order, nodes = targets[target])
#   }
#   
#   target_neighbour_list <- lapply(target_neighbour_list, function(x) igraph::as_ids(x[[1]]))
#   covered_nodes_overall <- unique(unlist(target_neighbour_list))
#   inform_rank <- median(which(names(graph_rank) %in% covered_nodes_overall))
#   coverage_overall <- length(covered_nodes_overall)/length(V(Graph))
#   toninos_drug_score <- coverage_overall*100/inform_rank #chiedi a Giovanni se va bene, altrimenti usiamo solo il coverage
#   return(toninos_drug_score)
# }



# get_drug_input_coverage <- function(Graph, drugs_input, drug_target_df, drugs_col, targets_col, graph_rank, order = 1, minCov=0, maxCov=1) {
#   toninos_drug_score <- c()
#   target_neighbour_list <- list()
#   targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% drugs_input)]))
#   for(target in 1:length(targets)){
#     target_neighbour_list[[target]] <- igraph::ego(Graph, order = order, nodes = targets[target])
#   }
#   target_neighbour_list <- lapply(target_neighbour_list, function(x) igraph::as_ids(x[[1]]))
#   covered_nodes_overall <- unique(unlist(target_neighbour_list))
#   # inform_rank <- median(which(names(graph_rank) %in% covered_nodes_overall))
#   coverage_overall <- length(covered_nodes_overall)/length(V(Graph))
#   toninos_drug_score <- coverage_overall #chiedi a Giovanni se va bene, altrimenti usiamo solo il coverage
#   return((toninos_drug_score - minCov)/ maxCov)
# }



random_targets <- function(targets, Graph){
  #1 prendo il vettore VT di N degrees in degree_bin per ogni target presente nel grafo
  #2 Divido i nodi di Graph bin in bins da almeno 200 nodi l'uno base al degree
  #3 Calcolo per ogni nodo in VT il bin dove cade
  #4 estraggio N nodi a cazzo dal bin dove il nodo cadeva
  deg_dist <- tibble(vertex=as_ids(V(Graph)), dist = igraph::degree(Graph)) %>% arrange(dist)
  # almeno 200 elementi a bin
  binned <-  binr::bins(deg_dist$dist,target.bins = 23,minpts = 200)
  deg_dist <- deg_dist %>% dplyr::mutate(BIN = findInterval(dist,binr::bins.getvals(binned))) %>%  arrange(BIN)
  deg_dist %>% group_by(BIN) %>% dplyr::summarise(n=n()) %>% as.data.frame()
  drugNodes_bins <- deg_dist %>% filter(vertex %in% targets) %>% group_by(BIN) %>% dplyr::summarise(n=n())
  rand_drug <- drugNodes_bins %>% pmap_df(function(BIN, n, deg_dist){deg_dist[deg_dist$BIN == BIN,] %>% sample_n(n)}, deg_dist)
  rand_drug$vertex
}


compute_cov_score <- function(targets, Graph, order, minCov, maxCov) {
  target_neighbour_list <- list()
  for(target in 1:length(targets)){
    target_neighbour_list[[target]] <- igraph::ego(Graph, order = order, nodes = targets[target])
  }
  target_neighbour_list <- lapply(target_neighbour_list, function(x) igraph::as_ids(x[[1]]))
  covered_nodes_overall <- unique(unlist(target_neighbour_list))
  # inform_rank <- median(which(names(graph_rank) %in% covered_nodes_overall))
  coverage_overall <- length(covered_nodes_overall)/length(V(Graph))
  toninos_drug_score <- coverage_overall #chiedi a Giovanni se va bene, altrimenti usiamo solo il coverage
  return((toninos_drug_score - minCov)/ maxCov)
}




get_drug_input_coverage <- function(Graph, drugs_input, drug_target_df, drugs_col, targets_col, graph_rank, order = 1, minCov=0, maxCov=1) {
  toninos_drug_score <- c()
  targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% drugs_input)]))
  observed <- compute_cov_score(targets, Graph, order, minCov, maxCov)
  simulated <- rep(0,100)
  for(i in 1:100){
    random_targ <- random_targets(targets, Graph)
    simulated[i] <- compute_cov_score(random_targ, Graph, order, minCov, maxCov)
  }
  z_score = (observed - mean(simulated))/sd(simulated)
  pvalue = pnorm(-abs(z_score))
  final_list <- (list(z_score, pvalue))
  return(final_list)
}




toninos_waterfall <- function (.data = NULL, values, labels, rect_text_labels = values, 
                               rect_text_size = 1, rect_text_color="white", rect_text_labels_anchor = "centre", put_rect_text_outside_when_value_below = 0.05 * 
                                 (max(cumsum(values)) - min(cumsum(values))), calc_total = FALSE, 
                               total_axis_text = "Total", total_rect_text = sum(values), 
                               total_rect_color = "black", total_rect_text_color = "white", 
                               fill_colours = NULL, fill_by_sign = TRUE, rect_width = 0.7, 
                               rect_border = "black", draw_lines = TRUE, lines_anchors = c("right", 
                                                                                           "left"), linetype = "dashed", draw_axis.x = "behind", 
                               theme_text_family = "", scale_y_to_waterfall = TRUE, print_plot = FALSE, 
                               ggplot_object_name = "mywaterfall") 
{
  if (!is.null(.data)) {
    if (!is.data.frame(.data)) {
      stop("`.data` was a ", class(.data)[1], ", but must be a data.frame.")
    }
    if (ncol(.data) < 2L) {
      stop("`.data` had fewer than two columns, yet two are required: labels and values.")
    }
    dat <- as.data.frame(.data)
    char_cols <- vapply(dat, is.character, FALSE)
    factor_cols <- vapply(dat, is.factor, FALSE)
    num_cols <- vapply(dat, is.numeric, FALSE)
    if (!xor(num_cols[1], num_cols[2]) || sum(char_cols[1:2], 
                                              factor_cols[1:2], num_cols[1:2]) != 2L) {
      const_width_name <- function(noms) {
        if (is.data.frame(noms)) {
          noms <- names(noms)
        }
        max_width <- max(nchar(noms))
        formatC(noms, width = max_width)
      }
      stop("`.data` did not contain exactly one numeric column and exactly one character or factor ", 
           "column in its first two columns.\n\t", "1st column: '", 
           const_width_name(dat)[1], "'\t", sapply(dat, 
                                                   class)[1], "\n\t", "2nd column: '", const_width_name(dat)[2], 
           "'\t", sapply(dat, class)[2])
    }
    if (num_cols[1L]) {
      .data_values <- .subset2(dat, 1L)
      .data_labels <- .subset2(dat, 2L)
    }
    else {
      .data_values <- .subset2(dat, 2L)
      .data_labels <- .subset2(dat, 1L)
    }
    if (!missing(values) && !missing(labels)) {
      warning(".data and values and labels supplied, .data ignored")
    }
    else {
      values <- .data_values
      labels <- as.character(.data_labels)
    }
  }
  if (!(length(values) == length(labels) && length(values) == 
        length(rect_text_labels))) {
    stop("values, labels, fill_colours, and rect_text_labels must all have same length")
  }
  if (rect_width > 1) 
    warning("rect_Width > 1, your chart may look terrible")
  number_of_rectangles <- length(values)
  north_edge <- cumsum(values)
  south_edge <- c(0, cumsum(values)[-length(values)])
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
  }
  if (fill_by_sign) {
    if (!is.null(fill_colours)) {
      warning("fill_colours is given but fill_by_sign is TRUE so fill_colours will be ignored.")
    }
    fill_colours <- ifelse(values >= 0, gg_color_hue(2)[2], 
                           gg_color_hue(2)[1])
  }
  else {
    if (is.null(fill_colours)) {
      fill_colours <- gg_color_hue(number_of_rectangles)
    }
  }
  if (!(grepl("^[lrc]", lines_anchors[1]) && grepl("^[lrc]", 
                                                   lines_anchors[2]))) 
    stop("lines_anchors must be a pair of any of the following: left, right, centre")
  if (grepl("^l", lines_anchors[1])) 
    anchor_left <- rect_width/2
  if (grepl("^c", lines_anchors[1])) 
    anchor_left <- 0
  if (grepl("^r", lines_anchors[1])) 
    anchor_left <- -1 * rect_width/2
  if (grepl("^l", lines_anchors[2])) 
    anchor_right <- -1 * rect_width/2
  if (grepl("^c", lines_anchors[2])) 
    anchor_right <- 0
  if (grepl("^r", lines_anchors[2])) 
    anchor_right <- rect_width/2
  if (!calc_total) {
    p <- if (scale_y_to_waterfall) {
      ggplot2::ggplot(data.frame(x = c(labels, labels), 
                                 y = c(south_edge, north_edge)), ggplot2::aes_string(x = "x", 
                                                                                     y = "y"))
    }
    else {
      ggplot2::ggplot(data.frame(x = labels, y = values), 
                      ggplot2::aes_string(x = "x", y = "y"))
    }
    p <- p + ggplot2::geom_blank() + ggplot2::theme(axis.title = ggplot2::element_blank())
  }
  else {
    p <- if (scale_y_to_waterfall) {
      ggplot2::ggplot(data.frame(x = c(labels, total_axis_text, 
                                       labels, total_axis_text), y = c(south_edge, north_edge, 
                                                                       south_edge[number_of_rectangles], north_edge[number_of_rectangles])), 
                      ggplot2::aes_string(x = "x", y = "y"))
    }
    else {
      ggplot2::ggplot(data.frame(x = c(labels, total_axis_text), 
                                 y = c(values, north_edge[number_of_rectangles])), 
                      ggplot2::aes_string(x = "x", y = "y"))
    }
    p <- p + ggplot2::geom_blank() + ggplot2::theme(axis.title = ggplot2::element_blank())
  }
  if (grepl("behind", draw_axis.x)) {
    p <- p + ggplot2::geom_hline(yintercept = 0)
  }
  for (i in seq_along(values)) {
    p <- p + ggplot2::annotate("rect", xmin = i - rect_width/2, 
                               xmax = i + rect_width/2, ymin = south_edge[i], ymax = north_edge[i], 
                               colour = rect_border, fill = fill_colours[i])
    if (i > 1 && draw_lines) {
      p <- p + ggplot2::annotate("segment", x = i - 1 - 
                                   anchor_left, xend = i + anchor_right, linetype = linetype, 
                                 y = south_edge[i], yend = south_edge[i])
    }
  }
  for (i in seq_along(values)) {
    if (abs(values[i]) > put_rect_text_outside_when_value_below) {
      p <- p + ggplot2::annotate("text", x = i, y = 0.5 * 
                                   (north_edge[i] + south_edge[i]), family = theme_text_family, 
                                 label = ifelse(rect_text_labels[i] == values[i], 
                                                ifelse(values[i] < 0, paste0("−", -1 * values[i]), 
                                                       values[i]), rect_text_labels[i]), size = rect_text_size/(5/14), colour=rect_text_color)
    }
    else {
      p <- p + ggplot2::annotate("text", x = i, y = north_edge[i], 
                                 family = theme_text_family, label = ifelse(rect_text_labels[i] == 
                                                                              values[i], ifelse(values[i] < 0, paste0("−", 
                                                                                                                      -1 * values[i]), values[i]), rect_text_labels[i]), 
                                 vjust = ifelse(values[i] >= 0, -0.2, 1.2), size = rect_text_size/(5/14), colour="black")
    }
  }
  if (calc_total) {
    p <- p + ggplot2::annotate("rect", xmin = number_of_rectangles + 
                                 1 - rect_width/2, xmax = number_of_rectangles + 1 + 
                                 rect_width/2, ymin = 0, ymax = north_edge[number_of_rectangles], 
                               colour = rect_border, fill = total_rect_color) + 
      ggplot2::annotate("text", x = number_of_rectangles + 
                          1, y = 0.5 * north_edge[number_of_rectangles], 
                        family = theme_text_family, label = ifelse(total_rect_text == 
                                                                     sum(values), ifelse(north_edge[number_of_rectangles] < 
                                                                                           0, paste0("−", -1 * north_edge[number_of_rectangles]), 
                                                                                         north_edge[number_of_rectangles]), total_rect_text), 
                        color = total_rect_text_color, size = rect_text_size/(5/14)) + 
      ggplot2::scale_x_discrete(labels = c(labels, total_axis_text))
    if (draw_lines) {
      p <- p + ggplot2::annotate("segment", x = number_of_rectangles - 
                                   anchor_left, xend = number_of_rectangles + 1 + 
                                   anchor_right, y = north_edge[number_of_rectangles], 
                                 yend = north_edge[number_of_rectangles], linetype = linetype)
    }
  }
  else {
    p <- p + ggplot2::scale_x_discrete(labels = labels)
  }
  if (grepl("front", draw_axis.x)) {
    p <- p + ggplot2::geom_hline(yintercept = 0)
  }
  if (print_plot) {
    if (ggplot_object_name %in% ls(.GlobalEnv)) 
      warning("Overwriting ", ggplot_object_name, " in global environment.")
    assign(ggplot_object_name, p, inherits = TRUE)
    print(p)
  }
  else {
    return(p)
  }
}



#' Creates a squared matrix of average jaccard indices between subsets of neighbors of a certain order of couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param Graph An object of class igraph.
#' @param drugs_rank Dataframe of at least two columns containing drug names and corresponding targets.
#' @param nDrugs Column index of the dataframe drug_target_df where the drugs are stored.
#' @param order Order of the neighborhood to be considered.
#' @return A rank of drugs based on the additive network coverage. This function also returns a plot showing the additive coverage of the network.
#' @examples additive_coverage <- get_additive_graph_coverage(Graph = conGraph, drugs_rank = finaldistrank, nDrugs = 5, order = 1)
#' 

get_additive_graph_coverage <- function(Graph, drugs_rank, nDrugs, order=1) {
  
  target_areas_list <- list()
  drugcoverage <- c()
  cumcov <- c()
  
  for(i in 1:nDrugs) {
    
    target_areas <- c()
    targets <- as.vector(unique(filter_by_L1000_drugs$dat.target.gene_info.symbol[which(filter_by_L1000_drugs$dat.drug.molecule_name==drugs_rank[i])]))
    for (drugs in 1:length(targets)) {
      firstdegree_cov <- igraph::ego(graph = Graph, order = order, nodes = targets[drugs])
      target_areas <- c(target_areas, firstdegree_cov[[1]])
      
    }
    
    target_areas_list[[i]] <- unique(target_areas)
    drugcoverage <- c(drugcoverage, length(unique(unlist(target_areas_list[1:i])))/length(igraph::V(Graph)))
    print(i)
    
  }
  
  names(drugcoverage) <- drugs_rank[1:i]
  add_cov <- data.frame(drugcoverage=drugcoverage, drug=names(drugcoverage))
  # plot(drugcoverage, type = "l")
  # p <- ggplot2::ggplot(data = add_cov, aes(drug, drugcoverage, colour="red", group=1)) +
  #   geom_line(stat = "identity") +
  #   scale_x_discrete(name ="Drug",
  #                      limits=rownames(add_cov))
  # print(p)
  
  add_cov_value <- c(drugcoverage[1])

  for (cov in 2:length(drugcoverage)){add_cov_value <- c(add_cov_value, drugcoverage[cov]-drugcoverage[cov-1])}

  add_cov_df <- data.frame(labels=names(drugcoverage), values=add_cov_value)

  rect_col <- rep("grey38", length(add_cov_df$values))

  toninos_waterfall(.data = add_cov_df, labels = rownames(add_cov_df), calc_total = FALSE, rect_text_labels = round(add_cov_df$values, digits = 2),
                               print_plot = TRUE, fill_by_sign = FALSE, fill_colours = rect_col, rect_text_size = 1, ggplot_object_name = "ggplot_obj", rect_text_color = "white")
  #waterfalls::waterfall(.data = add_cov_df, labels = rownames(add_cov_df), calc_total = FALSE, rect_text_labels = round(add_cov_df$values, digits = 2),
   #                     print_plot = TRUE, fill_by_sign = FALSE, fill_colours = rect_col, rect_text_size = 1, ggplot_object_name = "ggplot_obj")

  return(drugcoverage)
}






