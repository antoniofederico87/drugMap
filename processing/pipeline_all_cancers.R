rm(list=ls())

project_flag <- "TCGA-BRCA"
cell_type <- "MCF7"
mainDir <- "."
subDir <- paste0("./secondtest/", project_flag, "/")

if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))

}


options(stringsAsFactors = FALSE)
source("functions_pipeline.R")

devtools::source_url("https://raw.githubusercontent.com/Greco-Lab/INfORM/master/INfORM_functions.R")

query_NT_paired <- GDCquery(project = project_flag,
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM-UQ", sample.type = "Solid Tissue Normal")


samplesDown_NT <- getResults(query_NT_paired,cols=c("cases"))

newvec <- sapply(samplesDown_NT, function(x) paste0(substr(x, start = 1, stop = 13), "*"))


query_TP_paired <- GDCquery(project = project_flag,
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM-UQ", barcode = newvec, sample.type = "Primary solid Tumor")


samplesDown_TP <- getResults(query_TP_paired,cols=c("cases"))


### First check: if all the control samples have the tumoral counterpart, and viceversa ###

if(length(samplesDown_NT)!=length(samplesDown_TP)){
  short_barcodes_NT <- sapply(samplesDown_NT, function(x) substr(x, start = 1, stop = 13))
  short_barcodes_TP <- sapply(samplesDown_TP, function(x) substr(x, start = 1, stop = 13))
  if(length(setdiff(short_barcodes_NT, short_barcodes_TP))>0){
    samplesDown_NT <- samplesDown_NT[-which(short_barcodes_NT %in% setdiff(short_barcodes_NT, short_barcodes_TP))]
  }else if(length(setdiff(short_barcodes_TP, short_barcodes_NT))>0) {
    samplesDown_TP <- samplesDown_TP[-which(short_barcodes_TP %in% setdiff(short_barcodes_TP, short_barcodes_NT))]
  }
}


### Second check: remove all the duplicated samples in the tumor or normal counterparts ###

if(any(table(short_barcodes_TP)>1)){
  for(i in 1:length(names(which(table(short_barcodes_TP)>1)))){
    duplicated_samples <- samplesDown_TP[grep(samplesDown_TP, pattern = names(which(table(short_barcodes_TP)!=1))[i])]
    samplesDown_TP <- samplesDown_TP[-which(samplesDown_TP %in% duplicated_samples)]
    samplesDown_TP <- append(x = samplesDown_TP, values = duplicated_samples[1])
    }
     
  }else if(any(table(short_barcodes_NT)>1)){
    for(k in 1:length(names(which(table(short_barcodes_NT)>1)))){
      duplicated_samples <- samplesDown_NT[grep(samplesDown_NT, pattern = names(which(table(short_barcodes_NT)!=1))[i])]
      samplesDown_NT <- samplesDown_NT[-which(samplesDown_NT %in% duplicated_samples)]
      samplesDown_NT <- append(x = samplesDown_NT, values = duplicated_samples[1])
    }
  }

print(paste("The number of normal samples for", project_flag, "before checking for duplicates is:", length(samplesDown_NT)))
print(paste("The number of tumoral samples for", project_flag, "before checking for duplicates is:", length(samplesDown_TP)))


query_NT_paired <- GDCquery(project = project_flag,
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM-UQ", sample.type = "Solid Tissue Normal", barcode = samplesDown_NT)


samplesDown_NT <- getResults(query_NT_paired,cols=c("cases"))

newvec <- sapply(samplesDown_NT, function(x) paste0(substr(x, start = 1, stop = 13), "*"))


query_TP_paired <- GDCquery(project = project_flag,
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM-UQ", sample.type = "Primary solid Tumor", barcode = samplesDown_TP)


samplesDown_TP <- getResults(query_TP_paired,cols=c("cases"))

print(paste("The number of normal samples for", project_flag, "after checking for duplicates is:", length(samplesDown_NT)))
print(paste("The number of tumoral samples for", project_flag, "after checking for duplicates is:", length(samplesDown_TP)))


GDCdownload(query_TP_paired, directory = "Tumor")
GDCdownload(query_NT_paired, directory = "Normal")


pathnormal <- paste0("./secondtest/", project_flag, "/Normal/", project_flag, "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification")
pathtumor <- paste0("./secondtest/", project_flag, "/Tumor/", project_flag, "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification")



### File di conversione filenames/barcodes###

##TUMOR unzip##
barcodes_tumor <- data.frame(query_TP_paired$results[[1]]$file_name, query_TP_paired$results[[1]]$cases)
colnames(barcodes_tumor) <- c("filename", "barcode")


list1 <- list.files(path = pathtumor, recursive = TRUE)
for (i in 1:length(list1))
  R.utils::gunzip(paste0(pathtumor,"/", list1[i]), destname=paste0(pathtumor,"/",strsplit(gsub(".gz","",list1[i]),"/")[[1]][2]))

#setwd("./BRCA_data/") #uncomment in taito
list.files(".")


#unlink(list.dirs(pathtumor), recursive = TRUE) #elimina cartelle (equivalente al comando rmdir)


##NORMAL unzip##

barcodes_normal <- data.frame(query_NT_paired$results[[1]]$file_name, query_NT_paired$results[[1]]$cases)
colnames(barcodes_normal) <- c("filename", "barcode")


list1 <- list.files(path = pathnormal, recursive = TRUE)
for (i in 1:length(list1))
  R.utils::gunzip(paste0(pathnormal,"/", list1[i]), destname=paste0(pathnormal,"/",strsplit(gsub(".gz","",list1[i]),"/")[[1]][2]))

#setwd("./BRCA_data/") #uncomment in taito
list.files(pathnormal)


#unlink(list.dirs(pathnormal), recursive = TRUE) #elimina cartelle (equivalente al comando rmdir)


###TUMOR matrix###
setwd(pathtumor)
myfiles <- sort(list.files(pathtumor, pattern = "*.FPKM-UQ.txt"))
sample1 <- read.table(myfiles[1], header = F, sep = "\t", quote = "")
rownames(sample1) <- sample1$V1
tumor <- data.frame(rownames(sample1), sample1[,2])


for (file in 2:length(myfiles)) {
  sample <- read.table(myfiles[file], header = FALSE, sep = "\t", quote = "")
  tumor <- cbind(tumor, sample[,2])
}

tumor[,1] <- NULL
row.names(tumor) <- sample$V1
colnames(tumor) <- myfiles
head(tumor)


### Format filenames ###

dummy <- c()
newvec <- c()

for (i in 1:length(barcodes_tumor$filename)) {
  dummy <- strsplit(as.character(barcodes_tumor$filename[i]), ".", fixed = TRUE)[[1]][1]
  newvec <- c(newvec, dummy)
}

barcodes_tumor$filename <- newvec
head(barcodes_tumor$filename)



dummy <- c()
newvec2 <- c()


for (k in 1:length(colnames(tumor))) {
  dummy <- strsplit(colnames(tumor)[k], ".", fixed = TRUE)[[1]][1]
  newvec2 <- c(newvec2, dummy)
}


barcodes_tumor[match(newvec2, barcodes_tumor$filename),]$barcode

colnames(tumor) <- barcodes_tumor[match(newvec2, barcodes_tumor$filename),]$barcode


###NORMAL matrix###
setwd(pathnormal)
myfiles <- sort(list.files(pathnormal, pattern = "*.FPKM-UQ.txt"))
sample1 <- read.table(myfiles[1], header = F, sep = "\t", quote = "")
rownames(sample1) <- sample1$V1
normal <- data.frame(rownames(sample1), sample1[,2])


for (file in 2:length(myfiles)) {
  sample <- read.table(myfiles[file], header = FALSE, sep = "\t", quote = "")
  normal <- cbind(normal, sample[,2])
}

normal[,1] <- NULL
row.names(normal) <- sample$V1
colnames(normal) <- myfiles
head(normal)


### Format filenames ###

dummy <- c()
newvec <- c()

for (i in 1:length(barcodes_normal$filename)) {
  dummy <- strsplit(as.character(barcodes_normal$filename[i]), ".", fixed = TRUE)[[1]][1]
  newvec <- c(newvec, dummy)
}

barcodes_normal$filename <- newvec
head(barcodes_normal$filename)



dummy <- c()
newvec2 <- c()


for (k in 1:length(colnames(normal))) {
  dummy <- strsplit(colnames(normal)[k], ".", fixed = TRUE)[[1]][1]
  newvec2 <- c(newvec2, dummy)
}


barcodes_normal[match(newvec2, barcodes_normal$filename),]$barcode

colnames(normal) <- barcodes_normal[match(newvec2, barcodes_normal$filename),]$barcode




# setdiff(rownames(normal), rownames(tumor))
# ix <- which(c("ENSG00000270112.3", "ENSG00000167578.15") %in% rownames(normal))
# normal <- normal[-ix,]

### Wilcoxon and var.test between Tumor and Control samples ###

all_samples <- cbind(tumor, normal)

keep <- rowSums(all_samples)>0
all_samples_filtered <- all_samples[keep,]
all_samples_filtered <- all_samples_filtered+1

### Collapse at gene level ###

ensemblid <- rownames(all_samples_filtered)

ensemblid <- sapply(strsplit(ensemblid, ".", fixed = TRUE), "[", 1)

rownames(all_samples_filtered) <- ensemblid

genesym <- mapIds(org.Hs.eg.db, keys = ensemblid, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

all_samples_filtered$genesym <- genesym
all_samples_filtered <- aggregate.data.frame(all_samples_filtered, by = list(all_samples_filtered$genesym), FUN = "mean")
rownames(all_samples_filtered) <- all_samples_filtered$Group.1
all_samples_filtered[,c("Group.1", "genesym")] <- NULL


inputinform <- get_informative_genes(expr_mat = all_samples_filtered, group = factor(as.factor(c(rep("T", length(colnames(tumor))), rep("N", length(colnames(normal)))))), test = "wilcoxon + var.test", percentile = 10)

saveRDS(inputinform, file=paste0(mainDir, subDir, "inputinform_", project_flag, ".rds"))


generatematrices=get_ranked_consensus_matrix(gx_table = inputinform[[1]], iMethods = c("clr","aracne","mrnet"),
                                             iEst = c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
                                             iDisc=c("none","equalfreq","equalwidth","globalequalwidth"), ncores = 35, debug_output = TRUE, updateProgress = TRUE)


#Parse ranked matrix and get bin_mat and edge_rank
# Get edge rank list and binary inference matrix from edge rank matrix computed by get_ranked_consensus_matrix().
# parse_edge_rank_matrix parses the edge rank matrix created by using the internal function get_ranked_consensus_matrix_matrix() to get a ranked edge list and a binary matrix.

rankMat.parsed=parse_edge_rank_matrix(edge_rank_matrix = generatematrices, edge_selection_strategy = "default",
                                      mat_weights = "rank", topN = 10, debug_output = TRUE, updateProgress = TRUE)


conGraph <- get_iGraph(rankMat.parsed$bin_mat)
saveRDS(conGraph, file=paste0(mainDir, subDir, "conGraph_", project_flag, ".rds"))

#conGraph.modules <- get_modules(conGraph, method="walktrap")

graph_rank <- get_graph_named_ranklist(Graph = conGraph, ranked_vectors = list(inputinform[["var_score"]], inputinform[["wilcox_score"]]))

#json_opentargets_parsed <- parse_opentargets_json(json_file = "/home/antonio/19.02_evidence_data.json")
json_opentargets_parsed <- readRDS("/home/antonio/final_opentargets_parsed.rds")

map_opentargets <- map_opentargets_on_graph(Graph = conGraph, parsed_opentargets = json_opentargets_parsed)


L1000_GSE92742_pert_info <- read.delim("./BRCA_data_paired_wm/GSE92742_Broad_LINCS_pert_info.txt", header=TRUE, quote = "")
toRemove <- c(1,3,4,5,6,8)
L1000_GSE92742_pert_info <- L1000_GSE92742_pert_info[, -toRemove]
L1000_GSE70138_pert_info <- read.delim("./BRCA_data_paired_wm/GSE70138_Broad_LINCS_pert_info.txt", header = TRUE, quote = "")
toRemove <- c(1,3,5)
L1000_GSE70138_pert_info <- L1000_GSE70138_pert_info[, -toRemove]
L1000_GSE70138_pert_info <- L1000_GSE70138_pert_info[, c(2,1)]

L1000_all_drugs <- rbind(L1000_GSE92742_pert_info, L1000_GSE70138_pert_info)
L1000_all_drugs$pert_iname <- sapply(L1000_all_drugs$pert_iname, tolower)

map_opentargets$dat.drug.molecule_name <- gsub(map_opentargets$dat.drug.molecule_name, pattern = ".", replacement = "-", fixed = TRUE)
filter_by_L1000_drugs <- map_opentargets[which(map_opentargets$dat.drug.molecule_name %in% L1000_all_drugs$pert_iname),]
filter_by_L1000_drugs$canonical_smiles <- NA
filter_by_L1000_drugs$canonical_smiles <- L1000_all_drugs$canonical_smiles[match(filter_by_L1000_drugs$dat.drug.molecule_name, L1000_all_drugs$pert_iname)]
filter_by_L1000_drugs <- unique(filter_by_L1000_drugs)

col_meta_path <- "./LINCS_project/data/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt"
col_meta <- read.delim(col_meta_path, sep="\t", stringsAsFactors=F)
col_path <- "./LINCS_project/data/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx"
GSE70138_metadata <- read.gctx.meta(col_path, dimension = "column")

col_meta_filtered <- col_meta %>%
  filter(cell_id==cell_type, pert_idose=="10.0 um", pert_itime=="24 h") %>%
  filter(pert_iname %in% filter_by_L1000_drugs$dat.drug.molecule_name)

col_meta_filtered_CTRL <- col_meta %>%
  filter(cell_id==cell_type, pert_itime=="24 h", pert_iname=="DMSO")

col_meta_path_2 <- "./LINCS_project/data/GSE92742_Broad_LINCS_sig_info.txt"
col_meta_2 <- read.delim(col_meta_path_2, sep = "\t", stringsAsFactors = F)
col_path_2 <- "./LINCS_project/data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"

col_meta_filtered_2 <- col_meta_2 %>%
  filter(cell_id==cell_type, pert_idose=="10 µM", pert_itime=="24 h") %>%
  filter(pert_iname %in% filter_by_L1000_drugs$dat.drug.molecule_name) 

col_meta_filtered_2_CTRL <- col_meta_2 %>%
  filter(cell_id==cell_type, pert_itime=="24 h", pert_iname=="DMSO")

additional_drugs <- setdiff(filter_by_L1000_drugs$dat.drug.molecule_name, union(col_meta_filtered$pert_iname, col_meta_filtered_2$pert_iname))


col_meta_filtered_all <- rbind(col_meta_filtered, col_meta[which(col_meta$pert_iname %in% names(additional_drugs)),])
col_meta_filtered_all <- rbind(col_meta_filtered_all, col_meta[which(col_meta$pert_iname %in% additional_drugs),])

col_meta_filtered_all <- col_meta_filtered_all %>%
  filter(cell_id==cell_type, pert_idose=="10.0 um", pert_itime=="24 h")

col_meta_filtered_2_all <- rbind(col_meta_filtered_2, col_meta_2[which(col_meta_2$pert_iname %in% names(additional_drugs)),])
col_meta_filtered_2_all <- rbind(col_meta_filtered_2_all, col_meta_2[which(col_meta_2$pert_iname %in% additional_drugs),])

col_meta_filtered_2_all <- col_meta_filtered_2_all %>%
  filter(cell_id==cell_type, pert_idose=="10 µM", pert_itime=="24 h")

col_meta_filtered_2_all[, c(6,7,9,10)] <- NULL
col_meta_filtered_2_CTRL[, c(6,7,9,10)] <- NULL

col_meta_filtered_all <- tidyr::separate_rows(col_meta_filtered_all, distil_id, sep = "\\|")
col_meta_filtered_all <- unique(col_meta_filtered_all)
col_meta_filtered_2_all <- tidyr::separate_rows(col_meta_filtered_2_all, distil_id, sep = "\\|")
col_meta_filtered_2_all <- unique(col_meta_filtered_2_all)

col_meta_filtered_CTRL <- tidyr::separate_rows(col_meta_filtered_CTRL, distil_id, sep = "\\|")
col_meta_filtered_2_CTRL <- tidyr::separate_rows(col_meta_filtered_2_CTRL, distil_id, sep = "\\|")


GSE70138_ds <- parse.gctx(col_path, cid=col_meta_filtered_all$distil_id)
GSE70138_CTRL <- parse.gctx(col_path, cid=col_meta_filtered_CTRL$distil_id)
GSE92742_ds <- parse.gctx(col_path_2, cid = col_meta_filtered_2_all$distil_id)
GSE92742_CTRL <- parse.gctx(col_path_2, cid=unique(col_meta_filtered_2_CTRL$distil_id))

entire_gctx <- merge.gct(GSE70138_ds, GSE92742_ds, dimension = "column")
entire_gctx_CTRL <- merge.gct(GSE70138_CTRL, GSE92742_CTRL, dimension = "column")

col_meta_filtered_total <- rbind(col_meta_filtered_all, col_meta_filtered_2_all)
col_meta_filtered_total_CTRL <- rbind(col_meta_filtered_CTRL, col_meta_filtered_2_CTRL)

landmark_genes <- read.delim("./LINCS_project/data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", header=TRUE, quote="")
landmark_genes <- as.character(landmark_genes$pr_gene_id)

create_L1000_exprmat <- build_drug_expression_matrices(L1000_data = entire_gctx, L1000_metadata = col_meta_filtered_total, include_controls = FALSE, dest_path = paste0(mainDir, subDir,"Expression_matrices"), inputgenes = landmark_genes)

# infer_L1000_coexpr_networks <- get_multiple_coexpression_networks(path ="./Expression_matrices/", minSamples = 5, ncores = 40)
# 
# names(infer_L1000_coexpr_networks) <- sapply(list.files("./Expression_matrices"), function(x) strsplit(strsplit(x, "_")[[1]][3], ".", fixed=TRUE)[[1]][1])[1:414]
# 
# infer_L1000_coexpr_networks <- infer_L1000_coexpr_networks[-which(lapply(infer_L1000_coexpr_networks, class)=="NULL")]

expr_mat <- list.files(paste0(mainDir, subDir,"Expression_matrices"), pattern = "Expression_")

for (file in 1:length(expr_mat)){
  minSamples <- 5
  drugname <- strsplit(strsplit(expr_mat[file], "_")[[1]][3], ".", fixed=TRUE)[[1]][1]
  if(file.exists(paste0(mainDir, subDir, "Expression_matrices/", drugname, ".rds"))){
    print(paste("Adjacency matrix for", drugname, "already present!"))
    next
  }else{
    
    input <- data.frame()
    
      #for (exprmatrix in 1:5) {  
      print(paste("Processing: ", expr_mat[file]))
      print(paste("Drug", file, "of", length(expr_mat)))
      
      input <- read.delim(paste0(mainDir, subDir, "Expression_matrices/", expr_mat[file]), header = TRUE, quote = "", row.names = 1)
      if(length(which(apply(input, 1, var)==0)!=0)){input <- input[-which(apply(input, 1, var)==0),]}
      
      if (length(colnames(input))>=minSamples) {    # Filter out the matrices having less than minSamples samples
        
        #controls <- inputinform[, grep(pattern = "CTRL", x = colnames(inputinform))]
        
        generatematrices=get_ranked_consensus_matrix(gx_table = input, iMethods = c("clr","aracne","mrnet"),
                                                     iEst = c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
                                                     iDisc=c("none","equalfreq","equalwidth","globalequalwidth"), ncores = 40, debug_output = TRUE, updateProgress = TRUE)
        
        
        rankMat.parsed=parse_edge_rank_matrix(edge_rank_matrix = generatematrices, edge_selection_strategy = "default",
                                              mat_weights = "rank", topN = 10, debug_output = TRUE, updateProgress = TRUE)
        
        conGraph_drug <- get_iGraph(rankMat.parsed$bin_mat)
  
        saveRDS(conGraph_drug, file = paste0(mainDir, subDir, "Expression_matrices/", drugname, ".rds"))

      }else{print(paste("Number of samples is too low for ", expr_mat[file], ": skipping!"))}
      
    }
    
}


all_drugs_vec <- sapply(expr_mat, function(x) strsplit(strsplit(x, "_")[[1]][3], ".", fixed=TRUE)[[1]][1])
infer_L1000_coexpr_networks <- rep(list(NULL), length(all_drugs_vec))

for (net in 1:length(all_drugs_vec)){
  if (file.exists(paste0(mainDir, subDir,"Expression_matrices/", all_drugs_vec[net], ".rds"))){
    infer_L1000_coexpr_networks[[net]] <- readRDS(paste0(mainDir, subDir, "Expression_matrices/", all_drugs_vec[net], ".rds"))  
  }
  
}

names(infer_L1000_coexpr_networks) <- all_drugs_vec

infer_L1000_coexpr_networks <- infer_L1000_coexpr_networks[-which(sapply(infer_L1000_coexpr_networks, function(x) class(x))=="NULL")] ### This line gives error if no NULL is found
saveRDS(infer_L1000_coexpr_networks, file=paste0(mainDir, subDir, "infer_L1000_coexpr_networks_", project_flag, ".rds"))

infer_L1000_coexpr_networks <- readRDS(paste0(mainDir, subDir, "infer_L1000_coexpr_networks_", project_flag, ".rds"))
### Function from the nettools package to compute a distance measure among adjacency matrices

ddist <- nettools::netdist(infer_L1000_coexpr_networks, d="HIM", n.cores=40)
mat_dist_mat <- ddist$HIM
rownames(mat_dist_mat) <- colnames(mat_dist_mat) <- tolower(names(infer_L1000_coexpr_networks))


### !!!!! Change cancer flag!!!!!!!!!!!
saveRDS(mat_dist_mat, file = "MOA_BRCA_distance_matrix.rds")
mat_dist_mat <- readRDS("MOA_BRCA_distance_matrix.rds")
write.table(mat_dist_mat, file = "MOA_BRCA_distance_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col=NA)
#dir.create("./Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/XML_dir")
#dir.create("./Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/SDF_dir")

#drug_smiles <- download_SDF_and_cid_from_compounds_names(compounds_name = rownames(mat_dist_mat), XML_dest = "./Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/XML_dir", SDF_dest = "./Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/SDF_dir")



#dup_smiles <- unique(L1000_all_drugs$canonical_smiles[which(L1000_all_drugs$pert_iname %in% rownames(L1000_networks_ddist)),])
# dup_smiles <- unique(L1000_all_drugs$canonical_smiles[which(L1000_all_drugs$pert_iname %in% rownames(mat_dist_mat))])
# tmp <- L1000_all_drugs[which(L1000_all_drugs$canonical_smiles %in% dup_smiles),]
# tmp <- unique(tmp)
# drugs_good_smiles <- download_SDF_and_cid_from_compounds_names(compounds_name = names(infer_L1000_coexpr_networks), XML_dest = "./Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/XML_dir", SDF_dest = "./Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/XML_dir")
# notfound <- setdiff(names(infer_L1000_coexpr_networks), drugs_good_smiles[,1])
# 
# # Here I did some manual curation to find the smiles which are not found from Angela´s function
# setwd("./BRCA_data_paired_wm/")
# notfound <- read.table("notfound_smiles.txt", header = FALSE, sep = "\t", quote = "")
# notfound <- notfound[!duplicated(notfound$V1),]
# notfound$CID=NA
# notfound$CAS=NA
# notfound <- notfound[, c("V1", "CID", "CAS", "V2")]
# drugs_good_smiles <- rbind(drugs_good_smiles, as.matrix(notfound))
# drugs_good_smiles[,1] <- sapply(drugs_good_smiles[,1], function(x) gsub(x, pattern = " ", replacement = "-"))

#L1000_networks_ddist <- L1000_networks_ddist[which(rownames(L1000_networks_ddist) %in% drugs_good_smiles[,1]), which(colnames(L1000_networks_ddist) %in% drugs_good_smiles[,1])]

drugs_good_smiles <- data.frame()
smile_vec <- c()
drug <- 1
while (drug <= length(rownames(mat_dist_mat))){
  Sys.sleep(1)
  drug_string <- rownames(mat_dist_mat)[drug]
  drug_string <- URLencode(drug_string)
  print(paste(drug_string, "-", drug))
  
  if(isFALSE(url.exists(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drug_string, "/JSON")))) {
    drug_string <- gsub(rownames(mat_dist_mat)[drug], pattern = "-", replacement = " ")
    drug_string <- URLencode(drug_string)
    
    if(isFALSE(url.exists(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drug_string, "/JSON")))) {
      smile_vec <- c(rownames(mat_dist_mat)[drug], "URL not found")
      drugs_good_smiles <- rbind(drugs_good_smiles, smile_vec)
      drug=drug+1
      next
    }
    
  }
  
  tryCatch({
    
    y <- jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drug_string, "/JSON"))
    smile <- y$PC_Compounds$props[[1]][19,]$value$sval
    smile_vec <- c(rownames(mat_dist_mat)[drug], smile)
    drugs_good_smiles <- rbind(drugs_good_smiles, smile_vec)
    
  }, error=function(e) {
    print("Error")
    drug<<-drug-1
    Sys.sleep(5)
    
  })
  
  drug=drug+1
}

colnames(drugs_good_smiles) <- c("drug", "smile")

if(length(which(drugs_good_smiles$smile=="URL not found"))!=0)
drugs_good_smiles <- drugs_good_smiles[-which(drugs_good_smiles$smile=="URL not found"),]

mat_dist_mat <- mat_dist_mat[which(rownames(mat_dist_mat) %in% drugs_good_smiles[,1]), which(colnames(mat_dist_mat) %in% drugs_good_smiles[,1])]

# flist <- head(unique(filter_by_L1000_drugs$dat.drug.molecule_name))
# flist[3] <- "tonino"
# i=1
# while(i <= length(flist)) {
#   drug_string <- flist[i]
#   drug_string <- URLencode(drug_string)
#   print(paste(drug_string, "-", i))
#   if(isFALSE(url.exists(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drug_string, "/JSON")))) {
#     print("URL doesn't exist")
#     }
#   
#   tryCatch({
#     y <- jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drug_string, "/JSON"))
#     
#   }, error=function(e){
#     print("Error")
#     i<<-i-1
#     Sys.sleep(2)
#     next
#   })
#     i=i+1
# }
# 
# 
# for(i in 1:10){
#   print(i)
#   i <<- i-1
# }

#mat_dist_mat <- mat_dist_mat[which(colnames(mat_dist_mat) %in% tolower(drugs_good_smiles[,1])), which(colnames(mat_dist_mat) %in% tolower(drugs_good_smiles[,1]))]


a <- lapply(drugs_good_smiles[,2], utf8ToInt)
SMILES_dist <- seq_distmatrix(a = a, method = "lv", nthread=20)

SMILES_dist <- as.matrix(SMILES_dist)

#colnames(SMILES_dist) <- rownames(SMILES_dist) <- colnames(mat_dist_mat)


drugs_good_smiles[,1] <- tolower(drugs_good_smiles[,1])

colnames(SMILES_dist) <- rownames(SMILES_dist) <- drugs_good_smiles[,1]

#SMILES_dist <- SMILES_dist[which(rownames(SMILES_dist) %in% rownames(mat_dist_mat)), which(colnames(SMILES_dist) %in% colnames(mat_dist_mat))]


SMILES_dist2 = as.matrix(SMILES_dist)
SMILES_dist2 = SMILES_dist2/max(SMILES_dist2)

saveRDS(SMILES_dist2, file = "SMILES_BRCA_distance_matrix.rds")


avg_drug_shortest_path <- get_avg_shortest_path(Graph = conGraph, drug_target_df = unique(filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% tolower(colnames(SMILES_dist2))),]), drugs_col = 6, targets_col = 1)
saveRDS(avg_drug_shortest_path, file = "shortestpath_BRCA_distance_matrix.rds")


avg_drug_jacc_indices <- get_avg_jacc_indices(Graph = conGraph, drug_target_df = filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% tolower(colnames(SMILES_dist2))),], drugs_col = 6, targets_col = 1, order = 1)

#avg_drug_inform_rank <- get_drug_avg_rank_from_graph(drug_target_df = filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% colnames(SMILES_dist)),], drugs_col = 6, targets_col = 1, graph_rank = graph_rank)

drug_cov_score <- get_drug_coverage_score(Graph = conGraph, drug_target_df = filter_by_L1000_drugs[which(filter_by_L1000_drugs$dat.drug.molecule_name %in% tolower(colnames(SMILES_dist2))),], drugs_col = 6, targets_col = 1, graph_rank = graph_rank)
saveRDS(drug_cov_score, file = "effect_score_BRCA_distance_matrix.rds")

finaldistrank_coverage <- Borda(list(names(sort(drug_cov_score, decreasing = TRUE)), names(sort(rowSums(avg_drug_shortest_path), decreasing = TRUE))))
finaldistrank_coverage <- finaldistrank_coverage$TopK$median

finaldistrank <- Borda(list(finaldistrank_coverage, names(sort(rowSums(as.matrix(SMILES_dist2)), decreasing = TRUE)), names(sort(rowSums(mat_dist_mat), decreasing = TRUE))))
finaldistrank <- finaldistrank$TopK$median

finaldistrank[1:20]
saveRDS(finaldistrank, file=paste0(mainDir, subDir, "finaldistrank_", project_flag, ".rds"))

png(file=paste0("additive_coverage_", project_flag, ".png"))
additive_coverage <- get_additive_graph_coverage(Graph = conGraph, drugs_rank = finaldistrank, nDrugs = 5, order = 1)
dev.off()

super_final_rank <- c()

for (drugs in 1:length(additive_coverage)) {
  final_score <-  sum(sum(mat_dist_mat[names(additive_coverage[drugs]),]), sum(SMILES_dist2[names(additive_coverage[drugs]),]), additive_coverage[drugs])
  super_final_rank <- c(super_final_rank, final_score)
}

names(super_final_rank) <- names(additive_coverage)

super_final_rank <- sort(super_final_rank, decreasing = TRUE)

save.image(paste0(mainDir, subDir, "workspace_drugmap_", project_flag, ".RData"))

### Slope analysis ###

# abline(v=3)
# 
# y <- additive_coverage
# x <- seq_along(y)
# #plot(x, y)
# 
# fit <- nls(y ~ SSgompertz(x, Asym, b2, b3), data = data.frame(x, y))
# curve(predict(fit, newdata = data.frame(x)), add = TRUE)
# 
# #assign coefficients into global environment
# list2env(as.list(coef(fit)), .GlobalEnv)
# #create function that returns the gradient
# dGomp <- deriv((y ~ Asym*exp(-b2*b3^x)), "x", func = TRUE)
# 
# #the model slopes:
# c(attr(dGomp(x), "gradient"))
# 
# 
# df <- data.frame(drug=factor(names(additive_coverage), levels = names(additive_coverage)), coverage=additive_coverage)
# 
# ggplot(data=df[order(df$drug),], aes(x=drug, y=coverage, group=1)) +
#   geom_line(color="red")+
#   geom_point() 


#save(mat_dist_mat, SMILES_dist2, avg_drug_shortest_path, drug_cov_score, file = "distance_matrices_BRCA.RData")
#save.image("./Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/BRCA_all_controls.RData")

#load("/home/MOKA/antonio/master_script_folder/Normal/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/first_test.RData")
load("./secondtest/TCGA-BRCA/workspace_drugmap_TCGA-BRCA.RData")
