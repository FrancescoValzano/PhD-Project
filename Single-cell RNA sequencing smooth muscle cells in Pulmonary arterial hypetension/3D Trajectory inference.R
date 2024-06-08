if(!require(remotes)){
    install.packages("remotes")
    library(remotes)
}
if(!require(Seurat)){
    remotes::install_version("Seurat", "4.0.3")
    library(Seurat)
}
if(!require(devtools)){
    install.packages("devtools")
    library(devtools)
}
if (!requireNamespace("BiocManager")){
  install.packages("BiocManager")
  BiocManager::install(version = "3.14")
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
}
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
library(monocle3)
if(!require(ggplot2)){
   install.packages("ggplot2")    
   library(ggplot2)
}
if(!require(reshape2)){
   install.packages("reshape2")    
   library(reshape2)
}
if(!require(dplyr)){
   install.packages("dplyr")    
   library(dplyr)
}
if(!require(readxl)){
   install.packages("readxl")    
   library(readxl)
}

#Recreate cds input monocle3 file from Seurat object
scrna_object = readRDS("YourDirectory/scrna_object.rds")
scrna_subset = subset(scrna_object, idents = c("Fibroblasts", "SMC1", "SMC2"))
all.genes <- rownames(scrna_subset)
scrna_subset <- ScaleData(scrna_subset, features = all.genes)
scrna_subset <- RunPCA(scrna_subset, features = VariableFeatures(object = scrna_subset))
scrna_subset <- FindNeighbors(scrna_subset, dims = 1:10)
i = seq(0.1,1,by=0.1)
scrna_subset <- FindClusters(scrna_subset, resolution = i)
gene_annotation <- as.data.frame(rownames(scrna_subset@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(scrna_subset@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
cell_metadata <- as.data.frame(scrna_subset@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = scrna_subset@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
New_matrix <- scrna_subset@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(scrna_subset@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
list_cluster <- scrna_subset$Cluster_IDs
names(list_cluster) <- scrna_subset@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-scrna_subset@reductions[["umap"]]@cell.embeddings
cds_from_seurat@preprocess_aux$gene_loadings <- scrna_subset@reductions[["pca"]]@feature.loadings
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F, learn_graph_control = list(ncenter = 1000))
plot_cells_3d(cds_from_seurat)
#Select starting nodes - Fibroblast and SMC1
cds_from_seurat = order_cells(cds_from_seurat)
plot_cells(cds_from_seurat, 
           color_cells_by = 'pseudotime',
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=5, labels_per_group = F) 

