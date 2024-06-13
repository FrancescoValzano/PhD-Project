if(!require(Seurat)){
    remotes::install_version("Seurat", "4.0.3")
    library(Seurat)
}

scrna_object = readRDS("YourDirectory/scrna_object.rds")
#This will create a .h5seurat file in your directory
SaveH5Seurat(scrna_object, filename = "scrna_object.h5Seurat")
#This will convert the .h5seurat file in a .h5ad file
Convert("scrna_object.h5Seurat", dest = "h5ad")
