#Setup libraries, installing packages is only required once, afterwards skip to library calling
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
if(!require(scTalk)){
    devtools::install_github("VCCRI/scTalk", build = TRUE, build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
    library(scTalk)
}
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

#Compare how the interaction networks are affected in the two conditions of the dataset (Donor vs PAH) - samples conditions are stored in the "disease" metadata slot
##Load your single cell object
scrna_object = readRDS("YourDirectory/scrna_object.rds")
Idents(scrna_object) = "disease"
scrna_object_split = list()
##Split the object
for (i in unique(scrna_object$disease) {
     scrna_object_split[[i]] = subset(scrna_object, idents = i)
 } 
###Interaction network
for (i in unique(scrna_object$disease) {
  dir.create(paste0("YourDirectory/", i))
  populations.use <- names(table(Idents(scrna_object_split[[i]])))
  file.lab <- "LR"
  #This will generate the table of possible edges in YourDirectory
  GenerateEdgeWeights(seurat.object = scrna_object_split[[i]],
                    file.label = file.lab,
                    species = "human",
                    populations.use = populations.use)
  edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
  edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
  head(edge.table)
  #This will generate the network paths in YourDirectory
  GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5, #Default value, possible optimisazion
                     ncores = 4)
  paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
  paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
  results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4)
  head(results.table)
}
