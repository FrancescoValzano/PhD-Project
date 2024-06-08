library(Seurat)
library(scTalk)

setwd("C:/Users/FRANCESCO/Desktop/R/Pulmonary Hypertension/Smooth_Muscle_Cells_SC/SmoothMuscleCells/LR-SMCs/Donor")

populations.use <- names(table(Idents(scrna.smc.wo4.D)))
file.lab <- "LR"

GenerateEdgeWeights(seurat.object = scrna.smc.D,
                    file.label = file.lab,
                    species = "human",
                    populations.use = populations.use)
edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)
GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)
paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)
results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4)
head(results.table)


setwd("C:/Users/FRANCESCO/Desktop/R/Pulmonary Hypertension/Smooth_Muscle_Cells_SC/SmoothMuscleCells/LR-SMCs/IPAH")

populations.use <- names(table(Idents(scrna.smc.wo4.ipah)))
file.lab <- "LR"

GenerateEdgeWeights(seurat.object = scrna.smc.IPAH,
                    file.label = file.lab,
                    species = "human",
                    populations.use = populations.use)
edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)
GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)
paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)
results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4)
head(results.table)
