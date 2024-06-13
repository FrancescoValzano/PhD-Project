import scvelo as scv 

scrna_object = scv.read("YourDirectory/scrna_object.h5ad")
scv.pp.filter_genes_dispersion(scrna_object)
scv.tl.velocity(scrna_object, mode = "stochastic")
scv.tl.recover_dynamics(scrna_object)
scv.tl.velocity_graph(scrna_object)

scv.pl.velocity_embedding(scrna_object, basis = "umap", arrow_size = 1.5, arrow_length=5, color= "black", size =25, density = 0.5)
scv.pl.velocity_embedding_stream(scrna_object, basis = "umap", arrow_size = 1.5, color= "Clusters", size =25) 
scv.pl.velocity_embedding_grid(scrna_object, basis = "umap", arrow_size = 4, arrow_length=2, color= "Clusters", size =35, density = 0.5, alpha = 1, dpi= 150)
scv.pl.scatter(scrna_object, color = "velocity_pseudotime", cmap = "gnuplot")
