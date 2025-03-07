#this file contains the R script for SCcatch to label cell types on the mice data
install.packages("scCATCH")
library(scCATCH)
library(Seurat)

mice_counts_scCatch <- createscCATCH(data = mice_counts[['SCT']]@data, cluster = as.character(Idents(mice_counts)))

# find marker gene for each cluster
mice_counts_scCatch <- findmarkergene(mice_counts_scCatch, species = 'Mouse',marker = cellmatch,
                                      tissue = 'Brain', cancer = "Normal")

# find cell type for each cluster
mice_counts_scCatch <- findcelltype(mice_counts_scCatch)

save(mice_counts_scCatch, file = "./temp_seurat_output/mice_counts_scCatch")

#visualize the SCcatch results in a UMAP
load("temp_seurat_output/mice_clustered.RData")
mice_counts$sorted_cells = mice_counts_scCatch$celltype
mice_counts <- RunUMAP(object = mice_counts, dims = 1:19, n.neighbors = 50, min.dist = 0.2)
pdf("UMAP_scSorter.pdf")
DimPlot(mice_counts, group.by = c("sorted_cells"), raster = FALSE)
dev.off()

pdf("UMAP_scSorter_splitCond.pdf")
DimPlot(mice_counts, group.by = "sorted_cells", split.by = "Ident", raster = FALSE)
dev.off()

