library(Seurat)
#library(SeuratData)
load("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/human_seurat_workflow/suerat.merged.scsorter.umap.RData")

#Finding DEGs by condition for each celltype cluster
human.combined$celltype.condition <- paste(human.combined$sorted_cells, human.combined$orig.ident, sep = "_")
Idents(human.combined) <- "celltype.condition"

list_of_celltypes <- unique(human.combined$sorted_cells)

list_of_celltype_DEGs <- c()

#for each cell type:
for (cell_type in list_of_celltypes){
  
  # Create identifiers for each cell type and condition
  ident_aud <- paste0(cell_type, "_AUD")
  ident_control <- paste0(cell_type, "_control")
  
  # Find markers between conditions for each cell type
  mono.de <- FindMarkers(human.combined, ident.1 = ident_aud, ident.2 = ident_control, verbose = FALSE)
  
  # Store results in the list with the cell type as the key
  list_of_celltype_DEGs[[cell_type]] <- mono.de
  
}

save(list_of_celltype_DEGs, file = "list_of_celltype_DEGs.RData")
# Display the first few results for each cell type
lapply(list_of_celltype_DEGs, head, n = 10)


#===More stringent cutoff===

filtered_degs_human <- c()

for (i in 1:length(list_of_celltype_DEGs)) {
  temp_DEGs <- list_of_celltype_DEGs[[i]]
  temp_DEGs <- temp_DEGs[temp_DEGs$p_val_adj < 0.05, ]
  temp_DEGs <- temp_DEGs[abs(temp_DEGs$avg_log2FC) > 0.8, ]
  filtered_degs_human[[ names(list_of_celltype_DEGs)[[i]]]] <- temp_DEGs
}
rm(temp_DEGs)
rm(i)

#===saving signicant celltypes as csvs===#

write.csv(list_of_celltype_DEGs$Astrocyte, "human_astrocytes_DEs.csv")
write.csv(list_of_celltype_DEGs$`Mature Non-myelinating Oligodendrocyte`, "Mature_Non-myelinating_Oligodendrocyte_DEs.csv")


