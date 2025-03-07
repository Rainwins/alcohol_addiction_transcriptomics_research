library(Seurat)
#mast test https://github.com/RGLab/MAST
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("MAST")
#load("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/human_seurat_workflow/suerat.merged.scsorter.umap.RData")

#Finding DEGs by condition for each celltype cluster
human.combined$celltype.condition <- paste(human.combined$sorted_cells, human.combined$DSMIV, sep = "_")
Idents(human.combined) <- "celltype.condition"

# Alcohol subtypes to compare against Control
alcohol_conditions <- c("A: Harmful Use (alcohol)", 
                        "B: Substance Dependence (alcohol)", 
                        "C: Substance Abuse (alcohol)")

list_of_celltypes <- unique(human.combined$sorted_cells)
list_of_celltype_DEGs <- list()



# Loop through each cell type
for (cell_type in list_of_celltypes) {
  
  # Initialize a sub-list for each cell type
  list_of_celltype_DEGs[[cell_type]] <- list()
  
  # Define the control group identifier
  ident_control <- paste0(cell_type, "_Control")
  
  # Check if control exists in the dataset
  if (!(ident_control %in% levels(Idents(human.combined)))) {
    message(paste("Skipping", cell_type, "because Control group is missing."))
    next
  }
  
  # Loop through each alcohol condition
  for (condition in alcohol_conditions) {
    
    # Define the condition identifier
    ident_alcohol <- paste0(cell_type, "_", condition)
    
    # Ensure both groups exist before running FindMarkers
    if (ident_alcohol %in% levels(Idents(human.combined))) {
      
      # Find differentially expressed genes using MAST
      deg_results <- FindMarkers(human.combined, 
                                 ident.1 = ident_alcohol, 
                                 ident.2 = ident_control, 
                                 verbose = TRUE, 
                                 test.use = 'MAST')
      
      # Store results in the nested list
      list_of_celltype_DEGs[[cell_type]][[condition]] <- deg_results
      
    } else {
      message(paste("Skipping", condition, "for", cell_type, "because it is missing in the dataset."))
    }
  }
}

save(list_of_celltype_DEGs, file = "aud_subtypes_list_of_celltype_DEGs.RData")
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


