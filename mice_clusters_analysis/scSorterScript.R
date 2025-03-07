# SC Sorter
# install.packages("scSorter")
library(scSorter)
library(Seurat)
library(dplyr)
library(readxl)

load("temp_seurat_output/mice_highly_variable_features.RData")

# pre processing of Seurat object
topgenes <- head(VariableFeatures(mice_counts), 2000)
expr <- GetAssayData(mice_counts)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

# creating annotation marker file
anno <- data.frame(
  Marker = c(
    # Astrocyte markers
    "ALDH1L1", "ALDOC", "AQP4", "GJA1", "SLC1A3", "SLC1A2",
    "LGALS3", "GAP43", "GFAP", "GLUL", "HES1", "NDRG2",
    "NOTCH1", "PPP1R13B", "S100B", "SOX9", "BIRC5",
    
    # Neuron general markers
    "CALB1", "DCX", "ENO2", "MAP2", "RBFOX3", "NeuroD1",
    "NFL", "PSD-93", "PSD-95", "TUBB3", "UCHL1",
    
    # Dopaminergic neuronal markers
    "FSHB", "KCNJ6", "NR4A2", "SLC6A2", "SLC6A3", "TH",
    
    # Gabaergic neuronal markers
    "GABBR1", "GABBR2", "GAD1", "GAD2", "SLC6A1", "SLC32A1",
    
    # Glutamatergic Neuronal Markers
    "GRIN1", "GRIN2B", "SLC17A7", "SLC17A6",
    
    # Glycinergic Neuronal Markers
    "SLC6A5", "SLC32A1",
    
    # Serotonergic Neuronal Markers
    "DDC", "POU3F1", "SLC6A4", "TPH1", "SLC18A2",
    
    # Cholinergic Neuronal Markers
    "ACHE", "CHAT", "SLC18A3",
    
    # Oligodendrocyte precursor markers
    "A2B5", "CSPG4", "NKX2-2", "OLIG2", "PDGFRA", "SOX10",
    
    # Mature, non-myelinating Oligodendrocyte markers
    "APC", "CNP1", "O1M1", "O4", "PLP1", "MOG",
    
    # Mature myelinating Oligodendrocyte markers
    "APC", "CSPG4", "CNP1", "MAG", "MBP", "MOG",
    "PLP1", "NOGO", "OMGP", "PPP3CC",
    
    # Microglia markers
    "ITGAM", "CD14", "FCGR3", "CD40", "PTPRC", "CSF1R",
    "CX3CR1", "TMEM119", "ADGRE1", "FCER1G", "FCRLS",
    "SIRPA", "Siglec", "SLC2A5", "P2RY12",
    
    # Intercellular microglia markers
    "SPI1", "AIF1", "HEXB", "VIM", "FTL", "Sall1"
  ),
  Type = c(
    # Astrocyte markers
    rep("Astrocyte", 17),
    
    # Neuron general markers
    rep("Neuron", 11),
    
    # Dopaminergic neuronal markers
    rep("Dopaminergic Neuron", 6),
    
    # Gabaergic neuronal markers
    rep("Gabaergic Neuron", 6),
    
    # Glutamatergic Neuronal Markers
    rep("Glutamatergic Neuron", 4),
    
    # Glycinergic Neuronal Markers
    rep("Glycinergic Neuron", 2),
    
    # Serotonergic Neuronal Markers
    rep("Serotonergic Neuron", 5),
    
    # Cholinergic Neuronal Markers
    rep("Cholinergic Neuron", 3),
    
    # Oligodendrocyte precursor markers
    rep("Oligodendrocyte Precursor", 6),
    
    # Mature, non-myelinating Oligodendrocyte markers
    rep("Mature Non-myelinating Oligodendrocyte", 6),
    
    # Mature myelinating Oligodendrocyte markers
    rep("Mature Myelinating Oligodendrocyte", 10),
    
    # Microglia markers
    rep("Microglia", 15),
    
    # Intercellular microglia markers
    rep("Intercellular Microglia", 6)
  )
)

picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# now running SC Sorter
debug(scSorter)
rts_2 <- scSorter(expr, anno)

print(table(rts$Pred_Type))

# graphing cell types

load("mice_seurat_output/mice_clustered.RData")
mice_counts$sorted_cells = rts$Pred_Type
mice_counts <- RunUMAP(object = mice_counts, dims = 1:19, n.neighbors = 50, min.dist = 0.2)

save(mice_counts, file = "mice_umap.RData")

load("mice_seurat_output/mice_umap.RData")

pdf("UMAP_scSorter.pdf")
DimPlot(mice_counts, group.by = c("sorted_cells"), raster = FALSE)
dev.off()

pdf("UMAP_scSorter_splitCond.pdf")
DimPlot(mice_counts, group.by = "sorted_cells", split.by = "Ident", raster = FALSE)
dev.off()

# permutation plot
library("scProportionTest")
prop_test <- sc_utils(mice_counts)
prop_test <- permutation_test(prop_test, cluster_identity = "sorted_cells",
  sample_1 = "vapor", sample_2 = "control", sample_identity = "Ident")

pdf("proportion_cell_types_mice.pdf")
permutation_plot(prop_test)
dev.off()

# pie chart
percentages <- table(mice_counts$sorted_cells) %>% as.data.frame()
colnames(percentages) <- c("CellType","NumberOfCells")
# Calculate the percentages based on the total number of cells
percentages$Percentage <- 100 * percentages$NumberOfCells / sum(percentages$NumberOfCells)

pie(percentages$Percentage, labels = rep("", length(percentages$CellType)), 
    main = "Distribution of Cell Types", 
    col = brewer.pal(length(percentages$CellType), "Set3"))

# Add a legend
legend("topright", 
       legend = paste(percentages$CellType, "(", round(percentages$Percentage, 1), "%)", sep=""), 
       fill = brewer.pal(length(percentages$CellType), "Set3"), 
       title = "Cell Types", 
       cex = 0.8)


percentages_mice <- percentages

load("human_seurat_workflow/rts_2_scSorter_human.RData")
percentages_human <- table(rts_2$Pred_Type) %>% as.data.frame()
colnames(percentages_human) <- c("CellType","NumberOfCells")
percentages_human$Percentage <- 100 * percentages_human$NumberOfCells / sum(percentages_human$NumberOfCells)

comparison_data <- merge(percentages_mice, percentages_human, by="CellType", suffixes = c("_Mice", "_Humans"))

# Create a grouped bar chart
library(ggplot2)
ggplot(comparison_data, aes(x = CellType)) +
  geom_bar(aes(y = Percentage_Mice, fill = "Mice"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = Percentage_Humans, fill = "Humans"), stat = "identity", position = "dodge") +
  labs(title = "Comparison of Cell Type Percentages: Mice vs Humans", 
       x = "Cell Type", 
       y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Mice" = "blue", "Humans" = "red"))

# Reshape data to long format for ggplot
library(tidyr)
comparison_data_long <- comparison_data %>%
  gather(key = "Group", value = "Percentage", Percentage_Mice, Percentage_Humans) %>%
  mutate(Group = ifelse(Group == "Percentage_Mice", "Mice", "Humans"))

comparison_data_long$Percentage <- round(comparison_data_long$Percentage, 1)

# Create a grouped bar chart with percentage labels
ggplot(comparison_data_long, aes(x = CellType, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # Position bars side by side
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +  # Add labels with percentages
  labs(title = "Comparison of Cell Type Percentages: Humans vs Mice", 
       x = "Cell Type", 
       y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_fill_manual(values = c("Mice" = "coral", "Humans" = "coral3"))  # Set custom colors


# running scSorter with different set of markers
mice_markers <- data.frame(read_excel("mice_markers.xlsx"))
mice_markers_2 <- mice_markers %>% pivot_longer(cols = everything(), names_to = 'Type', values_to = 'Marker') %>%
  filter(!is.na(`Marker`)) %>% select(`Marker`, `Type`) %>% data.frame()

# todo taking into account how the same gene can be a marker for multiple cell types
picked_genes_2 = unique(c(mice_markers_2$Marker, topgenes))
expr_2 = expr[rownames(expr) %in% picked_genes_2, ]
rts_2 <- scSorter(expr_2, mice_markers_2) 
 
# cell type percentages n
(table(rts_2$Pred_Type) / length(rts_2$Pred_Type)) * 100

# rerun find markers within cell types, looking at cpa6 in each cell type
load("mice_umap.RData")

# should I run find markers on all cell types or just the ones with significantly 
# different proportions, should I run separately on alcohol and control cell types
# and compare biomarkers or all together

# should the list of cell type markers used for humans be different?

#todo look into qs library

# Todo: pie chart, run find all markers on cell types, run find markers deg - look at website
load("mice_umap.RData")
#mice_counts$seurat_clusters = mice_counts$sorted_cells
Idents(object = mice_counts) <- "sorted_cells"
markers_cell_types <- FindAllMarkers(mice_counts, min.pct = 0.1, logfc.threshold = 0.5)
save(markers_cell_types, file = "markers_cell_types.RData")
load("markers_cell_types.RData")

# deg analysis
library(DESeq2)
mice_counts$celltype_control <- paste(mice_counts$sorted_cells, mice_counts$Ident, sep = "_")
Idents(mice_counts) <- "celltype_control"
debug(FindMarkers)
deg_test <- FindMarkers(mice_counts, ident.1 = "Astrocyte_control", ident.2 = "Astrocyte_vapor", test.use = "DESeq2")


# other plots
top10 <- markers_cell_types %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmcTSNE, features = top10$gene, size = 3) + theme(axis.text.y = element_text(size = 5)) + NoLegend()

markers_cell_type_counts <- table(anno$Type) %>% as.data.frame()
colnames(markers_cell_type_counts) <- c("Cell Type", "Num. Markers")
markers_cell_type_counts$Genes <- anno$Marker %>% group_by(anno$Type)

library(dplyr)

genes_per_cell_type <- anno %>%
  group_by(Type) %>%
  summarise(Genes = paste(Marker, collapse = ", "), .groups = "drop")

markers_cell_type_counts$Genes <- genes_per_cell_type$Genes
colnames(markers_cell_type_counts) <- c("Cell Type", "Num. Markers", "Genes")

# Now, markers_cell_type_counts will have the Genes column populated correctly



