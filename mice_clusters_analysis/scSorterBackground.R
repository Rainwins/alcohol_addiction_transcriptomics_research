library(scSorter)
library(Seurat)
library(dplyr)
library(readxl)
library(tidyr)

load("temp_seurat_output/mice_clustered.RData")
topgenes <- head(VariableFeatures(mice_counts), 2000)
expr <- GetAssayData(mice_counts)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

mice_markers <- data.frame(read_excel("mice_markers.xlsx"))
mice_markers_2 <- mice_markers %>% pivot_longer(cols = everything(), names_to = 'Type', values_to = 'Marker') %>%
  filter(!is.na(`Marker`)) %>% select(`Marker`, `Type`) %>% data.frame()

picked_genes_2 = unique(c(mice_markers_2$Marker, topgenes))
expr_2 = expr[rownames(expr) %in% picked_genes_2, ]
rts_2 <- scSorter(expr_2, mice_markers_2)

# graphing cell types
mice_counts$sorted_cells = rts_2$Pred_Type
mice_counts <- RunUMAP(object = mice_counts, dims = 1:19, n.neighbors = 50, min.dist = 0.2)

save(mice_counts, file = "mice_umap_paper_markers.RData")

pdf("UMAP_scSorter_paper.pdf")
DimPlot(mice_counts, group.by = c("sorted_cells"), raster = FALSE)
dev.off()

pdf("UMAP_scSorter_splitCond_paper.pdf")
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

(table(rts_2$Pred_Type) / length(rts_2$Pred_Type)) * 100


