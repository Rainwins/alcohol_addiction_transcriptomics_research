#looking at and making markers 
library(devtools)
devtools::install_github('immunogenomics/presto')
library(Seurat)
library(tidyverse)
#

load("suerat.merged.scsorter.umap.RData")

Idents(object = human.combined) <- "sorted_cells"
markers_cell_types <- FindAllMarkers(human.combined, min.pct = 0.1, logfc.threshold = 0.5)
save(markers_cell_types, file = "markers_cell_types.findallamrkers.RData")


# human.markers <- FindAllMarkers(human.combined, only.pos = TRUE)
# human.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 2) -> filtered_human.markers
# 
# write.csv(filtered_human.markers, "human_markers.csv")
