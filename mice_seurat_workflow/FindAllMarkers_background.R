#Find all marker genes
library(Seurat)
mice_counts <- FindAllMarkers(mice_counts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mice_counts  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(mice_counts, file = "markersAll.csv")