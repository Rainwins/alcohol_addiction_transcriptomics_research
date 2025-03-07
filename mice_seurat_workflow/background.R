library(Seurat)
library(dplyr)
load("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/temp_seurat_output/mice_clustered.RData")
#Find all marker genes
mice_counts <- FindAllMarkers(mice_counts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(mice_counts, file = "/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/
     temp_seurat_output/mice_counts_markers.RData")
mice_counts  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(mice_counts, file = "markersAll.csv")

