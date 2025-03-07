
library(Seurat)

# making seurat object
mice_counts <- Read10X(data.dir = "/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/GEO_downloads_mice")
mice_counts <- CreateSeuratObject(counts = mice_counts)
#add metadata based on condition group for the mice data
meta_data <- read_csv("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/mice_data/GSE233763_meta_CIE_data.csv")

# check that the metadata and the seurat have the same alignment
meta_data <- as.data.frame(meta_data)
rownames(meta_data) <- meta_data$...1
if (all(rownames(meta_data) == rownames(mice_counts@meta.data))) {
  print("Row names are identical")
} else {
  print("Row names are not identical")
}

#add meta condition group after confirming that it aligns
mice_counts$Ident <- meta_data$group
save(mice_counts, file = "mice_seurat_counts.RData") # seurat object

# capitalize gene names
mice_counts@assays[["RNA"]]@counts@Dimnames[[1]] <- toupper(mice_counts@assays[["RNA"]]@counts@Dimnames[[1]])
mice_counts@assays[["RNA"]]@data@Dimnames[[1]] <- toupper(mice_counts@assays[["RNA"]]@data@Dimnames[[1]])

# Normalize Data Using SCTransform and Regress out Genes Related to Cell Cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
mice_counts <- CellCycleScoring(mice_counts, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mice_counts$CC.Difference <- mice_counts$S.Score - mice_counts$G2M.Score
mice_counts <- SCTransform(mice_counts, vars.to.regress = "CC.Difference")

save(mice_counts, file = "normalized_mice_counts.RData")

#load it in only if you didn't have the normalized data already on your environment
mice_counts <- load("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/temp_seurat_output/normalized_mice_counts.RData")

# Identify Highly Variable Features and Generate Plots
mice_counts <- FindVariableFeatures(mice_counts, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mice_counts), 10)
pdf("variableFeaturesPlot.pdf")
plot1 <- VariableFeaturePlot(mice_counts)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

save(mice_counts, file = "mice_highly_variable_features.RData")

####RECEIVED ERROR IGNORE FOR NOW#######
mice_counts <- FindVariableFeatures(mice_counts, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mice_counts), 10)
pdf("variableFeaturesPlot.pdf")
plot1 <- VariableFeaturePlot(mice_counts)
plot1 <- plot1 + scale_x_continuous(limits = c( )) + theme(axis.text.x = element_text(size = 7)) # to scale plot
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
########################################


# scale data
all.genes <- rownames(mice_counts)
mice_counts <- ScaleData(mice_counts, features = all.genes)

save(mice_counts, file = "mice_scaled.RData")

#load("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/temp_seurat_output/mice_scaled.RData")
# pca plots
mice_counts <- RunPCA(mice_counts, features = VariableFeatures(object = mice_counts))
pdf("findDimensions.pdf")
print(mice_counts[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mice_counts, dims = 1:2, reduction = "pca")
DimPlot(mice_counts, reduction = "pca")

#Generate Principal Component Heatmaps through PC_30
DimHeatmap(mice_counts, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mice_counts, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(mice_counts, dims = 16:30, cells = 500, balanced = TRUE)

#Create Elbow Plot of Principal Components through PC_30
ElbowPlot(mice_counts, ndims = 30, reduction = "pca")
dev.off()

save(mice_counts, file = "mice_pca.RData")

#Cluster the Cells for the selected number of Principal Components 
mice_counts <- FindNeighbors(mice_counts, dims = 1:19)
mice_counts <- FindClusters(mice_counts, resolution = 0.5)

save(mice_counts, file = "mice_clustered.RData")

#View the Cluster IDs of the First 5 Cells
head(Idents(mice_counts), 5)

#if in a new session, reload it and name it mice_clustered


#change the active identifier to Ident/condition group
mice_counts <- SetIdent(mice_counts, value = "Ident")

#run UMAP on first 9 PCs
#the default for the n_neighbors = 30 and min_dist = 0.1
#increase n_neighbors lead to larger cluster and increasing min_distance lead to fewer cluster
pdf("UMAP.pdf")
mice_counts <- RunUMAP(object = mice_counts, dims = 1:19, n.neighbors = 50, min.dist = 0.2)
# Plot results
DimPlot(object = mice_counts, reduction = 'umap', raster = FALSE, pt.size = 1)
dev.off()

#Find all marker genes
mice_counts <- FindAllMarkers(mice_counts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mice_counts  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(mice_counts, file = "markersAll.csv")

