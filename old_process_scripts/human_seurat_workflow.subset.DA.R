library(tidyr)
library(Seurat)
library(SeuratObject)
library(future)

#setting new global max size for memory to 50GiB:
plan(multisession, workers = 4)
options(future.globals.maxSize = 200 * 1024^3) 

#Loading in the Seurat Object:
load("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/human_seurat_workflow/seurat.merged_v4.Rdata")

#data normalization

#human.combined <- NormalizeData(human.combined)
#human.combined <- ScaleData(human.combined)

#getting the metadata into frame:
#human_meta_data <- human.combined@meta.data

#human_meta_data$cell_names <- rownames(human_meta_data)

# Split the cell_names column into 'sample_id' and 'cell_barcode'
#human_meta_data <- separate(human_meta_data, cell_names, into = c("sample_id", "cell_barcode"), sep = "_")

#human.combined@meta.data <- human_meta_data

#changing rownames in meta.data

cell_meta_rownames <- Cells(human.combined) #character class type

#assign this vector as rownames
rownames(human.combined@meta.data) <- cell_meta_rownames

####NORMALIZATION

#filtering out missing genes:
# Step 1: Check for feature names in the Seurat object
#head(rownames(human.combined))

# Step 2: Convert gene lists to uppercase if needed (adjust as necessary)
#s.genes <- toupper(cc.genes$s.genes)
#g2m.genes <- toupper(cc.genes$g2m.genes)

#s.genes <- s.genes[s.genes %in% rownames(human.combined)]
#g2m.genes <- g2m.genes[g2m.genes %in% rownames(human.combined)]

#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
#human.combined <- CellCycleScoring(human.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#human.combined$CC.Difference <- human.combined$S.Score - human.combined$G2M.Score

human.combined <- SCTransform(human.combined, vars.to.regress = "CC.Difference", ncells = 1000, variable.features.n = 1500)

save.image("seurat.merged.SCTransform.Rdata")

#load("seurat.merged.SCTransform.Rdata")

#####IDENTIFY HIGHLY VARIABLE FEATURES:
human.combined <- FindVariableFeatures(human.combined, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(human.combined), 10)
pdf("variableFeaturesPlot.pdf")
plot1 <- VariableFeaturePlot(human.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

###Scaling the data:
all.genes <- rownames(human.combined)
human.combined <- ScaleData(human.combined, features = VariableFeatures(object = human.combined))

#PCA for dimensionality reduction:
#Perform PCA for Dimensionality Reduction
human.combined <- RunPCA(human.combined, features = VariableFeatures(object = human.combined))

#Print and Examine PCA Results in Multiple Formats
pdf("findDimensions.pdf")
print(human.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(human.combined, dims = 1:2, reduction = "pca")
DimPlot(human.combined, reduction = "pca")

#Generate Principal Component Heatmaps through PC_30
DimHeatmap(human.combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(human.combined, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(human.combined, dims = 16:30, cells = 500, balanced = TRUE)

#Create Elbow Plot of Principal Components through PC_30
ElbowPlot(human.combined, ndims = 30, reduction = "pca")
dev.off()

#Cluster the Cells for the selected numeber of Principal Components
human.combined <- FindNeighbors(human.combined, dims = 1:14)
human.combined <- FindClusters(human.combined, resolution = 0.5)

#View the Cluster IDs of the First 5 Cells
head(Idents(human.combined), 5)

#Visualizing using TSNE
pdf("tsne.pdf")
human.combined <- RunTSNE(human.combined, dims = 1:14)
DimPlot(human.combined, reduction = "tsne")
dev.off()

save.image("seurat.merged.SCTransform.Clustered.RData")
