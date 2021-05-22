##### HOMEWORK AT THE VERY BOTTOM #######

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/jakub/Desktop/filtered_matrices_mex/hg19/")
pbmc.data

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc68k", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)

#pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
#cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
#head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

###################################
library(dplyr)
library(Seurat)
library(patchwork)

# pbmc.counts.homework <- Read10X(data.dir = "/Users/jakub/Desktop/filtered_matrices_mex/hg19/")
# pbmc_homework <- CreateSeuratObject(counts = pbmc.counts.homework, project = "pbmc68k", min.cells = 3, min.features = 200)
# pbmc_homework <- NormalizeData(object = pbmc_homework)
# pbmc_homework <- FindVariableFeatures(object = pbmc_homework)
# pbmc_homework <- ScaleData(object = pbmc_homework)
# pbmc_homework <- RunPCA(object = pbmc_homework)
# pbmc_homework <- FindNeighbors(object = pbmc_homework)
# pbmc_homework <- FindClusters(object = pbmc_homework)
# pbmc_homework <- RunTSNE(object = pbmc_homework)
# DimPlot(object = pbmc_homework, reduction = "tsne")

# # another attempt
# pbmc.data <- Read10X(data.dir = "/Users/jakub/Desktop/filtered_matrices_mex/hg19/")
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc68k", min.cells = 3, min.features = 200)
# pbmc <- NormalizeData(pbmc)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# pbmc <- ScaleData(pbmc)
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# pbmc <- FindNeighbors(pbmc, dims = 1:5)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# pmbc <- kmeans(pbmc, centers = 10, iter.max = 10)
# pbmc <- RunTSNE(object = pbmc, dims = 10)
# DimPlot(object = pbmc)
# above function plots the plot

###################
library(Rtsne)
library(data.table)
library(ggplot2)

pbmc.data <- Read10X(data.dir = "/Users/jakub/Desktop/filtered_matrices_mex/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc68k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc.embedet <- Embeddings(object = pbmc, reduction = "pca")[,1:20]
k.mean <- kmeans(pbmc.embedet, 10)
k.mean.clusters <- k.mean$cluster
set.seed(1)
tsne_out <- Rtsne(pbmc.embedet,pca=F,perplexity=30)
tsne_out_pos = data.table(tsne_out$Y)
tsne_out_pos$cluster <- k.mean$cluster
ggplot(tsne_out_pos) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))  + labs(color = "cluster")

# Second part of the homework
cluster.1 <- pbmc.embedet[k.mean.clusters == 1,]
cluster.2 <- pbmc.embedet[k.mean.clusters == 2,]
cluster.3 <- pbmc.embedet[k.mean.clusters == 3,]
cluster.4 <- pbmc.embedet[k.mean.clusters == 4,]
cluster.5 <- pbmc.embedet[k.mean.clusters == 5,]
cluster.6 <- pbmc.embedet[k.mean.clusters == 6,]
cluster.7 <- pbmc.embedet[k.mean.clusters == 7,]
cluster.8 <- pbmc.embedet[k.mean.clusters == 8,]
cluster.9 <- pbmc.embedet[k.mean.clusters == 9,]
cluster.10 <- pbmc.embedet[k.mean.clusters == 10,]

library(factoextra)
library(Rfast)

###cluster 1
cluster.1.kmean <- kmeans(cluster.1, 5)
clusters.cluster.1 <- cluster.1.kmean$cluster
cluster.1.tsne <- Rtsne(cluster.1, pca = F, perplexity = 30)
cluster.1.tsne.dt <- data.table(cluster.1.tsne$Y) 
cluster.1.tsne.dt$cluster <- cluster.1.kmean$cluster
ggplot(cluster.1.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 2
cluster.2.kmean <- kmeans(cluster.2, 6)
clusters.cluster.2 <- cluster.2.kmean$cluster
cluster.2.tsne <- Rtsne(cluster.2, pca = F, perplexity = 30)
cluster.2.tsne.dt <- data.table(cluster.2.tsne$Y) 
cluster.2.tsne.dt$cluster <- cluster.2.kmean$cluster
ggplot(cluster.2.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 3
cluster.3.kmean <- kmeans(cluster.3, 5)
clusters.cluster.3 <- cluster.3.kmean$cluster
cluster.3.tsne <- Rtsne(cluster.3, pca = F, perplexity = 30)
cluster.3.tsne.dt <- data.table(cluster.3.tsne$Y) 
cluster.3.tsne.dt$cluster <- cluster.3.kmean$cluster
ggplot(cluster.3.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 4
cluster.4.kmean <- kmeans(cluster.4, 5)
clusters.cluster.4 <- cluster.4.kmean$cluster
cluster.4.tsne <- Rtsne(cluster.4, pca = F, perplexity = 30)
cluster.4.tsne.dt <- data.table(cluster.4.tsne$Y) 
cluster.4.tsne.dt$cluster <- cluster.4.kmean$cluster
ggplot(cluster.4.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 5
cluster.5.kmean <- kmeans(cluster.5, 5)
clusters.cluster.5 <- cluster.5.kmean$cluster
cluster.5.tsne <- Rtsne(cluster.5, pca = F, perplexity = 30)
cluster.5.tsne.dt <- data.table(cluster.5.tsne$Y) 
cluster.5.tsne.dt$cluster <- cluster.5.kmean$cluster
ggplot(cluster.5.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 6
cluster.6.kmean <- kmeans(cluster.6, 5)
clusters.cluster.6 <- cluster.6.kmean$cluster
cluster.6.tsne <- Rtsne(cluster.6, pca = F, perplexity = 30)
cluster.6.tsne.dt <- data.table(cluster.6.tsne$Y) 
cluster.6.tsne.dt$cluster <- cluster.6.kmean$cluster
ggplot(cluster.6.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 7
cluster.7.kmean <- kmeans(cluster.7, 5)
clusters.cluster.7 <- cluster.7.kmean$cluster
cluster.7.tsne <- Rtsne(cluster.7, pca = F, perplexity = 30)
cluster.7.tsne.dt <- data.table(cluster.7.tsne$Y) 
cluster.7.tsne.dt$cluster <- cluster.7.kmean$cluster
ggplot(cluster.7.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 8
cluster.8.kmean <- kmeans(cluster.8, 5)
clusters.cluster.8 <- cluster.8.kmean$cluster
cluster.8.tsne <- Rtsne(cluster.8, pca = F, perplexity = 30)
cluster.8.tsne.dt <- data.table(cluster.8.tsne$Y) 
cluster.8.tsne.dt$cluster <- cluster.8.kmean$cluster
ggplot(cluster.8.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 9
cluster.9.kmean <- kmeans(cluster.9, 5)
clusters.cluster.9 <- cluster.9.kmean$cluster
cluster.9.tsne <- Rtsne(cluster.9, pca = F, perplexity = 30)
cluster.9.tsne.dt <- data.table(cluster.9.tsne$Y) 
cluster.9.tsne.dt$cluster <- cluster.9.kmean$cluster
ggplot(cluster.9.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

###cluster 10
cluster.10.kmean <- kmeans(cluster.10, 5)
clusters.cluster.10 <- cluster.10.kmean$cluster
cluster.10.tsne <- Rtsne(cluster.10, pca = F, perplexity = 30)
cluster.10.tsne.dt <- data.table(cluster.10.tsne$Y) 
cluster.10.tsne.dt$cluster <- cluster.10.kmean$cluster
ggplot(cluster.10.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))

# I am attaching 3 plots because R is shutting down after I try to save
# more plots. The script generates all 10 of them
