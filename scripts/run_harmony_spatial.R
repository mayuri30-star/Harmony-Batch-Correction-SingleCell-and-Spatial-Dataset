#Spatial
library(Seurat)
library(harmony)
library(ggplot2)

# LOAD DATA 
sample1_counts <- Read10X_h5("GSM6319696_Visium_2_raw_feature_bc_matrix.h5")
sample1 <- CreateSeuratObject(counts = sample1_counts)
sample1$batch <- "visium_2"

sample2_counts <- Read10X_h5("GSM6319697_Visium_3_raw_feature_bc_matrix.h5")
sample2 <- CreateSeuratObject(counts = sample2_counts)
sample2$batch <- "visium_3"

sample3_counts <- Read10X_h5("GSM6319698_Visium_4_raw_feature_bc_matrix.h5")
sample3 <- CreateSeuratObject(counts = sample3_counts)
sample3$batch <- "visium_4"

# MERGE
combined <- merge(sample1, y = list(sample2, sample3))

# PREPROCESSING
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

# ELBOW PLOT
ElbowPlot(combined)
ggsave("elbow_plot.png")

# BEFORE HARMONY
combined <- RunUMAP(combined, dims = 1:30)

p1 <- DimPlot(combined, group.by = "batch") + ggtitle("Before Harmony - Batch")
ggsave("umap_before_harmony.png", plot = p1)

# ADD THIS (BEFORE CLUSTERING)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined)

p_before_cluster <- DimPlot(combined, group.by = "seurat_clusters") + 
  ggtitle("Before Harmony - Clusters")
ggsave("clusters_before_harmony.png", plot = p_before_cluster)

# RUN HARMONY
combined <- RunHarmony(combined, group.by.vars = "batch")

# AFTER HARMONY
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

p2 <- DimPlot(combined, group.by = "batch") + ggtitle("After Harmony - Batch")
ggsave("umap_after_harmony.png", plot = p2)

# CLUSTERING (already present)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined)

p3 <- DimPlot(combined, group.by = "seurat_clusters") + 
  ggtitle("After Harmony - Clusters")
ggsave("clusters_after_harmony.png", plot = p3)

library(cluster)

# Evaluation after Harmony
emb <- Embeddings(combined, "harmony")[, 1:30]

clusters <- as.numeric(as.factor(combined$seurat_clusters))
batch <- as.numeric(as.factor(combined$batch))

set.seed(42)
n_sub <- min(2000, nrow(emb))
idx <- sample(seq_len(nrow(emb)), n_sub)

dist_sub <- dist(emb[idx, ])

sil_clusters <- silhouette(clusters[idx], dist_sub)
cluster_asw <- mean(sil_clusters[, 3])

sil_batch <- silhouette(batch[idx], dist_sub)
batch_asw <- mean(sil_batch[, 3])

cat("Spatial Cluster ASW:", round(cluster_asw, 4), "\n")
cat("Spatial Batch ASW:", round(batch_asw, 4), "\n")

metrics <- data.frame(
  Method = "Harmony_spatial",
  Cluster_ASW = cluster_asw,
  Batch_ASW = batch_asw
)

write.csv(metrics, "spatial_harmony_metrics.csv", row.names = FALSE)

saveRDS(combined, "spatial_harmony_final.rds")

# FINAL METRICS: ARI, NMI, ASW BIO, BATCH ASW

library(cluster)
library(mclust)
library(aricode)

emb_harmony <- Embeddings(combined, "harmony")[, 1:30]

batch_labels <- combined$batch

# If you have true cell-type annotation, use:
# bio_labels <- combined$cell_type
# For now, using Harmony clusters as biological proxy
bio_labels <- combined$seurat_clusters

set.seed(42)
n_sub <- min(2000, nrow(emb_harmony))
idx <- sample(seq_len(nrow(emb_harmony)), n_sub)

emb_sub <- emb_harmony[idx, , drop = FALSE]
batch_sub <- as.factor(batch_labels[idx])
bio_sub <- as.factor(bio_labels[idx])

dist_sub <- dist(emb_sub)

# ASW bio / biological silhouette
sil_bio <- silhouette(as.numeric(bio_sub), dist_sub)
ASW_bio <- mean(sil_bio[, 3])

# Batch ASW
sil_batch <- silhouette(as.numeric(batch_sub), dist_sub)
ASW_batch <- mean(sil_batch[, 3])

# ARI: cluster vs batch
ARI_cluster_vs_batch <- adjustedRandIndex(bio_sub, batch_sub)

# NMI: cluster vs batch
NMI_cluster_vs_batch <- NMI(as.character(bio_sub), as.character(batch_sub))

metrics <- data.frame(
  Method = "Harmony_spatial",
  ASW_bio = ASW_bio,
  Batch_ASW = ASW_batch,
  ARI_cluster_vs_batch = ARI_cluster_vs_batch,
  NMI_cluster_vs_batch = NMI_cluster_vs_batch
)

write.csv(
  metrics,
  "spatial_harmony_final_metrics_ARI_NMI_ASW.csv",
  row.names = FALSE
)

cat("ASW Bio:", round(ASW_bio, 4), "\n")
cat("Batch ASW:", round(ASW_batch, 4), "\n")
cat("ARI cluster-vs-batch:", round(ARI_cluster_vs_batch, 4), "\n")
cat("NMI cluster-vs-batch:", round(NMI_cluster_vs_batch, 4), "\n")

cat("\nSaved: spatial_harmony_final_metrics_ARI_NMI_ASW.csv\n")

cat("\nInterpretation:\n")
cat("ASW Bio higher = better biological cluster separation\n")
cat("Batch ASW near 0 or low = better batch mixing\n")
cat("ARI cluster-vs-batch low = clusters are not batch-driven\n")
cat("NMI cluster-vs-batch low = clusters are not batch-driven\n")

cat("\n========== METRICS COMPLETE ==========\n")
