
# HARMONY PIPELINE - Seurat 3.2.0 + Harmony 1.2.4

library(Seurat)
library(harmony)
library(ggplot2)
library(cluster)
library(aricode)

set.seed(42)

data_dir <- "scRNA_data/sc_sample"
samples <- list.dirs(data_dir, recursive = FALSE)
cat("Found samples:\n")
print(samples)

seurat_list <- list()

# STEP 1: LOAD + QC + PREPROCESS
for (x in samples[1:5]) {
  sample_name <- basename(x)
  cat("\nLoading:", sample_name, "\n")
  
  data <- Read10X(data.dir = x)
  
  obj <- CreateSeuratObject(
    counts       = data,
    project      = sample_name,
    min.cells    = 3,
    min.features = 200
  )
  obj$batch <- sample_name
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  obj <- subset(obj, subset =
                  nFeature_RNA > 200 &
                  nFeature_RNA < 6000 &
                  nCount_RNA   > 500  &
                  percent.mt   < 20)
  
  if (ncol(obj) > 3000) {
    obj <- subset(obj, cells = sample(colnames(obj), 3000))
  }
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  
  seurat_list[[sample_name]] <- obj
  rm(data); gc()
}

# STEP 2: MERGE
# Seurat 3: no JoinLayers, no merge.data argument
combined <- merge(
  seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

combined <- FindVariableFeatures(combined, nfeatures = 2000, verbose = FALSE)
saveRDS(combined, "combined_merge_final.rds")

# STEP 3: SCALE + PCA
combined <- ScaleData(combined,
                      vars.to.regress = "percent.mt",
                      verbose = FALSE)

combined <- RunPCA(combined,
                   features = VariableFeatures(combined),
                   npcs     = 50,
                   verbose  = FALSE)

png("ElbowPlot.png", width = 900, height = 700)
print(ElbowPlot(combined, ndims = 50))
dev.off()

# STEP 4: BEFORE HARMONY
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20,
                    reduction.name = "umap_before", verbose = FALSE)

combined$cluster_before <- Idents(combined)

png("UMAP_before_batch.png", width = 1200, height = 1000, res = 150)
print(DimPlot(combined, reduction = "umap_before", group.by = "batch") +
        ggtitle("Before Harmony - Batch"))
dev.off()

png("UMAP_before_clusters.png", width = 1200, height = 1000, res = 150)
print(DimPlot(combined, reduction = "umap_before",
              group.by = "cluster_before", label = TRUE) +
        ggtitle("Before Harmony - Clusters"))
dev.off()

# STEP 5: RUN HARMONY
# Seurat 3 compatible — use HarmonyMatrix directly
pca_embeddings <- Embeddings(combined, reduction = "pca")

harmony_embeddings <- HarmonyMatrix(
  data_mat  = pca_embeddings,
  meta_data = combined@meta.data,
  vars_use  = "batch",
  do_pca    = FALSE,
  verbose   = TRUE
)

# Store harmony as a new reduction manually
combined[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embeddings,
  key        = "harmony_",
  assay      = DefaultAssay(combined)
)

cat("Harmony done. Reductions:", names(combined@reductions), "\n")

# STEP 6: AFTER HARMONY
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20,
                    reduction.name = "umap_harmony", verbose = FALSE)

combined$cluster_after <- Idents(combined)

png("UMAP_after_batch.png", width = 1200, height = 1000, res = 150)
print(DimPlot(combined, reduction = "umap_harmony", group.by = "batch") +
        ggtitle("After Harmony - Batch"))
dev.off()

png("UMAP_after_clusters.png", width = 1200, height = 1000, res = 150)
print(DimPlot(combined, reduction = "umap_harmony",
              group.by = "cluster_after", label = TRUE) +
        ggtitle("After Harmony - Clusters"))
dev.off()

saveRDS(combined, file = "combined_harmony_final.rds")

# STEP 7: SingleR ANNOTATION
# Seurat 3: use GetAssayData, no layer argument
cat("\nRunning SingleR annotation...\n")

library(SingleR)

ref <- readRDS("data/HumanPrimaryCellAtlasData.rds")

# Seurat 3 — no layer argument
counts_mat <- as.matrix(GetAssayData(combined, assay = "RNA", slot = "data"))

singler_results <- SingleR::SingleR(
  test      = counts_mat,
  ref       = ref,
  labels    = ref$label.main,
  fine.tune = TRUE
)

combined$singler_labels <- singler_results$labels
combined$singler_scores <- apply(singler_results$scores, 1, max)

cat("Cell types found:\n")
print(table(combined$singler_labels))

png("UMAP_singler_labels.png", width = 1400, height = 1000, res = 150)
print(DimPlot(combined,
              reduction  = "umap_harmony",
              group.by   = "singler_labels",
              label      = TRUE,
              repel      = TRUE,
              label.size = 3) +
        ggtitle("SingleR Cell Type Annotation"))
dev.off()

png("SingleR_score_heatmap.png", width = 1400, height = 900, res = 150)
SingleR::plotScoreHeatmap(singler_results)
dev.off()

write.csv(
  as.data.frame(table(
    Cluster  = combined$cluster_after,
    CellType = combined$singler_labels
  )),
  "cluster_celltype_table.csv",
  row.names = FALSE
)

saveRDS(combined, file = "combined_harmony_final.rds")
cat("SingleR done.\n")

# STEP 8: EVALUATION METRICS
cat("\nComputing metrics...\n")

emb_pca     <- Embeddings(combined, "pca")[, 1:20]
emb_harmony <- Embeddings(combined, "harmony")[, 1:20]

set.seed(42)
n_sub <- min(2000, nrow(emb_harmony))
idx   <- sample(seq_len(nrow(emb_harmony)), n_sub)

emb_pca_sub  <- emb_pca[idx, ]
emb_harm_sub <- emb_harmony[idx, ]
cl_before    <- as.character(combined$cluster_before)[idx]
cl_after     <- as.character(combined$cluster_after)[idx]
batch_sub    <- as.character(combined$batch)[idx]
celltype_sub <- as.character(combined$singler_labels)[idx]

dist_pca  <- dist(emb_pca_sub)
dist_harm <- dist(emb_harm_sub)

# ASW cluster
asw_cl_before <- round(mean(silhouette(as.numeric(factor(cl_before)), dist_pca)[, 3]),  4)
asw_cl_after  <- round(mean(silhouette(as.numeric(factor(cl_after)),  dist_harm)[, 3]), 4)

# ASW batch
asw_batch_before <- round(mean(silhouette(as.numeric(factor(batch_sub)), dist_pca)[, 3]),  4)
asw_batch_after  <- round(mean(silhouette(as.numeric(factor(batch_sub)), dist_harm)[, 3]), 4)

# NMI
nmi_before <- round(NMI(cl_before, celltype_sub), 4)
nmi_after  <- round(NMI(cl_after,  celltype_sub), 4)

# ARI
ari_before <- round(ARI(cl_before, celltype_sub), 4)
ari_after  <- round(ARI(cl_after,  celltype_sub), 4)

# iLISI
compute_ilisi <- function(emb, batch_labels, k = 30) {
  n <- nrow(emb)
  scores <- numeric(n)
  d <- as.matrix(dist(emb))
  batch_labels <- as.factor(batch_labels)
  for (i in seq_len(n)) {
    neighbors <- order(d[i, ])[2:(k + 1)]
    freqs <- table(batch_labels[neighbors]) / k
    scores[i] <- 1 / sum(freqs^2)
  }
  return(scores)
}

# cLISI
compute_clisi <- function(emb, celltype_labels, k = 30) {
  n <- nrow(emb)
  scores <- numeric(n)
  d <- as.matrix(dist(emb))
  celltype_labels <- as.factor(celltype_labels)
  for (i in seq_len(n)) {
    neighbors <- order(d[i, ])[2:(k + 1)]
    freqs <- table(celltype_labels[neighbors]) / k
    scores[i] <- 1 / sum(freqs^2)
  }
  return(scores)
}

cat("Computing iLISI and cLISI...\n")

set.seed(42)
lisi_idx <- sample(seq_len(nrow(emb_harm_sub)), min(500, nrow(emb_harm_sub)))

ilisi_after  <- round(mean(compute_ilisi(emb_harm_sub[lisi_idx, ], batch_sub[lisi_idx])),    4)
clisi_after  <- round(mean(compute_clisi(emb_harm_sub[lisi_idx, ], celltype_sub[lisi_idx])), 4)
ilisi_before <- round(mean(compute_ilisi(emb_pca_sub[lisi_idx, ],  batch_sub[lisi_idx])),    4)
clisi_before <- round(mean(compute_clisi(emb_pca_sub[lisi_idx, ],  celltype_sub[lisi_idx])), 4)

# RESULTS TABLE
results_table <- data.frame(
  Metric         = c("ASW_cluster", "ASW_batch", "NMI", "ARI", "iLISI", "cLISI"),
  Before_Harmony = c(asw_cl_before, asw_batch_before, nmi_before,
                     ari_before, ilisi_before, clisi_before),
  After_Harmony  = c(asw_cl_after,  asw_batch_after,  nmi_after,
                     ari_after,  ilisi_after,  clisi_after),
  Direction      = c("higher=better", "lower=better", "higher=better",
                     "higher=better", "higher=better", "lower=better")
)

cat("\n========== HARMONY METRICS ==========\n")
print(results_table)

write.csv(results_table, "harmony_metrics_scrna.csv", row.names = FALSE)
cat("\nDone. Saved harmony_metrics_scrna.csv\n")
