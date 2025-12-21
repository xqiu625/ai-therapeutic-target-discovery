#!/usr/bin/env Rscript
# ============================================================
# Extract Seurat Data for Python/Scanpy
# Compatible with Seurat v4
# ============================================================

library(Seurat)
library(Matrix)

# ============================================================
# Configuration
# ============================================================
RDS_FILE <- "/bigdata/godziklab/shared/Xinru/302004/sepsis_singlecell_transcriptome.rds"
OUTPUT_DIR <- "/bigdata/godziklab/shared/Xinru/302004/302004_git/data/processed"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load Seurat object
# ============================================================
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(RDS_FILE)
print(seurat_obj)

cat("\nMetadata columns:\n")
print(colnames(seurat_obj@meta.data))

# ============================================================
# 1. Extract count matrix (Seurat v4 syntax)
# ============================================================
cat("\n============================================================\n")
cat("EXTRACTING COUNT MATRIX\n")
cat("============================================================\n")

# Use GetAssayData with slot (Seurat v4)
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
cat("Count matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")

# Save count matrix
counts_file <- file.path(OUTPUT_DIR, "counts.mtx")
cat("Writing counts to:", counts_file, "\n")
writeMM(counts, counts_file)
system(paste("gzip -f", counts_file))
cat("Saved counts.mtx.gz\n")

# ============================================================
# 2. Extract gene names
# ============================================================
cat("\n============================================================\n")
cat("EXTRACTING GENE NAMES\n")
cat("============================================================\n")

genes <- data.frame(gene_name = rownames(counts), stringsAsFactors = FALSE)
genes_file <- file.path(OUTPUT_DIR, "genes.csv")
write.csv(genes, genes_file, row.names = FALSE)
cat("Saved", nrow(genes), "genes to genes.csv\n")

# ============================================================
# 3. Extract cell metadata
# ============================================================
cat("\n============================================================\n")
cat("EXTRACTING CELL METADATA\n")
cat("============================================================\n")

metadata <- seurat_obj@meta.data
metadata$barcode <- rownames(metadata)
metadata_file <- file.path(OUTPUT_DIR, "metadata.csv")
write.csv(metadata, metadata_file, row.names = FALSE)
cat("Saved metadata for", nrow(metadata), "cells\n")
cat("Columns:", paste(colnames(metadata), collapse = ", "), "\n")

# ============================================================
# 4. Extract embeddings
# ============================================================
cat("\n============================================================\n")
cat("EXTRACTING EMBEDDINGS\n")
cat("============================================================\n")

# PCA
if ("pca" %in% names(seurat_obj@reductions)) {
    pca <- Embeddings(seurat_obj, "pca")
    write.csv(pca, file.path(OUTPUT_DIR, "pca.csv"))
    cat("Saved PCA:", ncol(pca), "components\n")
} else {
    cat("No PCA found\n")
}

# UMAP
if ("umap" %in% names(seurat_obj@reductions)) {
    umap_emb <- Embeddings(seurat_obj, "umap")
    write.csv(umap_emb, file.path(OUTPUT_DIR, "umap.csv"))
    cat("Saved UMAP\n")
} else {
    cat("No UMAP found\n")
}

# tSNE
if ("tsne" %in% names(seurat_obj@reductions)) {
    tsne <- Embeddings(seurat_obj, "tsne")
    write.csv(tsne, file.path(OUTPUT_DIR, "tsne.csv"))
    cat("Saved tSNE\n")
} else {
    cat("No tSNE found\n")
}

# ============================================================
# Summary
# ============================================================
cat("\n============================================================\n")
cat("EXTRACTION COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Files created:\n")
system(paste("ls -lh", OUTPUT_DIR))
