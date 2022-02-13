.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

library(SeuratDisk)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(RISC)
library(SingleCellExperiment)
library(celldex)
library(SingleR)
library(HGNChelper)
library(tidyr)
library(stringr)


df302004data15 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets/302004data15/302004data15_covid19.rds")
SingleR <- read.delim2("/bigdata/godziklab/shared/Xinru/302005/datasets_v2/processed_datasets_v3/dat05_annotation_SingleR.tsv")
meta <- df302004data15@meta.data
meta$Cell_ID <- row.names(meta)
SingleR2 <- SingleR %>% dplyr::select(Cell_ID, pruned.labels)
meta2 <- merge(meta, SingleR2, by = "Cell_ID")
row.names(meta2) <- meta2$Cell_ID
df302004data15@meta.data <- meta2

df302004data15_Tcell <- subset(df302004data15, pruned.labels == 'T_cells')
df302004data15_Mono <- subset(df302004data15, pruned.labels == 'Monocyte')
df302004data15_Bcell <- subset(df302004data15, pruned.labels == 'B_cell')

# Integration with corrected gene names ----------------------------------------------------
setwd("/bigdata/godziklab/shared/Xinru/302005/datasets_v2/processed_datasets_v2/correct_gene_name")
file <- list.files(pattern = "*.txt")
gene_name <- do.call("rbind",lapply(file,
                                    FUN=function(file){read.delim2(file)}))
gene_name <- gene_name %>% distinct(x, .keep_all = T)

# Process seurat objects ready for Seurat v4 map--------------------------------------------------
reference_pbmc <- LoadH5Seurat("/bigdata/godziklab/shared/Xinru/v2_013020_sepsis/seurat_v4/pbmc_multimodal.h5seurat")
process0 <- function(obj){
  meta <- obj@meta.data
  mat <- obj@assays$RNA@counts
  dat <- as.data.frame(row.names(mat))
  names(dat) <- c("x")
  dat<- left_join(dat, gene_name, by = "x")
  dat <- dat %>% filter(!is.na(Suggested.Symbol))
  mat <- mat[dat$x, ]
  row.names(mat) <- dat$Suggested.Symbol
  obj <- CreateSeuratObject(counts = mat)
  obj@meta.data <- meta
  obj <- SCTransform(obj, verbose = FALSE)
  anchors <- FindTransferAnchors(
    reference = reference_pbmc,
    query = obj,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50,
    recompute.residuals = FALSE
  )
  obj <- MapQuery(
    anchorset = anchors,
    query = obj,
    reference = reference_pbmc,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}

dat5_Tcell <- process0(df302004data15_Tcell)
write.table(dat5_Tcell@meta.data,
            file = "/bigdata/godziklab/shared/Xinru/302005/datasets_v2/processed_datasets_v3/dat05_Tcell_annotation_SeuratV4.tsv",
            sep="\t")

dat5_Mono <- process0(df302004data15_Mono)
write.table(dat5_Mono@meta.data,
            file = "/bigdata/godziklab/shared/Xinru/302005/datasets_v2/processed_datasets_v3/dat05_Mono_annotation_SeuratV4.tsv",
            sep="\t")

dat5_Bcell<- process0(df302004data15_Bcell)
write.table(dat5_Bcell@meta.data,
            file = "/bigdata/godziklab/shared/Xinru/302005/datasets_v2/processed_datasets_v3/dat05_Bcell_annotation_SeuratV4.tsv",
            sep="\t")