
.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
library(purrr)

setwd('/bigdata/godziklab/shared/Xinru/302005/datasets_v5/B_cell')
file <- list.files(pattern = "*.rds")
datls <- lapply(file, FUN=function(file){readRDS(file)})

### original mat
mats <- lapply(datls, FUN=function(datls){datls@assays$RNA@counts})

### intersect genes cropped mat
symbols <- lapply(mats, FUN=function(mats){row.names(mats)})
var0 = Reduce(intersect, symbols)
mats2 <- lapply(mats, FUN=function(mats){mats[var0,]})

### metadata 
metas <- do.call("rbind",lapply(datls, FUN=function(datls){datls@meta.data %>%
    select(Cell_ID,Sample_ID,Data_NO)}))


pbmc <- CreateSeuratObject(counts = cbind(mats2[[1]], 
                                          mats2[[2]],
                                          mats2[[3]],
                                          mats2[[4]],
                                          mats2[[5]],
                                          mats2[[6]],
                                          mats2[[7]],
                                          mats2[[8]],
                                          mats2[[9]],
                                          mats2[[10]],
                                          mats2[[11]],
                                          mats2[[12]]
                                          )) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = pbmc@var.genes, npcs = 30, verbose = FALSE)

metas2 <- cbind(pbmc@meta.data,metas)
pbmc@meta.data <- metas2

dpi = 300
png(file = "/bigdata/godziklab/shared/Xinru/302005/v5_output/Harmony/Harmony_convergence_Data_NO_B_cell.png", 
    width = dpi * 4,height = dpi * 4,
    units = "px",res = dpi,type = 'cairo')
pbmc <- pbmc %>% 
  RunHarmony("Data_NO", plot_convergence = TRUE)
dev.off()

harmony_embeddings <- Embeddings(pbmc, 'harmony')
#harmony_embeddings[1:5, 1:5]

saveRDS(pbmc, '/bigdata/godziklab/shared/Xinru/302005/datasets_v5/harmony_int/302005_B_cell.rds')

# Downstream analysis -----------------------------------------------------
pbmc <- pbmc %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

saveRDS(pbmc, '/bigdata/godziklab/shared/Xinru/302005/datasets_v5/harmony_int/302005_B_cell_cluster.rds')



