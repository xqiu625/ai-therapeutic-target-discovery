
.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

library(dplyr)
library(Seurat)
library(SeuratObject)
library(MAST)


dat0 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/harmony_int/302005_Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"

dat0@misc$seuratcluster_markers <- FindAllMarkers(object=dat0,
                                                   assay='RNA',only.pos=TRUE,test.use="MAST")

saveRDS(dat0@misc$seuratcluster_markers, '/bigdata/godziklab/shared/Xinru/302005/v5_output/Markers/302005_Platelets-Ery-Pre_cluster_marker.rds')

