.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )


library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)

genesets <- getGmt("/bigdata/godziklab/shared/Xinru/302005/v5_output/GSVA/h.all.v7.5.1.symbols.gmt")

dat0 <- readRDS("/bigdata/godziklab/shared/Xinru/302005/datasets_v5/harmony_int/302005_Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"
dat0@meta.data$seurat_clusters <- paste0("C", dat0@meta.data$seurat_clusters)

df.data <- GetAssayData(object = dat0, slot = "data")
df.group <- data.frame(umi = names(Idents(dat0)), 
                       cluster = as.character(dat0@meta.data$seurat_clusters), 
                       stringsAsFactors = F)

gsvascore <- gsva(data.matrix(df.data), genesets, parallel.sz = 2)
saveRDS(gsvascore, '/bigdata/godziklab/shared/Xinru/302005/v5_output/GSVA/gsvascore.RDS')

# 先用热图展示一下hallmarker中全部的50个geneset得到的结果
ha.t <- HeatmapAnnotation(Cluster = df.group$cluster)
dpi = 300
png(file = '/bigdata/godziklab/shared/Xinru/302005/v5_output/GSVA/seurat_cluster_gsvascore.png',
    width = dpi * 12,height = dpi * 6,units = "px",res = dpi,type = 'cairo')

Heatmap(as.matrix(gsvascore), 
        show_column_names = F, 
        cluster_rows = T, 
        cluster_columns = T, 
        top_annotation = ha.t, 
        column_split = df.group$cluster, 
        row_names_gp = gpar(fontsize = 8), 
        row_names_max_width = max_text_width(rownames(gsvascore), 
                                             gp = gpar(fontsize = 8)))

dev.off()


